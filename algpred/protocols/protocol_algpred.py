# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Enzo Sierra (enzogael57@gmail.com)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This protocol is used to predict allergenicity of a set of peptide
candidates with a local AlgPred2 installation.
"""

from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, String
from pyworkflow.protocol import params

from .. import Plugin as algpredPlugin
from ..utils.algpred import parse_output, write_fasta


class ProtAlgPred2Prediction(EMProtocol):
    """
    AI Generated:

    Predicts allergenicity of a set of peptide candidates using a local
    AlgPred2 installation, and annotates every input ROI with the resulting
    score/verdict (does NOT filter the set: downstream protocols decide
    what to do with the annotation. DECISION 2026-08-13: the multi-epitope
    construct assembly protocol used to exclude 'Allergen' B-cell
    candidates -- it no longer does, see
    ``epitopeconstruct.utils.assembly.select_bcell_candidates``, ported
    from the same-day decision in the standalone project's publication
    validation panel. HTL/CTL candidates were never filtered by this verdict
    to begin with, since an MHC-buried core is never a free circulating
    epitope for IgE recognition).

    This same protocol is meant to be reused twice in a workflow: once on
    a large SetOfSequenceROIs of per-peptide B-cell candidates (upstream of
    construct assembly), and once on a single-ROI SetOfSequenceROIs wrapping
    the final assembled multi-epitope construct (construct-level check).
    The batch-size-1 AlgPred2 bug (see utils/algpred.py) makes the second
    case the NORMAL code path for that use, not a rare edge case.

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, with each ROI
    annotated with ``_algpredScore`` (float) and ``_algpredVerdict``
    (``'Allergen'``/``'Non-Allergen'``, raw AlgPred2 text).
    """

    _label = 'algpred2 allergenicity'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates to evaluate for allergenicity.')
        form.addParam('threshold', params.FloatParam, default=0.3,
                       label='ML_Score threshold: ',
                       help='AlgPred2 ML_Score threshold used internally by the tool to classify '
                            'Allergen vs. Non-Allergen.')
        form.addParam('timeoutSeconds', params.IntParam, label='Timeout (s): ', default=300,
                       expertLevel=params.LEVEL_ADVANCED)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.algpredStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRawOutputPath(self):
        return self._getExtraPath('algpred_raw_output.csv')

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def algpredStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        fastaPath = self._getExtraPath('candidates.fasta')
        write_fasta(sequences, fastaPath)

        args = f'-i {fastaPath} -o {self._getRawOutputPath()} -t {self.threshold.get()} -d 2'
        algpredPlugin.runAlgPred2(self, args)

    def createOutputStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        padded = len(sequences) == 1
        resultDf = parse_output(self._getRawOutputPath(), n_expected=len(sequences), padded=padded)

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi, row in zip(rois, resultDf.itertuples(index=False)):
            roi._algpredScore = Float(row.algpred_score)
            roi._algpredVerdict = String(row.algpred_verdict)
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return algpredPlugin.validateInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nAllergen = sum(1 for roi in outROIs if roi._algpredVerdict.get() == 'Allergen')
                summary.append(f'{nAllergen}/{len(outROIs)} candidate(s) classified as Allergen.')
        return summary
