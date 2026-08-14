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
    score/verdict. Does NOT filter the set: downstream protocols decide
    what to do with the annotation (``epitopeconstruct.utils.assembly.
    select_bcell_candidates`` keeps 'Allergen' candidates in the construct,
    with the verdict left visible on the ROI for an informed decision.
    HTL/CTL candidates are never filtered by this verdict to begin with,
    since an MHC-buried core is never a free circulating epitope for IgE
    recognition).

    Model: runs AlgPred2 in hybrid mode (RF + BLAST + MERCI) by default,
    not the ML-only mode. The ML-only mode classifies purely from
    amino-acid composition, with no comparison against real IgE allergen
    sequences or motifs, and is prone to marking 'Allergen' candidates with
    no such evidence behind them. The evidence breakdown
    (``_algpredMlScore``/``_algpredMerciScore``/``_algpredBlastScore``) is
    always exposed so a verdict backed by real homology/motif evidence can
    be told apart from one that is composition statistics alone (the last
    two are always 0 if run in ML-only mode).

    This same protocol is meant to be reused twice in a workflow: once on
    a large SetOfSequenceROIs of per-peptide B-cell candidates (upstream of
    construct assembly), and once on a single-ROI SetOfSequenceROIs wrapping
    the final assembled multi-epitope construct (construct-level check).
    The batch-size-1 AlgPred2 bug (see utils/algpred.py) makes the second
    case the NORMAL code path for that use, not a rare edge case.

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, with each ROI
    annotated with ``_algpredScore`` (float, Hybrid Score in the default
    model), ``_algpredVerdict`` (``'Allergen'``/``'Non-Allergen'``, raw
    AlgPred2 text), and the ``_algpredMlScore``/``_algpredMerciScore``/
    ``_algpredBlastScore`` evidence breakdown (float each).
    """

    _label = 'algpred2 allergenicity'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates to evaluate for allergenicity.')
        form.addParam('model', params.EnumParam, choices=['Hybrid (RF+BLAST+MERCI)', 'ML only (RF)'],
                       default=0, label='Model: ',
                       help='Hybrid (default): combines the RF classifier with BLAST against a real '
                            'IgE allergen database and MERCI against documented IgE motifs (Sharma et '
                            'al. 2021, PMID 33201237). ML only: RandomForest over amino-acid '
                            'composition alone, no real homology/motif evidence -- prone to marking '
                            "'Allergen' candidates with no clinical correlate.")
        form.addParam('threshold', params.FloatParam, default=0.3,
                       label='Score threshold: ',
                       help='AlgPred2 Hybrid Score (or ML_Score, in ML-only model) threshold used '
                            'internally by the tool to classify Allergen vs. Non-Allergen.')
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

    def _getModelArg(self):
        # EnumParam index 0 ('Hybrid') -> AlgPred2 -m 2; index 1 ('ML only') -> -m 1.
        return 1 if self.model.get() == 1 else 2

    def algpredStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        fastaPath = self._getExtraPath('candidates.fasta')
        write_fasta(sequences, fastaPath)

        args = (f'-i {fastaPath} -o {self._getRawOutputPath()} '
                f'-t {self.threshold.get()} -m {self._getModelArg()} -d 2')
        algpredPlugin.runAlgPred2(self, args)

    def createOutputStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        padded = len(sequences) == 1
        resultDf = parse_output(self._getRawOutputPath(), sequences=sequences, padded=padded)

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi, row in zip(rois, resultDf.itertuples(index=False)):
            roi._algpredScore = Float(row.algpred_score)
            roi._algpredVerdict = String(row.algpred_verdict)
            roi._algpredMlScore = Float(row.algpred_ml_score)
            roi._algpredMerciScore = Float(row.algpred_merci_score)
            roi._algpredBlastScore = Float(row.algpred_blast_score)
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
