from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence
from pwchem.protocols import ProtDefineSeqROI

from ..protocols import ProtAlgPred2Prediction


class TestAlgPred2Prediction(BaseTest):
    NAME = 'ALGPRED_TEST_SEQ'
    DESCRIPTION = 'Synthetic concatenation of 4 real reference peptides'
    # 4 real peptides (2 known HLA-A*02:01 CD8+ reference epitopes, 2 real
    # HIV-1 Env GP120 windows also used in the netmhcpan test fixture),
    # concatenated with 'GGG' spacers purely to build one importable
    # sequence -- the AlgPred2 evaluation itself is per-ROI/per-peptide, not
    # position-dependent.
    PEPTIDES = ['NLVPMVATV', 'GILGFVFTL', 'RAIEAQQHL', 'NAKTIIVQL']
    SPACER = 'GGG'
    AMINOACIDSSEQ = SPACER.join(PEPTIDES)

    # Real AlgPred2 output in hybrid mode (score threshold 0.3, default),
    # from a direct local run of the real binary -- not estimated. 2
    # Allergen + 2 Non-Allergen, a genuine mixed result confirming both
    # verdict paths. None of the 4 peptides hit the BLAST allergen
    # database or a MERCI IgE motif, so Hybrid Score equals the RF
    # (ML) score alone here -- still exercises the hybrid-mode output
    # schema end to end.
    EXPECTED_VERDICT = {
        'NLVPMVATV': ('Allergen', 0.335),
        'GILGFVFTL': ('Non-Allergen', 0.286),
        'RAIEAQQHL': ('Allergen', 0.349),
        'NAKTIIVQL': ('Non-Allergen', 0.29),
    }

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls._runImportSeq()
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

        cls.protSeedROIs = cls._runDefSeqROIs(cls.protImportSeq)
        cls._waitOutput(cls.protSeedROIs, 'outputROIs', sleepTime=5)

    @classmethod
    def _runImportSeq(cls):
        kwargs = {
            'inputSequenceName': cls.NAME,
            'inputSequenceDescription': cls.DESCRIPTION,
            'inputRawSequence': cls.AMINOACIDSSEQ,
        }
        cls.protImportSeq = cls.newProtocol(ProtImportSequence, **kwargs)
        cls.proj.launchProtocol(cls.protImportSeq, wait=False)

    @classmethod
    def _getWindows(cls):
        windows = []
        cursor = 0
        for pep in cls.PEPTIDES:
            start = cls.AMINOACIDSSEQ.index(pep, cursor) + 1
            end = start + len(pep) - 1
            windows.append((start, end))
            cursor = end
        return windows

    @classmethod
    def _runDefSeqROIs(cls, inProt):
        windows = cls._getWindows()
        inROIs = '\n'.join(
            '{}) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}'.format(
                i, start, end, cls.AMINOACIDSSEQ[start - 1:end]
            )
            for i, (start, end) in enumerate(windows, 1)
        )
        protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
        protDefSeqROIs.inputSequence.set(inProt)
        protDefSeqROIs.inputSequence.setExtended('outputSequence')

        cls.proj.launchProtocol(protDefSeqROIs, wait=False)
        return protDefSeqROIs

    def test(self):
        protAlgPred = self.newProtocol(ProtAlgPred2Prediction)
        protAlgPred.inputROIs.set(self.protSeedROIs)
        protAlgPred.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protAlgPred, wait=True)

        outROIs = getattr(protAlgPred, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), len(self.PEPTIDES))

        for roi in outROIs:
            seq = roi.getROISequence()
            expectedVerdict, expectedScore = self.EXPECTED_VERDICT[seq]
            self.assertEqual(roi._algpredVerdict.get(), expectedVerdict)
            self.assertAlmostEqual(roi._algpredScore.get(), expectedScore, places=3)
            self.assertAlmostEqual(roi._algpredMlScore.get(), expectedScore, places=3)
            self.assertAlmostEqual(roi._algpredMerciScore.get(), 0.0, places=3)
            self.assertAlmostEqual(roi._algpredBlastScore.get(), 0.0, places=3)
