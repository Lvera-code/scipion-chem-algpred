=========
CHANGES
=========

0.2.0
=====
- ``ProtAlgPred2Prediction`` now runs AlgPred2 in hybrid mode (RF + BLAST +
  MERCI) by default instead of ML-only, exposed as a new ``model`` param.
  The ML-only mode classifies from amino-acid composition alone, with no
  comparison against real IgE allergen sequences or motifs, and is prone to
  false 'Allergen' verdicts with no such evidence behind them. The output
  now also exposes the ``_algpredMlScore``/``_algpredMerciScore``/
  ``_algpredBlastScore`` breakdown on every ROI, so a verdict backed by
  real homology/motif evidence can be told apart from one that is
  composition statistics alone.

0.1.0
=====
- Initial release: AlgPred2 allergenicity annotation protocol
  (``ProtAlgPred2Prediction``), auto-installed via conda/pip (no academic
  license gate, unlike the DTU tools in this project). Includes the
  real batch-size-1 sklearn reshape workaround, needed for both the
  per-peptide (Fase 4b) and construct-level (Fase 8) reuse of this
  protocol.
