================================
AlgPred2 Scipion plugin
================================

Scipion framework plugin wrapping AlgPred2 (Raghava group, GPL-3.0) for
allergenicity prediction of peptide candidates. Unlike BepiPred/NetMHCpan/
NetMHCIIpan/SignalP, AlgPred2 is a plain pip package with no academic
license restriction: it is installed automatically in its own conda
environment.

The plugin implements a single protocol, ``ProtAlgPred2Prediction``, which
takes an arbitrary ``SetOfSequenceROIs`` and annotates every ROI with its
AlgPred2 score/verdict (``_algpredScore``/``_algpredVerdict``). It does
**not** filter the set -- downstream protocols decide what to do with the
annotation.

This same protocol is meant to be reused twice in a workflow: once on a
large set of per-peptide B-cell candidates, and once on a single-ROI set
wrapping a final assembled multi-epitope construct (construct-level
allergenicity check). AlgPred2's own classifier has a real reshape bug for
single-sequence input, which the plugin works around transparently (see
``utils/algpred.py``) -- this makes the construct-level use the *normal*
code path for that bug, not a rare edge case.

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-algpred.git
            cd scipion-chem-algpred
            scipion3 installp -p . --devel
