"""Parsing of AlgPred2 output, including the real batch-size-1 workaround
and the two output schemas the underlying script can produce depending on
which model it ran.

AlgPred2's own sklearn RandomForest classifier throws
"ValueError: Expected 2D array, got 1D array" when the input FASTA has
exactly 1 sequence (a real upstream reshape bug -- confirmed by running the
real binary, not assumed). This matters for this plugin specifically
because the construct-level check (this same protocol reused on the final
assembled multi-epitope construct) always submits exactly 1 sequence, so
this is the NORMAL code path there, not a rare edge case.

Model (``-m``, hybrid by default): the underlying script supports two
modes with DIFFERENT csv output schemas, so each is parsed separately
instead of assuming one shape. The hybrid model is preferred because the
ML-only model classifies from amino-acid composition alone, with no
comparison against real IgE allergen sequences or motifs, and is prone to
marking 'Allergen' candidates with no such evidence behind them.

    Model 1 (ML only): columns ``ID,Sequence,ML_Score,Prediction`` --
        RandomForest over amino-acid composition (AAC) only, no comparison
        against real IgE allergens/motifs.
    Model 2 (hybrid, default): columns
        ``Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction``
        -- combines the same RF with BLAST against a real IgE allergen
        database and MERCI against documented IgE motifs (Sharma et al.
        2021, PMID 33201237). Has NO ``Sequence`` column: each row's
        sequence is reconstructed by POSITION against the submitted list
        (verified the order is preserved end-to-end by the upstream script,
        same ``seqid``/``seq`` lists reused at every internal step).
"""

from pathlib import Path
from typing import List

import pandas as pd

_MODEL1_COLUMNS = {"Sequence", "ML_Score", "Prediction"}
_MODEL2_COLUMNS = {"Subject", "ML Score", "MERCI Score", "BLAST Score", "Hybrid Score", "Prediction"}
_OUTPUT_COLUMNS = [
    'sequence', 'algpred_score', 'algpred_verdict',
    'algpred_ml_score', 'algpred_merci_score', 'algpred_blast_score',
]


class AlgPred2ParseError(Exception):
    """The AlgPred2 output CSV does not match either expected format."""


def write_fasta(sequences: List[str], fasta_path: Path) -> bool:
    """Write ``sequences`` to ``fasta_path``, duplicating a single sequence.

    Returns:
        True if the sequence list was padded (single-sequence workaround),
        so the caller knows to drop the extra row from the parsed result.
    """
    padded = len(sequences) == 1
    toWrite = sequences * 2 if padded else sequences
    with open(fasta_path, 'w') as fh:
        for i, seq in enumerate(toWrite):
            fh.write(f'>candidate_{i}\n{seq}\n')
    return padded


def parse_output(csv_path: Path, sequences: List[str], padded: bool) -> pd.DataFrame:
    """Parse AlgPred2's raw output CSV, whichever model produced it.

    Args:
        csv_path: Path to AlgPred2's raw output CSV.
        sequences: The real (non-padded) input sequences, in the order they
            were submitted.
        padded: Whether ``write_fasta`` duplicated a single sequence -- if
            so, the extra row is dropped before returning.

    Returns:
        DataFrame with columns ``sequence``, ``algpred_score`` (Hybrid Score
        in model 2, ML_Score in model 1), ``algpred_verdict`` (raw AlgPred2
        text: ``'Allergen'``/``'Non-Allergen'``), and the
        ``algpred_ml_score``/``algpred_merci_score``/``algpred_blast_score``
        breakdown (the last two always 0 in model 1, which does not compute
        them), one row per real input sequence, in submission order.

    Raises:
        AlgPred2ParseError: If the CSV matches neither known schema, or (for
            the hybrid model, which lacks a ``Sequence`` column) its row
            count does not match the submitted sequences, making a safe
            positional reconstruction impossible.
    """
    try:
        raw = pd.read_csv(csv_path)
    except Exception as exc:
        raise AlgPred2ParseError(f"Could not parse AlgPred2 output at '{csv_path}': {exc}") from exc

    cols = set(raw.columns)
    written = sequences * 2 if padded else sequences

    if _MODEL2_COLUMNS.issubset(cols):
        if len(raw) != len(written):
            raise AlgPred2ParseError(
                f"AlgPred2 (hybrid model) returned {len(raw)} row(s) for {len(written)} submitted "
                "sequence(s) -- cannot safely reconstruct sequence by position."
            )
        result = pd.DataFrame({
            'sequence': written,
            'algpred_score': raw['Hybrid Score'].to_numpy(),
            'algpred_verdict': raw['Prediction'].to_numpy(),
            'algpred_ml_score': raw['ML Score'].to_numpy(),
            'algpred_merci_score': raw['MERCI Score'].to_numpy(),
            'algpred_blast_score': raw['BLAST Score'].to_numpy(),
        })
    elif _MODEL1_COLUMNS.issubset(cols):
        result = pd.DataFrame({
            'sequence': raw['Sequence'].to_numpy(),
            'algpred_score': raw['ML_Score'].to_numpy(),
            'algpred_verdict': raw['Prediction'].to_numpy(),
            'algpred_ml_score': raw['ML_Score'].to_numpy(),
            'algpred_merci_score': 0.0,
            'algpred_blast_score': 0.0,
        })
    else:
        raise AlgPred2ParseError(
            "AlgPred2 output CSV format does not match either expected schema. "
            f"Columns found: {list(raw.columns)}."
        )

    if padded:
        result = result.iloc[:len(sequences)]

    return result[_OUTPUT_COLUMNS].reset_index(drop=True)
