"""Parsing of AlgPred2 output, including the real batch-size-1 workaround.

AlgPred2's own sklearn RandomForest classifier throws
"ValueError: Expected 2D array, got 1D array" when the input FASTA has
exactly 1 sequence (a real upstream reshape bug -- confirmed by running the
real binary, not assumed). This matters for this plugin specifically
because the construct-level check (this same protocol reused on the final
assembled multi-epitope construct) always submits exactly 1 sequence, so
this is the NORMAL code path there, not a rare edge case.
"""

from pathlib import Path
from typing import List

import pandas as pd

_RAW_REQUIRED_COLUMNS = {"Sequence", "ML_Score", "Prediction"}


class AlgPred2ParseError(Exception):
    """The AlgPred2 output CSV does not match the expected format."""


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


def parse_output(csv_path: Path, n_expected: int, padded: bool) -> pd.DataFrame:
    """Parse AlgPred2's raw output CSV.

    Args:
        csv_path: Path to AlgPred2's raw output CSV.
        n_expected: Number of real (non-padded) input sequences.
        padded: Whether ``write_fasta`` duplicated a single sequence -- if
            so, the extra row is dropped before returning.

    Returns:
        DataFrame with columns ``sequence``, ``algpred_score``,
        ``algpred_verdict`` (raw AlgPred2 text: ``'Allergen'``/
        ``'Non-Allergen'``), one row per real input sequence, in the same
        order they were submitted.
    """
    try:
        raw = pd.read_csv(csv_path)
    except Exception as exc:
        raise AlgPred2ParseError(f"Could not parse AlgPred2 output at '{csv_path}': {exc}") from exc

    if not _RAW_REQUIRED_COLUMNS.issubset(raw.columns):
        raise AlgPred2ParseError(
            f"AlgPred2 output CSV format does not match what was expected: missing columns "
            f"{_RAW_REQUIRED_COLUMNS - set(raw.columns)}. Columns found: {list(raw.columns)}."
        )

    if padded:
        raw = raw.iloc[:n_expected]

    return pd.DataFrame({
        'sequence': raw['Sequence'].to_numpy(),
        'algpred_score': raw['ML_Score'].to_numpy(),
        'algpred_verdict': raw['Prediction'].to_numpy(),
    })
