"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/CalledProcessError escape to the Scipion GUI without an
actionable message.
"""


class AlgPred2ExecutionError(Exception):
    """Failed to run AlgPred2 locally: missing installation, failed/timed-out
    subprocess, or the output CSV was not generated / does not match the
    expected format."""
