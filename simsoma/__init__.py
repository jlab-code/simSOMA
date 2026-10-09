"""simSOMA: cell-lineage simulator of somatic VAF spectra in plants.

The simulator modules live in the ``simsoma_corefunc`` package (repository folder
``simSOMA_corefunc/``); they use flat imports, so ``simsoma.core_path()`` is put on sys.path
by the CLI. The shared sequencing observation model is the separate ``plantsoma_obs`` package.
"""
from pathlib import Path

__version__ = "0.2.0"


def core_path() -> Path:
    """Directory containing the simulator modules (installed or repository checkout)."""
    try:
        import simsoma_corefunc
        return Path(simsoma_corefunc.__file__).resolve().parent
    except ImportError:
        return Path(__file__).resolve().parents[1] / "simSOMA_corefunc"


def use_core() -> Path:
    """Make the flat simulator modules importable (``import pipeline_wrapper`` etc.)."""
    import sys
    p = str(core_path())
    if p not in sys.path:
        sys.path.insert(0, p)
    return Path(p)
