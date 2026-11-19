"""simSOMA: cell-lineage simulator of somatic VAF spectra in plants.

The simulator modules are submodules of this package (``from simsoma import pipeline_wrapper``).
Modules with a command-line interface run as ``python -m simsoma.<module>``; the main entry point is
the ``simsoma`` command (``simsoma.cli``). The shared sequencing observation model is the separate
``plantsoma_obs`` package, which is installed together with simSOMA.
"""
__version__ = "0.3.0"
