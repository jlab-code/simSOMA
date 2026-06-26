"""plugin_api.py

Lightweight plug-in interface for pipeline_v7 modules.

Design goal
-----------
Allow modules (especially organ formation) to be swapped without changing the pipeline wrapper.
The wrapper only orchestrates events and passes standardized primitive inputs.

A module plug-in can be specified either as:
- an importable module path (e.g. "my_pkg.my_organ_plugin"), or
- a filesystem path to a .py file.

In both cases, the plug-in must define:

    def get_module() -> object

The returned object may implement one or more adapter methods used by the wrapper.
Missing methods are delegated to the default adapter for that module.

No file I/O here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional
import importlib
import importlib.util
import os
import hashlib


def _load_from_file(py_path: str) -> Any:
    py_path = os.path.abspath(py_path)
    if not os.path.isfile(py_path):
        raise FileNotFoundError(py_path)
    # Deterministic unique module name to avoid collisions.
    h = hashlib.blake2b(py_path.encode("utf-8"), digest_size=8).hexdigest()
    mod_name = f"pipeline_v7_plugin_{h}"
    spec = importlib.util.spec_from_file_location(mod_name, py_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load plug-in from file: {py_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


def load_plugin(module_path_or_file: str) -> Any:
    """Import a plug-in and return get_module()."""
    s = str(module_path_or_file)
    mod = None
    if s.endswith(".py") or os.path.sep in s:
        mod = _load_from_file(s)
    else:
        mod = importlib.import_module(s)
    if not hasattr(mod, "get_module"):
        raise AttributeError(f"Plug-in {module_path_or_file} must define get_module()")
    return mod.get_module()


class _DelegatingAdapter:
    def __init__(self, default: Any, override: Any):
        self._default = default
        self._override = override

    def __getattr__(self, name: str) -> Any:
        if hasattr(self._override, name):
            return getattr(self._override, name)
        return getattr(self._default, name)


@dataclass
class ModuleRegistry:
    """Container for module adapters used by the wrapper."""
    self_renewal: Any
    pre_branching: Any
    branching: Any
    organ: Any

    @staticmethod
    def from_paths(defaults: "ModuleRegistry", module_paths: Optional[Dict[str, str]] = None) -> "ModuleRegistry":
        if not module_paths:
            return defaults
        reg = ModuleRegistry(
            self_renewal=defaults.self_renewal,
            pre_branching=defaults.pre_branching,
            branching=defaults.branching,
            organ=defaults.organ,
        )
        for key, path in module_paths.items():
            key_l = str(key).lower().strip()
            obj = load_plugin(str(path))
            if key_l == "self_renewal":
                reg.self_renewal = _DelegatingAdapter(defaults.self_renewal, obj)
            elif key_l == "pre_branching":
                reg.pre_branching = _DelegatingAdapter(defaults.pre_branching, obj)
            elif key_l == "branching":
                reg.branching = _DelegatingAdapter(defaults.branching, obj)
            elif key_l == "organ":
                reg.organ = _DelegatingAdapter(defaults.organ, obj)
            else:
                raise KeyError(f"Unknown module key '{key}'. Allowed: self_renewal, pre_branching, branching, organ.")
        return reg
