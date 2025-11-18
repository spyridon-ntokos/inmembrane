"""
Predictor registry for SerraPHIM–inmembrane
-------------------------------------------

Every predictor module must:
    - live in this directory
    - be named <predictor>.py
    - define an `annotate(params, proteins)` function

This registry automatically imports all predictor modules so that:
    from inmembrane.predictors import get_predictor
    p = get_predictor("tmbed")
    p.annotate(params, proteins)
"""

import os
import importlib
from inmembrane.helpers import log_stderr

# Directory where this file lives
PRED_DIR = os.path.dirname(__file__)

# Registered predictor modules
_registry = {}


def _auto_register_predictors():
    """
    Auto-import <name>.py predictor modules in this directory.
    """
    for fname in os.listdir(PRED_DIR):

        # Only python files, skip private and package files
        if not fname.endswith(".py"):
            continue
        if fname.startswith("_"):
            continue
        if fname in ("__init__.py",):
            continue

        pred_name = fname[:-3]            # strip ".py"
        module_name = f"inmembrane.predictors.{pred_name}"

        try:
            module = importlib.import_module(module_name)

            # must define `annotate()`
            if hasattr(module, "annotate"):
                _registry[pred_name] = module
            else:
                log_stderr(f"# WARNING: predictor '{pred_name}' has no annotate() function — skipping.")

        except Exception as e:
            log_stderr(f"# WARNING: failed to import predictor '{pred_name}': {e}")


def get_predictor(name):
    """
    Retrieve a predictor module by name.
    """
    if not _registry:
        _auto_register_predictors()

    if name not in _registry:
        raise KeyError(f"Unknown predictor: {name}")

    return _registry[name]


# Initialize registry at import time
_auto_register_predictors()
