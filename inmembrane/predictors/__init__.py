"""
Predictor registry for SerraPHIM–inmembrane
-------------------------------------------

Each predictor module:

  - lives in this directory,
  - is named <predictor>.py,
  - defines an `annotate(params, proteins)` function.

Usage
-----

    from inmembrane.predictors import get_predictor
    predictor = get_predictor("tmbed")
    predictor.annotate(params, proteins)
"""

import os
import importlib

from inmembrane.helpers import log_stderr

PRED_DIR = os.path.dirname(__file__)
_registry = {}


def _auto_register_predictors():
    """
    Auto-import <name>.py predictor modules in this directory.
    """
    for fname in sorted(os.listdir(PRED_DIR)):
        if not fname.endswith(".py"):
            continue
        if fname.startswith("_") or fname == "__init__.py":
            continue

        pred_name = fname[:-3]  # strip ".py"
        module_name = f"inmembrane.predictors.{pred_name}"

        try:
            module = importlib.import_module(module_name)
        except Exception as e:
            log_stderr(f"# WARNING: failed to import predictor '{pred_name}': {e}")
            continue

        if hasattr(module, "annotate"):
            _registry[pred_name] = module
        else:
            log_stderr(
                f"# WARNING: predictor '{pred_name}' has no annotate() function — skipping."
            )


def get_predictor(name):
    """
    Retrieve a predictor module by name (e.g. 'signalp6', 'tmbed', 'deeplocpro', 'hmmer').
    """
    if not _registry:
        _auto_register_predictors()

    try:
        return _registry[name]
    except KeyError:
        raise KeyError(f"Unknown predictor: {name!r}")
