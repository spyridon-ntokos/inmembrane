"""
Unified predictor entrypoint — lazy importer (safe for circular imports)
"""

import importlib
import sys

PREDICTOR_NAMES = ["signalp", "deeptmhmm", "deeplocpro", "massp", "hmmer"]

def get_predictor(name):
    """
    Dynamically imports a predictor module by name (lazy import).
    Example:
        predictor = get_predictor('signalp')
        predictor.annotate(params, proteins)
    """
    if name not in PREDICTOR_NAMES:
        raise ValueError(f"Unknown predictor: {name}")

    module_name = f"inmembrane.predictors.{name}"
    if module_name in sys.modules:
        return sys.modules[module_name]

    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as e:
        raise ImportError(f"Predictor module '{name}' not found ({e})")
    except Exception as e:
        raise RuntimeError(f"Error importing predictor '{name}': {e}")
