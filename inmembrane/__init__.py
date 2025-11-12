"""
Modernized inmembrane core — plugin-free, Python 3+
"""
__version__ = "0.96.0-dev"

import os, sys, shutil, importlib, codecs, textwrap, json
from collections import OrderedDict
from .helpers import dict_get, create_proteins_dict, log_stdout, log_stderr

# Optional YAML config support
try:
    import yaml
except ImportError:
    yaml = None

module_dir = os.path.abspath(os.path.dirname(__file__))

# Default config path
DEFAULT_CONFIG_PATH = os.path.expanduser(
    "~/SerraPHIM_v2/tools/inmembrane/inmembrane.config"
)

def get_params(config_path=None):
    """
    Loads configuration parameters from an existing config file.
    If no config is found, creates a new default one at the default path.
    """

    # Determine config path
    config = config_path or DEFAULT_CONFIG_PATH

    if os.path.isfile(config):
        print()
        log_stderr(f"# Loading configuration from {config}")
    else:
        log_stderr(f"# No inmembrane.config found, generating a default one at {config}")
        os.makedirs(os.path.dirname(config), exist_ok=True)
        fh = open(config, 'w')
        fh.write("""{
'fasta': '',
'csv': 'out_default.csv',
'json': 'out_default.json',
'out_dir': 'out_default',
'protocol': 'gram_neg_modern',

'signalp_bin': 'signalp6',
'deeptmhmm_bin': 'deeptmhmm',
'deeplocpro_bin': 'deeplocpro',
'massp_bin': 'massp',
'hmmer_bin': 'hmmer',

'predictors': ['signalp', 'deeptmhmm', 'deeplocpro', 'massp', 'hmmer'],

'signalp_organism': 'gram-',
'massp_threshold_exposed': 0.30,
'internal_exposed_loop_min': 30,
'bomp_clearly_cutoff': 3.0,
'bomp_maybe_cutoff': 1.0,
'tmbetadisc_positive_required_signal': True,

'citations': 'citations.txt'
}""")
        fh.close()

    params = eval(open(config).read())
    params["_config_path"] = config  # store for reference
    return params


def init_output_dir(params):
    """
    Modernized output directory initialization.

    - Ensures 'out_dir' exists.
    - Does NOT copy input FASTA or config.
    - Does NOT change the working directory.
    - Normalizes CSV/JSON paths inside 'out_dir'.
    """
    out_dir = os.path.abspath(params.get("out_dir", "inmembrane_out"))
    os.makedirs(out_dir, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(params["fasta"]))[0]
    params["csv"] = os.path.join(out_dir, f"{base_name}_out.csv")
    params["json"] = os.path.join(out_dir, f"{base_name}_out.json")
    params["citations"] = os.path.join(out_dir, "citations.txt")

    print()
    log_stderr(f"# Output directory: {out_dir}")
    log_stderr(f"# CSV will be saved to: {params['csv']}")
    log_stderr(f"# JSON will be saved to: {params['json']}")
    log_stderr(f"# Citations will be saved to: {params['citations']}\n")

    return out_dir


def import_protocol(params):
    """
    Import the selected Gram+ or Gram– protocol module.
    """
    name = params.get("protocol", "gram_pos")
    module_name = f"inmembrane.protocols.{name}"
    try:
        return importlib.import_module(module_name)
    except Exception as e:
        raise ImportError(f"Cannot import protocol {module_name}: {e}")


def write_json(proteins, out_path):
    with open(out_path, "w") as fh:
        json.dump(proteins, fh, indent=2)


def process(params):
    """
    Core execution loop. Loads protocol, runs each predictor in the registry.
    """
    from .predictors import get_predictor

    protocol = import_protocol(params)
    init_output_dir(params)
    seqids, proteins = create_proteins_dict(params["fasta"])

    predictor_ids = protocol.get_annotations(params)
    log_stdout(f"\n# Predictors to run: {', '.join(predictor_ids)}\n")

    for pid in predictor_ids:
        try:
            predictor = get_predictor(pid)
        except Exception as e:
            log_stderr(f"# WARNING: predictor '{pid}' not found or failed to import ({e}) — skipping.\n")
            continue
        
        try:
            predictor.annotate(params, proteins)
        except Exception as e:
            log_stderr(f"# ERROR running predictor '{pid}': {e}")

    # Protocol post-processing
    for seqid in seqids:
        protocol.post_process_protein(params, proteins[seqid])
        #log_stdout(protocol.protein_output_line(seqid, proteins))

    log_stdout(protocol.summary_table(params, proteins))

    # Write CSV (with header)
    with open(params["csv"], "w") as f:
        # Write header line
        f.write("SeqID,Category,LoopExtent,Details,Name\n")
        # Write each protein line
        for sid in seqids:
            f.write(protocol.protein_csv_line(sid, proteins))
    print()
    log_stderr(f"# CSV written to {params['csv']}")

    # Write JSON
    write_json(proteins, params["json"])
    log_stderr(f"# JSON written to {params['json']}")

    # Write simple citations / run summary
    with open(params['citations'], 'w', encoding='utf-8') as f:
        f.write("# SerraPHIM-inmembrane modern run\n")
        f.write("# Predictors used:\n")
        for pred in params.get("predictors", []):
            f.write(f"- {pred}\n")
        f.write("# (References available in SerraPHIM documentation.)\n")
    log_stderr(f"Citations summary written to {params['citations']}")

    return proteins
