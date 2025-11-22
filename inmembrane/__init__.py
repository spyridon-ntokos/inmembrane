"""
Modernized inmembrane Python 3+ core
"""

__version__ = "0.96.0-dev"

import os
import importlib
import json

from .helpers import create_proteins_dict, log_stdout, log_stderr

module_dir = os.path.abspath(os.path.dirname(__file__))

# Default config path (can be overridden when calling get_params)
DEFAULT_CONFIG_PATH = os.path.expanduser(
    "~/SerraPHIM_v2/tools/inmembrane/inmembrane.config"
)


def _write_default_config(config_path: str) -> None:
    """
    Write a modern default inmembrane.config compatible with the updated
    SerraPHIM–inmembrane workflow (SignalP 6, TMbed, DeepLocPro, HMMER).

    The file is a Python dict with comments, evaluated via eval().
    """
    os.makedirs(os.path.dirname(config_path), exist_ok=True)
    with open(config_path, "w") as fh:
        fh.write(
            """{
    # === Core ===
    'fasta': '',                         # Path to input FASTA (Bakta .faa)
    'csv': 'out_surfaceome.csv',         # CSV output
    'json': 'out_surfaceome.json',       # JSON output (optional)
    'out_dir': 'inmembrane_out',         # Output directory (will be created if needed)

    # Choose protocol: 'gram_neg_modern' or 'gram_pos_modern'
    'protocol': 'gram_neg_modern',
    'phage_receptor_veto': True,         # Annotation-based veto for obvious non-receptors

    # === SignalP 6.0 ===
    'signalp_bin': 'signalp6',
    'signalp_organism': 'other',         # 'other' is appropriate for bacteria (Gram+/−)
    'signalp_mode': 'fast',
    'signalp_batch': 10,
    'signalp_write_procs': 4,
    'signalp_threads': 4,
    # Set to True only if you want a dry-run without actually calling SignalP
    'signalp_skip_cmd': False,

    # === TMbed ===
    'tmbed_bin': 'tmbed',
    'tmbed_device': 'cpu',               # 'cpu' or 'cuda'
    'tmbed_use_gpu': False,              # set True if you have a compatible GPU
    'tmbed_cpu_fallback': True,
    'tmbed_threads': 4,
    'tmbed_batch_size': 2000,
    'tmbed_batch_size_long': 200,
    'tmbed_max_gpu_len': 4000,           # aa threshold for "long" sequences
    'tmbed_chunk_overlap': 500,          # aa overlap when chunking very long sequences
    'tmbed_skip_cmd': False,             # debugging / dry-run only
    'tmbed_debug': False,
    'tmbed_delete_embedding': False,

    # TMbed-derived topology / exposure thresholds (used in both Gram− and Gram+)
    'tmbed_outside_loop_min_pos': 50,    # aa; strict threshold for Gram+ bacteria
    'tmbed_outside_loop_min_neg': 14,    # aa; min contiguous outside loop for Gram−
    'tmbed_outside_frac_thr': 0.30,      # fraction 'o' to call IM(peri+cyto)/PSE-Membrane
    'tmbed_secreted_outside_frac': 0.80, # fraction 'o' to call SECRETED when signal present

    # === DeepLocPro ===
    'deeplocpro_bin': 'deeplocpro',
    'deeplocpro_group': 'negative',      # 'negative' (Gram−) or 'positive' (Gram+)
    'deeplocpro_max_gpu_len': 4000,
    'deeplocpro_truncate_len': 6000,
    'deeplocpro_device_short': 'cpu',    # 'cpu', 'cuda' or 'mps'
    'deeplocpro_device_long': 'cpu',
    'deeplocpro_skip_cmd': False,        # debugging / dry-run only

    # === HMMER / pyHMMER ===
    'hmmer_bin': 'hmmscan',
    'hmmer_db_root': 'hmmer_db',         # Root directory for Pfam/custom HMM DBs
    'hmmer_cpu': 4,
    'hmmer_evalue_cutoff': 1e-5,
    'hmmer_skip_cmd': False,

    # Custom HMM download URLs (optional; can be left empty)
    'hmmer_custom_hmm_urls': [],

    # === PhReD ===
    # If True, very conservative PhReD-based upgrades for pilins/flagellin are applied
    'phred_appendage_pse': True,

    # === Predictor list (order of execution) ===
    'predictors': ['signalp6', 'tmbed', 'deeplocpro', 'hmmer'],

    # === Metadata ===
    'citations': 'citations.txt',
}
"""
        )


def get_params(config_path=None):
    """
    Load configuration parameters from an existing config file.

    If no config is found at the chosen path, a new default one is created
    (Python dict with comments, evaluated via eval) and then loaded.
    """
    config = config_path or DEFAULT_CONFIG_PATH

    if os.path.isfile(config):
        print()
        log_stderr(f"# Loading configuration from {config}")
    else:
        log_stderr(f"# No inmembrane.config found, generating a default one at {config}")
        _write_default_config(config)

    with open(config, "r") as fh:
        config_src = fh.read()

    # NOTE: we intentionally keep eval() here because the config is a Python dict
    # with comments. If you later move to pure JSON/YAML, this can be replaced.
    try:
        params = eval(config_src, {}, {})
    except Exception as e:
        raise RuntimeError(f"Failed to parse configuration file {config}: {e}")

    # Keep track of where the config came from
    params["_config_path"] = config
    return params


def init_output_dir(params):
    """
    Modernized output directory initialization.

    - Ensures 'out_dir' exists.
    - Does NOT copy input FASTA or config.
    - Does NOT change the working directory.
    - Normalizes CSV/JSON paths inside 'out_dir' using the FASTA basename.
    """
    out_dir = os.path.abspath(params.get("out_dir", "inmembrane_out"))
    os.makedirs(out_dir, exist_ok=True)

    fasta_path = params.get("fasta", "")
    if fasta_path:
        base_name = os.path.splitext(os.path.basename(fasta_path))[0]
    else:
        base_name = "inmembrane"

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
    Import the selected Gram+ or Gram− protocol module.

    Expected values in params['protocol']:
        - 'gram_neg_modern'
        - 'gram_pos_modern'
    """
    name = params.get("protocol", "gram_neg_modern")
    module_name = f"inmembrane.protocols.{name}"
    try:
        return importlib.import_module(module_name)
    except Exception as e:
        raise ImportError(f"Cannot import protocol {module_name}: {e}")


def write_json(proteins, out_path):
    """
    Write full protein annotations (including per-protein fields like
    category, phage_receptor_candidate, details, etc.) as JSON.
    """
    with open(out_path, "w") as fh:
        json.dump(proteins, fh, indent=2)


def process(params):
    """
    Core execution loop.

    - Loads the selected protocol.
    - Runs each predictor in the protocol's predictor list.
    - Invokes protocol.post_process_protein() per sequence.
    - Writes CSV, JSON and a simple citations summary.
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
            log_stderr(
                f"# WARNING: predictor '{pid}' not found or failed to import ({e}) — skipping.\n"
            )
            continue

        try:
            predictor.annotate(params, proteins)
        except Exception as e:
            log_stderr(f"# ERROR running predictor '{pid}': {e}")

    # Protocol post-processing
    for seqid in seqids:
        protocol.post_process_protein(params, proteins[seqid])
        # If you want per-protein logging, uncomment:
        # log_stdout(protocol.protein_output_line(seqid, proteins))

    log_stdout(protocol.summary_table(params, proteins))

    # -----------------------------------------------------------------
    # Write CSV (with header)
    # -----------------------------------------------------------------
    with open(params["csv"], "w") as f:
        # Header is kept simple; PRC/PSE etc. live in the Details column.
        f.write("SeqID,Category,LoopExtent,Details,Name\n")
        for sid in seqids:
            f.write(protocol.protein_csv_line(sid, proteins))

    print()
    log_stderr(f"# CSV written to {params['csv']}")

    # -----------------------------------------------------------------
    # Write JSON
    # -----------------------------------------------------------------
    write_json(proteins, params["json"])
    log_stderr(f"# JSON written to {params['json']}")

    # -----------------------------------------------------------------
    # Write simple citations / run summary
    # -----------------------------------------------------------------
    with open(params["citations"], "w", encoding="utf-8") as f:
        f.write("# SerraPHIM-inmembrane modern run\n")
        f.write("# Predictors used:\n")
        for pred in params.get("predictors", []):
            f.write(f"- {pred}\n")
        f.write("# (References available in SerraPHIM documentation.)\n")

    log_stderr(f"Citations summary written to {params['citations']}")

    return proteins
