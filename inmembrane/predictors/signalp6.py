"""
SignalP 6.0 predictor wrapper
-----------------------------
Predicts Sec/Tat signal peptides and cleavage sites.

Sets per-protein fields:
    - signalp_is_sp (bool)
    - signalp_type (str)
    - signalp_cleave_position (int)
"""

import os
import subprocess
from inmembrane.helpers import log_stdout, log_stderr


def annotate(params, proteins):
    signalp_bin = params.get("signalp_bin", "signalp6")
    fasta = params["fasta"]
    skip_cmd = True

    out_dir = os.path.abspath(os.path.join(params.get("out_dir", "."), "signalp6_out"))
    os.makedirs(out_dir, exist_ok=True)

    log_stdout(f"~ Running SignalP 6.0 → {signalp_bin}")
    log_stdout(f"# Output directory: {out_dir}")

    cmd = [
        signalp_bin,
        "--fastafile", fasta,
        "--output_dir", out_dir,
        "--format", "txt",
        "--organism", params.get("signalp_organism", "other"),
        "--mode", params.get("signalp_mode", "fast"),
        "--bsize", str(params.get("signalp_batch", 10)),
        "--write_procs", str(params.get("signalp_write_procs", 8)),
        "--torch_num_threads", str(params.get("signalp_threads", 8)),
    ]

    if not skip_cmd:
        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError:
            log_stderr(f"# ERROR: SignalP 6.0 binary not found: {signalp_bin}")
            raise
        except subprocess.CalledProcessError as e:
            log_stderr(f"# ERROR: SignalP 6.0 run failed with exit code {e.returncode}")
            raise

    # Expected summary file: prediction_results.txt inside out_dir
    result_path = os.path.join(out_dir, "prediction_results.txt")
    if not os.path.isfile(result_path):
        log_stderr(f"# WARNING: SignalP 6.0 result file not found at {result_path}")
        return proteins

    log_stdout("# Parsing SignalP 6.0 prediction_results.txt ...")

    # Mapping short codes → descriptive signal peptide classes
    type_map = {
        "OTHER": None,
        "SP": "Sec/SPI",
        "LIPO": "Sec/SPII",
        "TAT": "Tat/SPI",
        "TATLIPO": "Tat/SPII",
        "PILIN": "Sec/SPIII",
    }

    with open(result_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 2:
                parts = line.strip().split()

            # The first column is "ID" but often includes the whole header.
            # Use only the first token to match proteins dict keys.
            raw_id = parts[0]
            seqid = raw_id.split()[0]

            pred_label = parts[1].strip().upper()

            is_sp = pred_label != "OTHER"
            type_map = {
                "OTHER": None,
                "SP": "Sec/SPI",
                "LIPO": "Sec/SPII",
                "TAT": "Tat/SPI",
                "TATLIPO": "Tat/SPII",
                "PILIN": "Sec/SPIII",
            }
            sp_type = type_map.get(pred_label, pred_label)
            cleave_pos = 0
            cleave_prob = 0.0

            if "CS pos:" in line:
                try:
                    frag = line.split("CS pos:")[1]
                    coords = frag.split(".")[0].strip().split()[0]  # "X-Y"
                    cleave_pos = int(coords.split("-")[1])          # take right bound
                    if "Pr:" in frag:
                        cleave_prob = float(frag.split("Pr:")[1].split()[0])
                except Exception:
                    pass

            pdata = proteins.get(seqid)
            if pdata is None:
                # If you want a hint while debugging:
                log_stderr(f"# WARNING: SignalP entry '{raw_id}' did not match any sequence id")
                continue

            pdata["signalp_is_sp"] = is_sp
            pdata["signalp_type"] = sp_type
            pdata["signalp_cleave_position"] = cleave_pos
            pdata["signalp_cleave_prob"] = cleave_prob

    log_stdout("# SignalP 6.0 annotations complete.\n")
    return proteins

