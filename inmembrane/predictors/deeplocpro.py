"""
DeepLocPro predictor wrapper
----------------------------
Predicts subcellular localization categories for bacterial proteins.

Sets per-protein fields:
    - deeplocpro_label  (str)   → main predicted class
    - deeplocpro_probs  (dict)  → per-class probabilities
    - deeplocpro_conf   (float) → probability of the top class
"""

import os
import re
import csv
import subprocess
from inmembrane.helpers import log_stdout, log_stderr


def annotate(params, proteins):
    deeplocpro_bin = params.get("deeplocpro_bin", "deeplocpro")
    fasta = params["fasta"]
    skip_cmd = True  # set True for dry-run testing

    out_dir = os.path.abspath(os.path.join(params.get("out_dir", "."), "deeplocpro_out"))
    os.makedirs(out_dir, exist_ok=True)

    log_stdout(f"~ Running DeepLocPro → {deeplocpro_bin}")
    log_stdout(f"# Output directory: {out_dir}")

    # -----------------------------------------------------------
    # Build CLI command
    # -----------------------------------------------------------
    cmd = [
        deeplocpro_bin,
        "-f", fasta,
        "-o", out_dir,
        "-d", params.get("deeplocpro_device", "cpu"),
    ]
    if "deeplocpro_group" in params:
        cmd += ["--group", params["deeplocpro_group"]]

    # -----------------------------------------------------------
    # Run DeepLocPro
    # -----------------------------------------------------------
    if not skip_cmd:
        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError:
            log_stderr(f"# ERROR: DeepLocPro binary not found: {deeplocpro_bin}")
            raise
        except subprocess.CalledProcessError as e:
            log_stderr(f"# ERROR: DeepLocPro failed (exit {e.returncode})")
            raise

    # -----------------------------------------------------------
    # Locate and parse the output CSV
    # -----------------------------------------------------------

    # Try to locate a timestamped results CSV like results_YYYYMMDD-HHMMSS.csv
    pattern = re.compile(r"^results_\d{8}-\d{6}\.csv$")
    csv_candidates = [f for f in os.listdir(out_dir) if pattern.match(f)]
    if csv_candidates:
        csv_candidates.sort(reverse=True)  # newest first
        result_path = os.path.join(out_dir, csv_candidates[0])
        log_stdout(f"# Found DeepLocPro results file: {os.path.basename(result_path)}")
    else:
        # fallback to any .csv in out_dir if pattern not matched
        all_csvs = [f for f in os.listdir(out_dir) if f.endswith(".csv")]
        if all_csvs:
            all_csvs.sort(reverse=True)
            result_path = os.path.join(out_dir, all_csvs[0])
            log_stdout(f"# Found DeepLocPro results file: {os.path.basename(result_path)}")
        else:
            log_stderr(f"# WARNING: DeepLocPro output CSV not found in {out_dir}")
            return proteins

    log_stdout("# Parsing DeepLocPro prediction CSV ...")

    # Expected header:
    # ,ACC,Localization,Cell wall & surface,Extracellular,Cytoplasmic,
    #  Cytoplasmic Membrane,Outer Membrane,Periplasmic
    class_names = [
        "Cell wall & surface",
        "Extracellular",
        "Cytoplasmic",
        "Cytoplasmic Membrane",
        "Outer Membrane",
        "Periplasmic",
    ]

    with open(result_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            seqid = row.get("ACC") or row.get("Accession") or row.get("ID")
            if not seqid:
                continue
            seqid = seqid.split()[0]

            # Normalized localization label
            label = row.get("Localization", "").replace(" ", "_").lower()

            # Extract per-class probabilities
            probs = {
                cls.replace(" ", "_").lower(): float(row.get(cls, 0.0) or 0.0)
                for cls in class_names if cls in row
            }
            conf = max(probs.values()) if probs else 0.0

            pdata = proteins.get(seqid)
            if pdata is None:
                log_stderr(f"# WARNING: DeepLocPro entry '{seqid}' not found in proteins dict")
                continue

            pdata["deeplocpro_label"] = label
            pdata["deeplocpro_probs"] = probs
            pdata["deeplocpro_confidence"] = conf


    log_stdout("# DeepLocPro annotations complete.\n")
    return proteins
