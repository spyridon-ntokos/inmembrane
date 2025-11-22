"""
DeepLocPro predictor wrapper
----------------------------
Predicts subcellular localization categories for bacterial proteins.

This wrapper:

  * Splits proteins into "short" and "long" sets based on sequence length.
  * Runs short sequences on a chosen device (GPU/CPU) for speed.
  * Runs long sequences on (usually) CPU, on truncated versions of the sequences,
    to avoid ESM2/DeepLocPro issues with very long proteins.
  * Parses the resulting CSVs and annotates the `proteins` dict.

Per-protein fields set in proteins[seqid]:

    - deeplocpro_label       (str)   → main predicted class (normalized, e.g. "outer_membrane")
    - deeplocpro_probs       (dict)  → per-class probabilities
    - deeplocpro_confidence  (float) → probability of the top class

Config keys used (with typical defaults):

    - deeplocpro_bin            : "deeplocpro"
    - deeplocpro_group          : "negative" / "positive" / "any" / "archaea"
    - deeplocpro_device_short   : device for short set ("cuda", "cpu", "mps")
    - deeplocpro_device_long    : device for long set
    - deeplocpro_max_gpu_len    : length threshold for short vs long (default 4000)
    - deeplocpro_truncate_len   : truncation length for long sequences (default = max_gpu_len)
    - deeplocpro_skip_cmd       : True to skip running DeepLocPro (dry-run)
"""

import os
import re
import csv
import subprocess

from inmembrane.helpers import log_stdout, log_stderr


def annotate(params, proteins):
    """
    Run DeepLocPro on the input FASTA (via split short/long FASTAs)
    and annotate `proteins` in place.
    """
    deeplocpro_bin = params.get("deeplocpro_bin", "deeplocpro")
    fasta = params["fasta"]
    deeplocpro_skip_cmd = bool(params.get("deeplocpro_skip_cmd", False))

    # Device and length parameters
    base_device = params.get("deeplocpro_device", "cuda")
    device_short = params.get("deeplocpro_device_short", base_device)
    device_long = params.get("deeplocpro_device_long", "cpu")

    max_gpu_len = int(params.get("deeplocpro_max_gpu_len", 4000))
    truncate_len = int(params.get("deeplocpro_truncate_len", max_gpu_len))

    group = params.get("deeplocpro_group", None)

    out_root = os.path.abspath(
        os.path.join(params.get("out_dir", "."), "deeplocpro_out")
    )
    os.makedirs(out_root, exist_ok=True)

    # ------------------------------------------------------------------
    # 1) Split sequences into short / long FASTA files
    # ------------------------------------------------------------------
    short_fasta = os.path.join(out_root, "deeplocpro_short.faa")
    long_fasta = os.path.join(out_root, "deeplocpro_long.faa")

    short_ids = []
    long_ids = []

    with open(short_fasta, "w") as fh_short, open(long_fasta, "w") as fh_long:
        for seqid, pdata in proteins.items():
            seq = pdata.get("seq") or pdata.get("sequence", "")
            if not seq:
                continue
            seqlen = pdata.get("sequence_length", len(seq))

            if seqlen <= max_gpu_len:
                short_ids.append(seqid)
                fh_short.write(f">{seqid}\n{seq}\n")
            else:
                long_ids.append(seqid)
                trunc = seq[:truncate_len]
                fh_long.write(f">{seqid}\n{trunc}\n")

    log_stdout(f"~ Running DeepLocPro → {deeplocpro_bin}")
    log_stdout(f"# Output root directory: {out_root}")
    log_stdout(
        f"# DeepLocPro: {len(short_ids)} sequences in short set (≤ {max_gpu_len} aa), "
        f"{len(long_ids)} in long set (> {max_gpu_len} aa; truncated to {truncate_len} aa)"
    )

    # ------------------------------------------------------------------
    # 2) Run DeepLocPro on SHORT set
    # ------------------------------------------------------------------
    short_out_dir = os.path.join(out_root, "short")
    os.makedirs(short_out_dir, exist_ok=True)

    csv_short = None
    if short_ids:
        csv_short = _run_deeplocpro_one(
            fasta_path=short_fasta,
            out_dir=short_out_dir,
            device_label=device_short,
            force_cpu=(device_short == "cpu"),
            deeplocpro_bin=deeplocpro_bin,
            group=group,
            skip_cmd=deeplocpro_skip_cmd,
            tag="short",
        )

    # ------------------------------------------------------------------
    # 3) Run DeepLocPro on LONG set (usually CPU, truncated sequences)
    # ------------------------------------------------------------------
    long_out_dir = os.path.join(out_root, "long")
    os.makedirs(long_out_dir, exist_ok=True)

    csv_long = None
    if long_ids:
        csv_long = _run_deeplocpro_one(
            fasta_path=long_fasta,
            out_dir=long_out_dir,
            device_label=device_long,
            force_cpu=(device_long == "cpu"),
            deeplocpro_bin=deeplocpro_bin,
            group=group,
            skip_cmd=deeplocpro_skip_cmd,
            tag="long_truncated",
        )

    # ------------------------------------------------------------------
    # 4) Parse CSVs and annotate proteins
    # ------------------------------------------------------------------
    if not csv_short and not csv_long:
        log_stderr("# WARNING: DeepLocPro produced no result CSVs; skipping.\n")
        return proteins

    if csv_short:
        _parse_deeplocpro_csv(csv_short, proteins)
    if csv_long:
        _parse_deeplocpro_csv(csv_long, proteins)

    log_stdout("# DeepLocPro annotations complete.\n")
    return proteins


# ======================================================================
# Helper functions
# ======================================================================


def _run_deeplocpro_one(
    fasta_path,
    out_dir,
    device_label,
    force_cpu,
    deeplocpro_bin="deeplocpro",
    group=None,
    skip_cmd=False,
    tag="",
):
    """
    Run DeepLocPro once on a given FASTA and output directory.

    Returns:
        path to the result CSV (str) or None if not found.
    """
    if not os.path.isfile(fasta_path) or os.path.getsize(fasta_path) == 0:
        return None

    cmd = [
        deeplocpro_bin,
        "-f",
        fasta_path,
        "-o",
        out_dir,
        "-d",
        device_label,
    ]
    if group:
        cmd += ["--group", group]

    env = os.environ.copy()
    if force_cpu:
        # Hide GPUs for this DeepLocPro run only
        env["CUDA_VISIBLE_DEVICES"] = ""

    if not skip_cmd:
        tag_str = f"{tag or device_label}"
        log_stdout(f"# DeepLocPro command [{tag_str}]: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, env=env)
        except FileNotFoundError:
            log_stderr(f"# ERROR: DeepLocPro binary not found: {deeplocpro_bin}")
            raise
        except subprocess.CalledProcessError as e:
            log_stderr(f"# ERROR: DeepLocPro failed (exit {e.returncode})")
            return None
    else:
        log_stdout(f"# DeepLocPro run skipped ({tag or device_label}; deeplocpro_skip_cmd=True)")

    # Locate result CSV (DeepLocPro names it results_YYYYMMDD-HHMMSS.csv)
    pattern = re.compile(r"^results_\d{8}-\d{6}\.csv$")
    csv_candidates = [f for f in os.listdir(out_dir) if pattern.match(f)]
    if csv_candidates:
        csv_candidates.sort(reverse=True)  # newest first
        result_path = os.path.join(out_dir, csv_candidates[0])
        log_stdout(
            f"# Found DeepLocPro results file ({tag or device_label}): "
            f"{os.path.basename(result_path)}"
        )
        return result_path

    # Fallback: any .csv inside out_dir
    all_csvs = [f for f in os.listdir(out_dir) if f.endswith(".csv")]
    if all_csvs:
        all_csvs.sort(reverse=True)
        result_path = os.path.join(out_dir, all_csvs[0])
        log_stdout(
            f"# Found DeepLocPro results file ({tag or device_label}): "
            f"{os.path.basename(result_path)}"
        )
        return result_path

    log_stderr(f"# WARNING: DeepLocPro output CSV not found in {out_dir}")
    return None


def _parse_deeplocpro_csv(csv_path, proteins):
    """
    Parse a DeepLocPro results CSV and update proteins[seqid] in place.

    Expected columns (based on original DeepLocPro output):

        ACC,Localization,Cell wall & surface,Extracellular,Cytoplasmic,
        Cytoplasmic Membrane,Outer Membrane,Periplasmic

    We normalize class names by replacing spaces with underscores and lowercasing.
    """
    log_stdout(f"# Parsing DeepLocPro prediction CSV: {os.path.basename(csv_path)}")

    class_names = [
        "Cell wall & surface",
        "Extracellular",
        "Cytoplasmic",
        "Cytoplasmic Membrane",
        "Outer Membrane",
        "Periplasmic",
    ]

    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            seqid = row.get("ACC") or row.get("Accession") or row.get("ID")
            if not seqid:
                continue
            seqid = seqid.split()[0]

            label_raw = row.get("Localization", "") or ""
            label = label_raw.replace(" ", "_").lower() if label_raw else "unknown"

            # Per-class probabilities
            probs = {}
            for cls in class_names:
                if cls in row:
                    key_norm = cls.replace(" ", "_").lower()
                    try:
                        val = float(row.get(cls, 0.0) or 0.0)
                    except ValueError:
                        val = 0.0
                    probs[key_norm] = val

            conf = max(probs.values()) if probs else 0.0

            pdata = proteins.get(seqid)
            if pdata is None:
                log_stderr(
                    f"# WARNING: DeepLocPro entry '{seqid}' not found in proteins dict"
                )
                continue

            prev_conf = pdata.get("deeplocpro_confidence")
            if prev_conf is None:
                prev_conf = -1.0

            # Only overwrite if this result is at least as confident as what we had
            if conf >= prev_conf:
                pdata["deeplocpro_label"] = label
                pdata["deeplocpro_probs"] = probs
                pdata["deeplocpro_confidence"] = conf
