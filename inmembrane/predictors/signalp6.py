"""
SignalP 6.0 predictor wrapper
-----------------------------
Predicts Sec/Tat signal peptides and cleavage sites.

Per-protein fields added to proteins[seqid]:

    - signalp_is_sp           (bool)
    - signalp_type            (str or None)   # e.g. "Sec/SPI", "Sec/SPII", "Tat/SPII"
    - signalp_cleave_position (int, 0 if unknown; 1-based)
    - signalp_cleave_prob     (float, 0.0 if unknown)
"""

import os
import subprocess

from inmembrane.helpers import log_stdout, log_stderr


# Mapping SignalP 6.0 short labels → descriptive type
_SIG6_TYPE_MAP = {
    "OTHER": None,
    "SP": "Sec/SPI",
    "LIPO": "Sec/SPII",
    "TAT": "Tat/SPI",
    "TATLIPO": "Tat/SPII",
    "PILIN": "Sec/SPIII",
}


def annotate(params, proteins):
    """
    Run SignalP 6.0 on the input FASTA and annotate `proteins` in place.
    """
    signalp_bin = params.get("signalp_bin", "signalp6")
    fasta = params["fasta"]
    signalp_skip_cmd = bool(params.get("signalp_skip_cmd", False))

    out_dir = os.path.abspath(
        os.path.join(params.get("out_dir", "."), "signalp6_out")
    )
    os.makedirs(out_dir, exist_ok=True)

    log_stdout(f"~ Running SignalP 6.0 → {signalp_bin}")
    log_stdout(f"# Output directory: {out_dir}")

    cmd = [
        signalp_bin,
        "--fastafile",
        fasta,
        "--output_dir",
        out_dir,
        "--format",
        "txt",
        "--organism",
        params.get("signalp_organism", "other"),
        "--mode",
        params.get("signalp_mode", "fast"),
        "--bsize",
        str(params.get("signalp_batch", 10)),
        "--write_procs",
        str(params.get("signalp_write_procs", 8)),
        "--torch_num_threads",
        str(params.get("signalp_threads", 8)),
    ]

    if not signalp_skip_cmd:
        try:
            log_stdout(f"# SignalP 6.0 command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
        except FileNotFoundError:
            log_stderr(f"# ERROR: SignalP 6.0 binary not found: {signalp_bin}")
            raise
        except subprocess.CalledProcessError as e:
            log_stderr(
                f"# ERROR: SignalP 6.0 run failed with exit code {e.returncode}"
            )
            raise
    else:
        log_stdout("# SignalP 6.0 run skipped (signalp_skip_cmd=True)")

    # Expected summary file: prediction_results.txt inside out_dir
    result_path = os.path.join(out_dir, "prediction_results.txt")
    if not os.path.isfile(result_path):
        log_stderr(f"# WARNING: SignalP 6.0 result file not found at {result_path}")
        return proteins

    log_stdout("# Parsing SignalP 6.0 prediction_results.txt ...")

    with open(result_path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 2:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue

            raw_id = parts[0]
            seqid = raw_id.split()[0]

            pred_label = parts[1].strip().upper()
            is_sp = pred_label != "OTHER"
            sp_type = _SIG6_TYPE_MAP.get(pred_label, pred_label)

            cleave_pos = 0
            cleave_prob = 0.0

            # Heuristic parsing of "CS pos: X-Y" and "Pr: <prob>"
            if "CS pos:" in line:
                try:
                    frag = line.split("CS pos:", 1)[1]
                    coords = frag.split(".")[0].strip().split()[0]  # "X-Y"
                    cleave_pos = int(coords.split("-")[1])          # right bound (Y)
                    if "Pr:" in frag:
                        cleave_prob = float(
                            frag.split("Pr:", 1)[1].split()[0].strip()
                        )
                except Exception:
                    # Keep defaults if parsing fails
                    pass

            pdata = proteins.get(seqid)
            if pdata is None:
                log_stderr(
                    f"# WARNING: SignalP entry '{raw_id}' did not match any sequence id"
                )
                continue

            pdata["signalp_is_sp"] = bool(is_sp)
            pdata["signalp_type"] = sp_type
            pdata["signalp_cleave_position"] = int(cleave_pos)
            pdata["signalp_cleave_prob"] = float(cleave_prob)

    log_stdout("# SignalP 6.0 annotations complete.\n")
    return proteins
