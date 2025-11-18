"""
TMbed predictor wrapper
-----------------------
Predicts transmembrane alpha-helices, beta-strands, signal peptides and
non-TM regions using the TMbed tool (https://github.com/BernhoferM/TMbed).

This wrapper uses a two-step workflow:
    1. tmbed embed   → create ESM2 embeddings (.h5)
    2. tmbed predict → run topology prediction with --out-format 4

We currently support only out-format 4, which encodes per-residue states as:

    B: TM beta strand (IN → OUT)
    b: TM beta strand (OUT → IN)
    H: TM alpha helix (IN → OUT)
    h: TM alpha helix (OUT → IN)
    S: Signal peptide
    i: Non-TM, inside
    o: Non-TM, outside

Per-protein fields set in proteins[seqid]:

    - tmbed_raw:              full label string (length = sequence_length)
    - tmbed_has_signal:       bool
    - tmbed_signal_start:     int or None (1-based)
    - tmbed_signal_end:       int or None (1-based)

    - tmbed_helices:          list of {start, end, orientation} (1-based)
    - tmbed_n_helices:        int

    - tmbed_beta_strands:     list of {start, end, orientation} (1-based)
    - tmbed_n_strands:        int
    - tmbed_is_barrel_candidate: bool  (heuristic: ≥8 strands)

    - tmbed_loops:            list of {start, end, location} (inside/outside)
    - tmbed_longest_outside_loop: int
    - tmbed_longest_inside_loop:  int
    - tmbed_outside_fraction: float in [0, 1]

    - tmbed_has_tm:           bool
    - tmbed_class:            one of {'beta_barrel','alpha_tm','signal_only','soluble'}

Author: Spyridon Ntokos (SerraPHIM project)
"""

import os
import subprocess

from inmembrane.helpers import log_stdout, log_stderr


def annotate(params, proteins):
    """
    Run TMbed embed + predict (format 4), then parse results and
    annotate the proteins dict in place.
    """

    tmbed_bin = params.get("tmbed_bin", "tmbed")
    fasta = params["fasta"]

    # Performance / behavior options
    tmbed_threads = int(params.get("tmbed_threads", 8))
    tmbed_use_gpu = bool(params.get("tmbed_use_gpu", True))
    tmbed_delete_embedding = bool(params.get("tmbed_delete_embedding", True))
    tmbed_debug = bool(params.get("tmbed_debug", False))
    tmbed_skip_cmd = bool(params.get("tmbed_skip_cmd", False))  # for dry-run/testing

    out_dir = os.path.abspath(os.path.join(params.get("out_dir", "."), "tmbed_out"))
    os.makedirs(out_dir, exist_ok=True)

    embed_path = os.path.join(out_dir, "tmbed_embed.h5")
    pred_path = os.path.join(out_dir, "tmbed_pred_format4.pred")

    log_stdout(f"~ Running TMbed → {tmbed_bin}")
    log_stdout(f"# Output directory: {out_dir}")
    log_stdout(f"# Embedding file:   {embed_path}")
    log_stdout(f"# Prediction file:  {pred_path}")

    # ------------------------------------------------------------------
    # Generate embeddings (tmbed embed)
    # ------------------------------------------------------------------
    embed_cmd = [
        tmbed_bin,
        "embed",
        "-f", fasta,
        "-e", embed_path,
        "-t", str(tmbed_threads),
    ]
    if tmbed_use_gpu:
        embed_cmd += ["--use-gpu", "--cpu-fallback"]

    # ------------------------------------------------------------------
    # Predict topology (tmbed predict, out-format 4)
    # ------------------------------------------------------------------
    predict_cmd = [
        tmbed_bin,
        "predict",
        "-f", fasta,
        "-e", embed_path,
        "-p", pred_path,
        "--out-format", "4",
        "-t", str(tmbed_threads),
    ]
    if tmbed_use_gpu:
        predict_cmd += ["--use-gpu", "--cpu-fallback"]

    if not tmbed_skip_cmd:
        # ---- run embed ----
        try:
            log_stdout(f"# TMbed embed command: {' '.join(embed_cmd)}")
            subprocess.run(embed_cmd, check=True)
        except FileNotFoundError:
            log_stderr(f"# ERROR: TMbed binary not found: {tmbed_bin}")
            raise
        except subprocess.CalledProcessError as e:
            log_stderr(f"# ERROR: TMbed embed failed (exit {e.returncode})")
            raise

        # ---- run predict ----
        try:
            log_stdout(f"# TMbed predict command: {' '.join(predict_cmd)}")
            subprocess.run(predict_cmd, check=True)
        except subprocess.CalledProcessError as e:
            log_stderr(f"# ERROR: TMbed predict failed (exit {e.returncode})")
            raise

    # ------------------------------------------------------------------
    # Check prediction file
    # ------------------------------------------------------------------
    if not os.path.isfile(pred_path):
        log_stderr(f"# WARNING: TMbed prediction file not found at {pred_path}")
        return proteins

    log_stdout("# Parsing TMbed prediction file (format 4) ...")

    # Parse per-sequence topology labels
    labels_by_id = _parse_tmbed_format4(pred_path)

    n_matched = 0
    for seqid, pdata in proteins.items():
        labels = labels_by_id.get(seqid)
        if not labels:
            continue

        # sanity check: length match
        seq_len = pdata.get("sequence_length", len(pdata.get("seq", "")))
        if len(labels) != seq_len:
            log_stderr(
                f"# WARNING: TMbed labels length ({len(labels)}) != sequence_length "
                f"({seq_len}) for {seqid}"
            )

        pdata["tmbed_raw"] = labels

        # analyze label string & update pdata with derived metrics
        metrics = _analyze_tmbed_labels(labels)
        pdata.update(metrics)
        n_matched += 1

    log_stdout(f"# TMbed annotations complete. Matched sequences: {n_matched}\n")

    # ------------------------------------------------------------------
    # Cleanup embedding file (optional)
    # ------------------------------------------------------------------
    if os.path.exists(embed_path) and tmbed_delete_embedding and not tmbed_debug:
        try:
            os.remove(embed_path)
            log_stdout(f"# Removed TMbed embedding file: {embed_path}")
        except OSError as e:
            log_stderr(f"# WARNING: Could not remove TMbed embedding file: {e}")

    return proteins


# ======================================================================
# Helper functions
# ======================================================================

def _parse_tmbed_format4(pred_path):
    """
    Parse TMbed --out-format 4 file.

    Expected structure (per sequence):

        >seqid optional description
        SEQUENCE
        LABELS

    Returns:
        dict: seqid -> label_string
    """
    labels_by_id = {}
    with open(pred_path, encoding="utf-8") as fh:
        lines = [l.rstrip("\n") for l in fh]

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        if line.startswith(">"):
            header = line[1:].strip()
            seqid = header.split()[0]

            # next two non-empty lines should be sequence + labels
            seq_line = ""
            lab_line = ""

            # sequence
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                seq_line = lines[j].strip()

            # labels
            k = j + 1
            while k < len(lines) and not lines[k].strip():
                k += 1
            if k < len(lines):
                lab_line = lines[k].strip()

            if lab_line:
                labels_by_id[seqid] = lab_line

            # continue after labels line
            i = k + 1
        else:
            i += 1

    return labels_by_id


def _analyze_tmbed_labels(labels):
    """
    Given a TMbed label string (format 4), compute a series of
    segment-level and summary metrics.

    Returns a dict suitable to be merged into proteins[seqid].
    """

    n = len(labels)
    if n == 0:
        return {
            "tmbed_has_signal": False,
            "tmbed_signal_start": None,
            "tmbed_signal_end": None,
            "tmbed_helices": [],
            "tmbed_n_helices": 0,
            "tmbed_beta_strands": [],
            "tmbed_n_strands": 0,
            "tmbed_is_barrel_candidate": False,
            "tmbed_loops": [],
            "tmbed_longest_outside_loop": 0,
            "tmbed_longest_inside_loop": 0,
            "tmbed_outside_fraction": 0.0,
            "tmbed_has_tm": False,
            "tmbed_class": "soluble",
        }

    segments = []  # list of (char, start_idx, end_idx) with 0-based indices
    cur_char = labels[0]
    seg_start = 0
    for i in range(1, n):
        ch = labels[i]
        if ch != cur_char:
            segments.append((cur_char, seg_start, i - 1))
            cur_char = ch
            seg_start = i
    segments.append((cur_char, seg_start, n - 1))

    # -------------------------
    # Signal peptide (S)
    # -------------------------
    signal_segments = [(c, s, e) for (c, s, e) in segments if c == "S"]
    has_signal = len(signal_segments) > 0
    if has_signal:
        # merge all S segments into one range
        s_start = min(s for (_, s, _) in signal_segments) + 1  # 1-based
        s_end = max(e for (_, _, e) in signal_segments) + 1
    else:
        s_start = None
        s_end = None

    # -------------------------
    # Helices (H, h)
    # -------------------------
    helices = []
    for c, s, e in segments:
        if c not in ("H", "h"):
            continue
        orientation = "in→out" if c == "H" else "out→in"
        helices.append({
            "start": s + 1,
            "end": e + 1,
            "orientation": orientation,
        })
    n_helices = len(helices)

    # -------------------------
    # Beta strands (B, b)
    # -------------------------
    beta_strands = []
    for c, s, e in segments:
        if c not in ("B", "b"):
            continue
        orientation = "in→out" if c == "B" else "out→in"
        beta_strands.append({
            "start": s + 1,
            "end": e + 1,
            "orientation": orientation,
        })
    n_strands = len(beta_strands)
    is_barrel_candidate = n_strands >= 8  # simple heuristic

    # -----------------------------------------
    # Loops (i, o) with small-segment filtering
    # -----------------------------------------
    MIN_LOOP_LEN = 5  # increase limit to be stricter

    loops = []
    longest_outside = 0
    longest_inside = 0

    for c, s, e in segments:
        if c not in ("i", "o"):
            continue
        length = e - s + 1

        # Always count residues toward outside_fraction,
        # but only treat as a proper loop if length >= MIN_LOOP_LEN
        if c == "i":
            location = "inside"
            if length >= MIN_LOOP_LEN and length > longest_inside:
                longest_inside = length
        else:
            location = "outside"
            if length >= MIN_LOOP_LEN and length > longest_outside:
                longest_outside = length

        if length >= MIN_LOOP_LEN:
            loops.append({
                "start": s + 1,
                "end": e + 1,
                "location": location,
            })

    # fractions (still count all i/o residues)
    count_o = labels.count("o")
    count_i = labels.count("i")
    outside_fraction = float(count_o) / float(len(labels)) if labels else 0.0

    # -------------------------
    # Simplified overall class
    # -------------------------
    has_tm = (n_helices > 0 or n_strands > 0)
    if n_strands >= 4:
        overall_class = "beta_barrel"
    elif has_tm and n_helices > 0:
        overall_class = "alpha_tm"
    elif (not has_tm) and has_signal:
        overall_class = "signal_only"
    else:
        overall_class = "soluble"

    return {
        "tmbed_has_signal": has_signal,
        "tmbed_signal_start": s_start,
        "tmbed_signal_end": s_end,
        "tmbed_helices": helices,
        "tmbed_n_helices": n_helices,
        "tmbed_beta_strands": beta_strands,
        "tmbed_n_strands": n_strands,
        "tmbed_is_barrel_candidate": is_barrel_candidate,
        "tmbed_loops": loops,
        "tmbed_longest_outside_loop": longest_outside,
        "tmbed_longest_inside_loop": longest_inside,
        "tmbed_outside_fraction": outside_fraction,
        "tmbed_has_tm": has_tm,
        "tmbed_class": overall_class,
    }
