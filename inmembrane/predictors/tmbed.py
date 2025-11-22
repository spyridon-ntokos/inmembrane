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

    Strategy:
      - Split sequences into "short" (≤ tmbed_max_gpu_len) and "long" (> threshold)
      - For short sequences: one standard TMbed run on GPU
      - For long sequences: split into overlapping chunks, run TMbed on chunks,
        then stitch labels back together by majority vote.
    """

    tmbed_bin = params.get("tmbed_bin", "tmbed")

    # Performance / behavior options
    tmbed_threads = int(params.get("tmbed_threads", 8))
    tmbed_use_gpu = bool(params.get("tmbed_use_gpu", True))
    tmbed_cpu_fallback = bool(params.get("tmbed_cpu_fallback", True))
    tmbed_delete_embedding = bool(params.get("tmbed_delete_embedding", True))
    tmbed_batch_size = int(params.get("tmbed_batch_size", 4000))
    tmbed_batch_size_long = int(params.get("tmbed_batch_size_long", 200))
    tmbed_max_gpu_len = int(params.get("tmbed_max_gpu_len", 4000))
    tmbed_chunk_overlap = int(params.get("tmbed_chunk_overlap", 500))

    # Debugging settings
    tmbed_debug = bool(params.get("tmbed_debug", False))
    tmbed_skip_cmd = bool(params.get("tmbed_skip_cmd", False))

    out_dir = os.path.abspath(os.path.join(params.get("out_dir", "."), "tmbed_out"))
    os.makedirs(out_dir, exist_ok=True)

    short_fasta = os.path.join(out_dir, "tmbed_short.faa")
    chunk_fasta = os.path.join(out_dir, "tmbed_long_chunks.faa")

    short_ids = []
    chunk_map = {}  # chunk_id -> (parent_seqid, offset0)

    # Ensure step is positive
    if tmbed_chunk_overlap >= tmbed_max_gpu_len:
        log_stderr(
            f"# WARNING: tmbed_chunk_overlap ({tmbed_chunk_overlap}) "
            f">= tmbed_max_gpu_len ({tmbed_max_gpu_len}); adjusting overlap."
        )
        tmbed_chunk_overlap = tmbed_max_gpu_len // 2 or 1

    step = tmbed_max_gpu_len - tmbed_chunk_overlap
    if step <= 0:
        step = tmbed_max_gpu_len

    # ------------------------------------------------------------------
    # Split sequences into short and long-chunk FASTAs
    # ------------------------------------------------------------------
    with open(short_fasta, "w") as fh_short, open(chunk_fasta, "w") as fh_long:
        for seqid, pdata in proteins.items():
            seq = pdata.get("seq") or pdata.get("sequence", "")
            if not seq:
                continue
            seqlen = pdata.get("sequence_length", len(seq))

            if seqlen <= tmbed_max_gpu_len:
                short_ids.append(seqid)
                fh_short.write(f">{seqid}\n{seq}\n")
            else:
                # Break into overlapping windows
                chunk_idx = 0
                pos = 0
                while pos < seqlen:
                    end = min(pos + tmbed_max_gpu_len, seqlen)
                    subseq = seq[pos:end]
                    chunk_id = f"{seqid}__chunk{chunk_idx}"
                    fh_long.write(f">{chunk_id}\n{subseq}\n")
                    chunk_map[chunk_id] = (seqid, pos)
                    chunk_idx += 1
                    if end == seqlen:
                        break
                    pos += step

    log_stdout(
        f"# TMbed: {len(short_ids)} sequences in short set, "
        f"{len(set(parent for parent, _ in chunk_map.values()))} long sequences "
        f"chunked into {len(chunk_map)} windows (> {tmbed_max_gpu_len} aa)"
    )

    log_stdout(f"~ Running TMbed → {tmbed_bin}")
    log_stdout(f"# Output directory: {out_dir}")

    # Helper to run TMbed
    def run_tmbed_once(fasta_in, embed_out, pred_out, batch_size, use_gpu=True, allow_cpu_fallback=True):
        """
        Run 'tmbed embed' + 'tmbed predict' on a single FASTA file.
        """
        embed_cmd = [
            tmbed_bin,
            "embed",
            "-f", fasta_in,
            "-e", embed_out,
            "--batch-size", str(batch_size),
        ]
        predict_cmd = [
            tmbed_bin,
            "predict",
            "-f", fasta_in,
            "-e", embed_out,
            "-p", pred_out,
            "--out-format", "4",
        ]

        if use_gpu:
            embed_cmd.append("--use-gpu")
            predict_cmd.append("--use-gpu")
            if allow_cpu_fallback:
                embed_cmd.append("--cpu-fallback")
                predict_cmd.append("--cpu-fallback")
        else:
            embed_cmd += ["--no-use-gpu", "-t", str(tmbed_threads)]
            predict_cmd += ["--no-use-gpu", "-t", str(tmbed_threads)]

        if tmbed_skip_cmd:
            log_stdout(f"# TMbed embed (skipped): {' '.join(embed_cmd)}")
            log_stdout(f"# TMbed predict (skipped): {' '.join(predict_cmd)}")
            return

        try:
            log_stdout(f"# TMbed embed command: {' '.join(embed_cmd)}")
            subprocess.run(embed_cmd, check=True)
        except FileNotFoundError:
            log_stderr(f"# ERROR: TMbed binary not found: {tmbed_bin}")
            raise
        except subprocess.CalledProcessError as e:
            log_stderr(f"# ERROR: TMbed embed failed (exit {e.returncode})")
            raise

        try:
            log_stdout(f"# TMbed predict command: {' '.join(predict_cmd)}")
            subprocess.run(predict_cmd, check=True)
        except subprocess.CalledProcessError as e:
            log_stderr(f"# ERROR: TMbed predict failed (exit {e.returncode})")
            raise

    labels_by_id = {}

    # ------------------------------------------------------------------
    # 1) Run TMbed on short sequences
    # ------------------------------------------------------------------
    if short_ids:
        embed_short = os.path.join(out_dir, "tmbed_embed_short.h5")
        pred_short = os.path.join(out_dir, "tmbed_pred_short_format4.pred")

        try:
            run_tmbed_once(
                fasta_in=short_fasta,
                embed_out=embed_short,
                pred_out=pred_short,
                batch_size=tmbed_batch_size,
                use_gpu=tmbed_use_gpu,
                allow_cpu_fallback=tmbed_cpu_fallback,
            )
        except Exception as e:
            log_stderr(f"# ERROR: TMbed short-set run failed: {e}")

        if os.path.isfile(pred_short):
            log_stdout("# Parsing TMbed short-set prediction file (format 4) ...")
            labels_by_id.update(_parse_tmbed_format4(pred_short))

        if os.path.exists(embed_short) and tmbed_delete_embedding and not tmbed_debug:
            try:
                os.remove(embed_short)
                log_stdout(f"# Removed TMbed embedding file: {embed_short}")
            except OSError as e:
                log_stderr(f"# WARNING: Could not remove TMbed embedding file: {e}")

    # ------------------------------------------------------------------
    # 2) Run TMbed on chunked long sequences
    # ------------------------------------------------------------------
    chunk_labels = {}
    if chunk_map:
        embed_long = os.path.join(out_dir, "tmbed_embed_chunks.h5")
        pred_long = os.path.join(out_dir, "tmbed_pred_chunks_format4.pred")

        try:
            run_tmbed_once(
                fasta_in=chunk_fasta,
                embed_out=embed_long,
                pred_out=pred_long,
                batch_size=tmbed_batch_size_long,
                use_gpu=tmbed_use_gpu,
                allow_cpu_fallback=tmbed_cpu_fallback,
            )
        except Exception as e:
            log_stderr(f"# ERROR: TMbed chunk-set run failed: {e}")

        if os.path.isfile(pred_long):
            log_stdout("# Parsing TMbed chunk-set prediction file (format 4) ...")
            chunk_labels = _parse_tmbed_format4(pred_long)

        if os.path.exists(embed_long) and tmbed_delete_embedding and not tmbed_debug:
            try:
                os.remove(embed_long)
                log_stdout(f"# Removed TMbed embedding file: {embed_long}")
            except OSError as e:
                log_stderr(f"# WARNING: Could not remove TMbed embedding file: {e}")

    # Optional: cleanup temp FASTAs
    if not tmbed_debug:
        for tmp_fa in (short_fasta, chunk_fasta):
            try:
                if os.path.exists(tmp_fa):
                    os.remove(tmp_fa)
            except OSError:
                pass

    # ------------------------------------------------------------------
    # 3) Reconstruct long-sequence labels from chunks
    # ------------------------------------------------------------------
    if chunk_labels and chunk_map:
        long_full_labels = _reconstruct_from_chunks(chunk_labels, chunk_map, proteins)
        labels_by_id.update(long_full_labels)

    # ------------------------------------------------------------------
    # Optional: write a merged TMbed format-4 file for all sequences
    # ------------------------------------------------------------------
    merged_pred_path = os.path.join(out_dir, "tmbed_pred_format4.merged.pred")
    try:
        with open(merged_pred_path, "w", encoding="utf-8") as fh:
            for seqid, pdata in proteins.items():
                labels = labels_by_id.get(seqid)
                if not labels:
                    continue
                seq = pdata.get("seq") or pdata.get("sequence", "")
                if not seq:
                    continue
                fh.write(f">{seqid}\n")
                fh.write(seq + "\n")
                fh.write(labels + "\n")
        log_stdout(f"# Wrote merged TMbed predictions to {merged_pred_path}")
    except OSError as e:
        log_stderr(f"# WARNING: could not write merged TMbed .pred file: {e}")


    # ------------------------------------------------------------------
    # 4) Attach metrics to proteins dict
    # ------------------------------------------------------------------
    n_matched = 0
    for seqid, pdata in proteins.items():
        labels = labels_by_id.get(seqid)
        if not labels:
            continue

        seq_len = pdata.get("sequence_length", len(pdata.get("seq", "")))
        if len(labels) != seq_len:
            log_stderr(
                f"# WARNING: TMbed labels length ({len(labels)}) != sequence_length "
                f"({seq_len}) for {seqid}"
            )

        pdata["tmbed_raw"] = labels

        metrics = _analyze_tmbed_labels(labels)
        pdata.update(metrics)
        n_matched += 1

    log_stdout(f"# TMbed annotations complete. Matched sequences: {n_matched}\n")
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

            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                seq_line = lines[j].strip()

            k = j + 1
            while k < len(lines) and not lines[k].strip():
                k += 1
            if k < len(lines):
                lab_line = lines[k].strip()

            if lab_line:
                labels_by_id[seqid] = lab_line

            i = k + 1
        else:
            i += 1

    return labels_by_id


def _reconstruct_from_chunks(chunk_labels, chunk_map, proteins):
    """
    Given per-chunk labels and chunk→(parent, offset) mapping, reconstruct
    full-length label strings for each parent sequence by majority vote
    across overlapping chunks.
    """
    # parent_id -> list of per-position vote lists
    votes_per_seq = {}

    for chunk_id, labels in chunk_labels.items():
        mapping = chunk_map.get(chunk_id)
        if mapping is None:
            continue
        parent_id, offset0 = mapping
        pdata = proteins.get(parent_id)
        if pdata is None:
            continue
        full_len = pdata.get("sequence_length", len(pdata.get("seq", "")))
        if full_len <= 0:
            continue

        per_pos_votes = votes_per_seq.setdefault(
            parent_id, [list() for _ in range(full_len)]
        )

        for i, ch in enumerate(labels):
            pos = offset0 + i
            if pos >= full_len:
                break
            per_pos_votes[pos].append(ch)

    # Turn votes into final label strings
    full_labels = {}
    for parent_id, votes in votes_per_seq.items():
        chars = []
        for v in votes:
            if not v:
                # no prediction for this residue → assume inside
                chars.append("i")
            elif len(v) == 1:
                chars.append(v[0])
            else:
                chars.append(_majority_label(v))
        full_labels[parent_id] = "".join(chars)

    return full_labels


def _majority_label(labels):
    """
    Pick the most frequent label; break ties using a biologically motivated
    priority: TM segments > signal peptide > outside > inside.
    """
    from collections import Counter

    counts = Counter(labels)
    # Priority order (higher is stronger surface/TM evidence)
    priority = {
        "B": 4,
        "b": 4,
        "H": 3,
        "h": 3,
        "S": 2,
        "o": 1,
        "i": 0,
    }

    best_label = None
    best_count = -1
    best_priority = -1

    for label, cnt in counts.items():
        pr = priority.get(label, 0)
        if cnt > best_count or (cnt == best_count and pr > best_priority):
            best_label = label
            best_count = cnt
            best_priority = pr

    return best_label or "i"


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

    segments = []  # list of (char, start_idx, end_idx)
    cur_char = labels[0]
    seg_start = 0
    for i in range(1, n):
        ch = labels[i]
        if ch != cur_char:
            segments.append((cur_char, seg_start, i - 1))
            cur_char = ch
            seg_start = i
    segments.append((cur_char, seg_start, n - 1))

    # Signal peptide (S)
    signal_segments = [(c, s, e) for (c, s, e) in segments if c == "S"]
    has_signal = len(signal_segments) > 0
    if has_signal:
        s_start = min(s for (_, s, _) in signal_segments) + 1
        s_end = max(e for (_, _, e) in signal_segments) + 1
    else:
        s_start = None
        s_end = None

    # Helices (H, h)
    helices = []
    for c, s, e in segments:
        if c not in ("H", "h"):
            continue
        orientation = "in→out" if c == "H" else "out→in"
        helices.append(
            {
                "start": s + 1,
                "end": e + 1,
                "orientation": orientation,
            }
        )
    n_helices = len(helices)

    # Beta strands (B, b)
    beta_strands = []
    for c, s, e in segments:
        if c not in ("B", "b"):
            continue
        orientation = "in→out" if c == "B" else "out→in"
        beta_strands.append(
            {
                "start": s + 1,
                "end": e + 1,
                "orientation": orientation,
            }
        )
    n_strands = len(beta_strands)
    is_barrel_candidate = n_strands >= 8

    # Loops (i, o)
    MIN_LOOP_LEN = 5
    loops = []
    longest_outside = 0
    longest_inside = 0

    for c, s, e in segments:
        if c not in ("i", "o"):
            continue
        length = e - s + 1

        if c == "i":
            location = "inside"
            if length >= MIN_LOOP_LEN and length > longest_inside:
                longest_inside = length
        else:
            location = "outside"
            if length >= MIN_LOOP_LEN and length > longest_outside:
                longest_outside = length

        if length >= MIN_LOOP_LEN:
            loops.append(
                {
                    "start": s + 1,
                    "end": e + 1,
                    "location": location,
                }
            )

    count_o = labels.count("o")
    outside_fraction = float(count_o) / float(len(labels)) if labels else 0.0

    has_tm = (n_helices > 0 or n_strands > 0)
    if n_strands >= 8:
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
