"""
HMMER predictor wrapper
-----------------------
Runs hmmscan against small, curated Pfam panels (beta-barrels, TonB receptors,
cell-wall anchors, veto domains, etc.) and annotates per-protein HMM hits.

This wrapper:

  1. Builds/updates small Pfam HMM libraries on the fly (via InterPro downloads),
  2. Concatenates them into Gram-specific positive and veto HMM databases,
  3. Runs hmmscan on the input FASTA,
  4. Parses the plain-text hmmscan output (stdout) to collect hits.

Per-protein fields set in proteins[seqid]:

    - hmmer_hits          : list of positive Pfam accessions (e.g. ["PF00267", "PF13505"])
    - hmmer_veto_hits     : list of veto Pfam accessions
    - hmmer_raw_hits      : list of dicts, one per significant hit:
                            {
                              "pfam": "PF00267",
                              "name": "Porin_1",
                              "evalue": 1.2e-149,
                              "score": 486.7,
                              "category": "pos" | "veto"
                            }

Configuration keys expected in inmembrane.config:

    - hmmer_bin            : path to hmmscan binary (default: "hmmscan")
    - hmmer_db_root        : base directory for Pfam databases
                             (default: "~/SerraPHIM_v2/tools/hmmer_db")
    - hmmer_cpu            : number of threads for hmmscan (default: 8)
    - hmmer_evalue_cutoff  : max full-sequence E-value to keep a hit (default: 1e-5)
    - hmmer_skip_cmd       : True to skip running hmmscan (dry-run / testing)

    # Pfam panels (lists of Pfam accessions as strings)
    - hmmer_pfam_gram_neg  : positive Pfams for Gram- protocol
    - hmmer_pfam_gram_pos  : positive Pfams for Gram+ protocol
    - hmmer_pfam_veto_neg  : veto Pfams for Gram-
    - hmmer_pfam_veto_pos  : veto Pfams for Gram+

Which panels are used depends on params["protocol"]:
    - "gram_neg_*" → uses hmmer_pfam_gram_neg + hmmer_pfam_veto_neg
    - "gram_pos_*" → uses hmmer_pfam_gram_pos + hmmer_pfam_veto_pos

Author: Spyridon Ntokos (SerraPHIM project)
"""

import os
import gzip
import subprocess
from urllib.request import urlopen, Request

from inmembrane.helpers import log_stdout, log_stderr


INTERPRO_PFAM_URL = (
    "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}?annotation=hmm"
)


def annotate(params, proteins):
    """
    Main entry point: build DB(s) if needed, run hmmscan, parse output,
    and attach HMM hit information to proteins[seqid].
    """
    fasta = params["fasta"]

    hmmer_bin = params.get("hmmer_bin", "hmmscan")
    hmmer_db_root = os.path.expanduser(
        params.get("hmmer_db_root", "~/SerraPHIM_v2/tools/hmmer_db")
    )
    hmmer_cpu = int(params.get("hmmer_cpu", 8))
    e_cut = float(params.get("hmmer_evalue_cutoff", 1e-5))
    skip_cmd = bool(params.get("hmmer_skip_cmd", False))

    protocol_name = params.get("protocol", "").lower()
    is_gram_neg = protocol_name.startswith("gram_neg")
    is_gram_pos = protocol_name.startswith("gram_pos")

    # Select Pfam panels based on protocol
    if is_gram_neg:
        pfam_pos = params.get("hmmer_pfam_gram_neg", []) or []
        pfam_veto = params.get("hmmer_pfam_veto_neg", []) or []
        db_tag = "gramneg"
    elif is_gram_pos:
        pfam_pos = params.get("hmmer_pfam_gram_pos", []) or []
        pfam_veto = params.get("hmmer_pfam_veto_pos", []) or []
        db_tag = "grampos"
    else:
        # Fallback: treat everything as Gram- style
        pfam_pos = params.get("hmmer_pfam_gram_neg", []) or []
        pfam_veto = params.get("hmmer_pfam_veto_neg", []) or []
        db_tag = "generic"

    pfam_pos = list(pfam_pos)
    pfam_veto = list(pfam_veto)

    if not pfam_pos and not pfam_veto:
        log_stdout("# HMMER: no Pfam accessions configured; skipping.\n")
        return proteins

    out_dir = os.path.abspath(os.path.join(params.get("out_dir", "."), "hmmer_out"))
    os.makedirs(out_dir, exist_ok=True)

    log_stdout("~ Running HMMER / hmmscan")
    log_stdout(f"# Output directory: {out_dir}")
    log_stdout(f"# Using protocol:   {protocol_name or 'N/A'}")

    # ------------------------------------------------------------------
    # 1) Ensure HMM databases exist (download Pfams, concatenate, hmmpress)
    # ------------------------------------------------------------------
    pos_hmm_path, pos_name_map = None, {}
    veto_hmm_path, veto_name_map = None, {}

    if pfam_pos:
        pos_hmm_path, pos_name_map = _ensure_hmm_db(
            db_root=hmmer_db_root,
            db_name=f"{db_tag}_pos",
            pfam_list=pfam_pos,
        )

    if pfam_veto:
        veto_hmm_path, veto_name_map = _ensure_hmm_db(
            db_root=hmmer_db_root,
            db_name=f"{db_tag}_veto",
            pfam_list=pfam_veto,
        )

    # If nothing could be built, bail out gracefully
    if not pos_hmm_path and not veto_hmm_path:
        log_stderr("# WARNING: HMMER: no HMM DBs available; skipping.\n")
        return proteins

    # Combined name → accession map (for both pos and veto)
    name_to_acc = {}
    name_to_acc.update(pos_name_map)
    name_to_acc.update(veto_name_map)

    # ------------------------------------------------------------------
    # 2) Run hmmscan on the FASTA for each available DB
    # ------------------------------------------------------------------
    seq_hits = {}  # seqid → { "pos": [...], "veto": [...], "raw": [...] }

    def _run_and_collect(hmm_path, category):
        if not hmm_path:
            return
        txt_out = os.path.join(out_dir, f"{db_tag}_{category}.hmmscan.txt")

        if not skip_cmd:
            cmd = [
                hmmer_bin,
                "--cpu",
                str(hmmer_cpu),
                hmm_path,
                fasta,
            ]
            log_stdout(f"# hmmscan ({category}) command: {' '.join(cmd)}")
            try:
                with open(txt_out, "w", encoding="utf-8") as fh_out:
                    subprocess.run(cmd, check=True, stdout=fh_out)
            except FileNotFoundError:
                log_stderr(f"# ERROR: hmmscan binary not found: {hmmer_bin}")
                raise
            except subprocess.CalledProcessError as e:
                log_stderr(
                    f"# ERROR: hmmscan ({category}) failed with exit code {e.returncode}"
                )
                raise
        else:
            log_stdout(f"# hmmscan ({category}) skipped (hmmer_skip_cmd=True)")

        if os.path.isfile(txt_out):
            hits = _parse_hmmscan_txt(txt_out, name_to_acc, e_cut, category)
            # merge into seq_hits
            for seqid, seq_list in hits.items():
                entry = seq_hits.setdefault(seqid, {"pos": [], "veto": [], "raw": []})
                entry[category].extend(seq_list)
                entry["raw"].extend(seq_list)

    if pos_hmm_path:
        _run_and_collect(pos_hmm_path, "pos")
    if veto_hmm_path:
        _run_and_collect(veto_hmm_path, "veto")

    # ------------------------------------------------------------------
    # 3) Attach hits to proteins dict
    # ------------------------------------------------------------------
    n_with_hits = 0
    for seqid, pdata in proteins.items():
        entry = seq_hits.get(seqid)
        if not entry:
            continue

        # Unique Pfam accessions for pos / veto
        pos_pfams = sorted({h["pfam"] for h in entry["pos"] if h.get("pfam")})
        veto_pfams = sorted({h["pfam"] for h in entry["veto"] if h.get("pfam")})

        pdata["hmmer_hits"] = pos_pfams
        pdata["hmmer_veto_hits"] = veto_pfams
        pdata["hmmer_raw_hits"] = entry["raw"]
        n_with_hits += 1

    log_stdout(f"# HMMER annotations complete. Sequences with hits: {n_with_hits}\n")
    return proteins


# ======================================================================
# Helper functions
# ======================================================================

def _ensure_hmm_db(db_root, db_name, pfam_list):
    """
    Ensure a Pfam HMM database exists:

        db_root/db_name/
            PFxxxxx.hmm   (one per accession)
            db_name.hmm   (concatenated; NAME fields deduplicated)
            db_name.hmm.h3* (hmmpress indices, if hmmpress succeeds)

    We deduplicate by HMM NAME to avoid hmmpress errors like:
        "primary keys not unique: 'Ribosomal_L6' occurs more than once"

    Returns:
        (hmm_path, name_map) where:

            hmm_path : path to concatenated HMM file (str) or None if build failed
            name_map : dict {HMM NAME → Pfam accession}
    """
    db_dir = os.path.join(db_root, db_name)
    os.makedirs(db_dir, exist_ok=True)

    # 1) Ensure individual PFxxxxx.hmm files exist
    for acc in pfam_list:
        acc = acc.strip()
        if not acc:
            continue
        hmm_file = os.path.join(db_dir, f"{acc}.hmm")
        if not os.path.isfile(hmm_file):
            _download_pfam_hmm(acc, hmm_file)

    # 2) Concatenate into db_name.hmm with NAME-based deduplication
    hmm_path = os.path.join(db_dir, f"{db_name}.hmm")

    seen_names = set()
    name_map = {}

    with open(hmm_path, "w", encoding="utf-8") as out_fh:
        first = True
        for acc in pfam_list:
            acc = acc.strip()
            if not acc:
                continue

            hmm_file = os.path.join(db_dir, f"{acc}.hmm")
            if not os.path.isfile(hmm_file):
                continue

            try:
                with open(hmm_file, "r", encoding="utf-8") as in_fh:
                    content = in_fh.read()
            except OSError:
                continue

            # Extract NAME from this HMM
            hmm_name = None
            for line in content.splitlines():
                if line.startswith("NAME"):
                    parts = line.split()
                    if len(parts) >= 2:
                        hmm_name = parts[1].strip()
                    break

            if not hmm_name:
                log_stderr(f"# WARNING: HMM file {hmm_file} missing NAME; skipping.")
                continue

            if hmm_name in seen_names:
                # This is exactly the Ribosomal_L6 situation
                log_stderr(
                    f"# WARNING: duplicate HMM NAME '{hmm_name}' in {hmm_file}; "
                    f"skipping this entry."
                )
                continue

            seen_names.add(hmm_name)
            name_map[hmm_name] = acc

            # Write with a separating blank line between models
            if not first:
                if not content.startswith("\n"):
                    out_fh.write("\n")
            else:
                first = False

            out_fh.write(content)
            if not content.endswith("\n"):
                out_fh.write("\n")

    # If nothing was successfully written, bail out
    if not seen_names:
        log_stderr(
            f"# WARNING: no valid HMMs were written for {db_name}; "
            "skipping this database."
        )
        return None, {}

    # 3) Build hmmpress indices; treat failure as fatal for this DB
    h3m = hmm_path + ".h3m"
    if not os.path.isfile(h3m):
        try:
            subprocess.run(["hmmpress", hmm_path], check=True)
        except FileNotFoundError:
            log_stderr(
                "# WARNING: hmmpress not found; cannot create indices "
                f"for {hmm_path}. This DB will be skipped."
            )
            return None, {}
        except subprocess.CalledProcessError as e:
            log_stderr(
                f"# WARNING: hmmpress failed for {hmm_path} "
                f"with exit code {e.returncode}"
            )
            return None, {}

    return hmm_path, name_map



def _download_pfam_hmm(acc, out_path):
    """
    Download a single Pfam HMM from InterPro and write it to out_path.
    The InterPro endpoint currently returns a .hmm.gz; we transparently
    decompress if needed.
    """
    url = INTERPRO_PFAM_URL.format(acc=acc)
    log_stdout(f"# Downloading Pfam HMM for {acc} from {url}")
    try:
        req = Request(url, headers={"User-Agent": "SerraPHIM-inmembrane"})
        with urlopen(req) as resp:
            data = resp.read()
    except Exception as e:
        log_stderr(f"# ERROR: failed to download {acc} from InterPro: {e}")
        return

    # Detect gzip magic number
    if data.startswith(b"\x1f\x8b"):
        try:
            data = gzip.decompress(data)
        except Exception as e:
            log_stderr(f"# ERROR: failed to decompress HMM for {acc}: {e}")
            return

    try:
        with open(out_path, "wb") as fh:
            fh.write(data)
    except OSError as e:
        log_stderr(f"# ERROR: cannot write HMM file {out_path}: {e}")


def _extract_hmm_name(hmm_path):
    """
    Extract the NAME field from a HMM file (first 'NAME' line).

    Returns:
        str or None
    """
    try:
        with open(hmm_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("NAME"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return parts[1].strip()
                    break
    except OSError:
        return None
    return None


def _parse_hmmscan_txt(txt_path, name_to_acc, evalue_cutoff, category):
    """
    Parse the plain-text hmmscan output (stdout) to collect full-sequence hits.

    We look at the "Scores for complete sequence" table for each Query block.

    Example block:

        Query:       LCBGCM_01883  [L=690]
        Description: TonB-dependent copper receptor
        Scores for complete sequence (score includes all domains):
           --- full sequence ---   --- best 1 domain ---    -#dom-
            E-value  score  bias    E-value  score  bias    exp  N  Model                 Description
            ------- ------ -----    ------- ------ -----   ---- --  --------              -----------
            1.8e-48  154.8  14.4    2.4e-48  154.3  14.4    1.1  1  TonB_dep_Rec_b-barrel  TonB dependent receptor-like, beta-bar

        Domain annotation for each model:
        ...

    We capture for each row with E-value <= evalue_cutoff:

        seqid, model_name, full_evalue, full_score

    and map model_name → Pfam accession via name_to_acc if possible.

    Returns:
        dict: seqid → list of hit dicts:
              { "pfam", "name", "evalue", "score", "category" }
    """
    hits_by_seq = {}
    current_seq = None
    in_scores = False
    header_lines_to_skip = 0

    with open(txt_path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")

            # New query
            if line.startswith("Query:"):
                parts = line.split()
                if len(parts) >= 2:
                    current_seq = parts[1]
                else:
                    current_seq = None
                in_scores = False
                header_lines_to_skip = 0
                continue

            if current_seq is None:
                continue

            if line.startswith("Scores for complete sequence"):
                # Next 3 lines are header; after that, data rows (or [No hits...])
                in_scores = True
                header_lines_to_skip = 3
                continue

            if in_scores:
                if header_lines_to_skip > 0:
                    header_lines_to_skip -= 1
                    continue

                stripped = line.strip()
                if not stripped:
                    # blank line → end of table
                    in_scores = False
                    continue

                if stripped.startswith("[No hits detected"):
                    # explicit no-hit message for this query
                    in_scores = False
                    continue

                if stripped.startswith("Domain annotation for each model"):
                    in_scores = False
                    continue

                # Try to parse a data row
                parts = stripped.split()
                if len(parts) < 9:
                    # too short to be a valid hmmscan summary row
                    continue

                # Columns (0-based):
                #   0: E-value (full sequence)
                #   1: score
                #   2: bias
                #   3: E-value (best 1 domain)
                #   4: score (best 1 domain)
                #   5: bias
                #   6: exp
                #   7: N
                #   8: Model (HMM NAME)
                #   9+: Description (ignored)
                try:
                    full_eval = float(parts[0])
                    full_score = float(parts[1])
                except ValueError:
                    continue

                if full_eval > evalue_cutoff:
                    continue

                model_name = parts[8]
                pfam_acc = name_to_acc.get(model_name, None)

                hit = {
                    "pfam": pfam_acc or model_name,
                    "name": model_name,
                    "evalue": full_eval,
                    "score": full_score,
                    "category": category,
                }

                hits_by_seq.setdefault(current_seq, []).append(hit)

    return hits_by_seq
