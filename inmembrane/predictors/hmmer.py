"""
HMMER predictor wrapper
-----------------------
Runs hmmscan against:
  - small, curated Pfam panels (beta-barrels, TonB receptors, cell-wall anchors,
    veto domains, etc.), and
  - optional *custom HMM collections* (e.g. CapsuleFinder profiles) provided as
    archive URLs with *.hmm files inside.

This wrapper:

  1. Builds/updates small Pfam HMM libraries on the fly (via InterPro downloads),
  2. Concatenates them into Gram-specific positive and veto HMM databases,
  3. For each custom HMM URL:
       - downloads the archive to hmmer_db_root,
       - extracts it,
       - recursively collects all *.hmm files,
       - concatenates them into a single HMM DB,
       - runs hmmpress on that DB,
  4. Runs hmmscan on the input FASTA against:
       - Pfam-positive DB (if any),
       - all custom DBs (treated as positive),
       - Pfam-veto DB (if any),
  5. Parses the plain-text hmmscan output (stdout) to collect hits.

Per-protein fields set in proteins[seqid]:

    - hmmer_hits          : list of positive Pfam/custom model IDs
                            (e.g. ["PF00267", "OMPP1", "CapsuleFinder_modelX"])
    - hmmer_veto_hits     : list of veto Pfam accessions
    - hmmer_raw_hits      : list of dicts, one per significant hit:
                            {
                              "pfam": "PF00267" or "CustomModelName",
                              "name": "HMM_NAME",
                              "evalue": 1.2e-149,
                              "score": 486.7,
                              "category": "pos" | "veto"
                            }

Configuration keys expected in inmembrane.config:

    - hmmer_bin             : path to hmmscan binary (default: "hmmscan")
    - hmmer_db_root         : base directory for HMM databases
                              (default: "~/SerraPHIM_v2/tools/hmmer_db")
    - hmmer_cpu             : number of threads for hmmscan (default: 8)
    - hmmer_evalue_cutoff   : max full-sequence E-value to keep a hit (default: 1e-5)
    - hmmer_skip_cmd        : True to skip running hmmscan (dry-run / testing)

    # Pfam panels (lists of Pfam accessions as strings)
    - hmmer_pfam_gram_neg   : positive Pfams for Gram- protocol
    - hmmer_pfam_gram_pos   : positive Pfams for Gram+ protocol
    - hmmer_pfam_veto_neg   : veto Pfams for Gram-
    - hmmer_pfam_veto_pos   : veto Pfams for Gram+

    # Custom HMM collections:
    # list of URLs pointing to archives (.tar, .tar.gz, .tgz, .zip) that contain
    # one or more *.hmm files in any subdirectory tree. These will be built into
    # separate HMM DBs under hmmer_db_root and treated as additional positive
    # (receptor-like) models for all protocols.
    - hmmer_custom_hmm_urls : list of archive URLs (optional)

Which Pfam panels are used depends on params["protocol"]:
    - "gram_neg_*" → uses hmmer_pfam_gram_neg + hmmer_pfam_veto_neg
    - "gram_pos_*" → uses hmmer_pfam_gram_pos + hmmer_pfam_veto_pos

Author: Spyridon Ntokos (SerraPHIM project)
"""

import os
import re
import gzip
import tarfile
import zipfile
import shutil
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
    log_stdout("~ Running HMMER / hmmscan")
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

    # Custom HMM URLs (e.g. CapsuleFinder profiles)
    custom_urls = params.get("hmmer_custom_hmm_urls", []) or []

    if not pfam_pos and not pfam_veto and not custom_urls:
        log_stdout("# HMMER: no Pfam accessions or custom HMM URLs configured; skipping.\n")
        return proteins

    out_dir = os.path.abspath(os.path.join(params.get("out_dir", "."), "hmmer_out"))
    os.makedirs(out_dir, exist_ok=True)

    log_stdout(f"# Output directory: {out_dir}")
    log_stdout(f"# Using protocol:   {protocol_name or 'N/A'}")
    if custom_urls:
        log_stdout(f"# Custom HMM archives configured: {len(custom_urls)}")

    # ------------------------------------------------------------------
    # 1) Ensure HMM databases exist (Pfam + custom archives)
    # ------------------------------------------------------------------
    pos_hmm_path, pos_name_map = None, {}
    veto_hmm_path, veto_name_map = None, {}
    custom_dbs = []       # list of (hmm_path, label)
    custom_name_maps = [] # list of {HMM_NAME -> model_id}

    # Pfam positive DB
    if pfam_pos:
        pos_hmm_path, pos_name_map = _ensure_hmm_db(
            db_root=hmmer_db_root,
            db_name=f"{db_tag}_pos",
            pfam_list=pfam_pos,
        )

    # Pfam veto DB
    if pfam_veto:
        veto_hmm_path, veto_name_map = _ensure_hmm_db(
            db_root=hmmer_db_root,
            db_name=f"{db_tag}_veto",
            pfam_list=pfam_veto,
        )

    # Custom HMM DBs (each URL → separate DB under hmmer_db_root)
    for idx, url in enumerate(custom_urls, start=1):
        hmm_path, name_map = _ensure_custom_hmm_db(
            db_root=hmmer_db_root,
            url=url,
            tag=f"custom{idx}",
        )
        if hmm_path and name_map:
            custom_dbs.append((hmm_path, f"custom{idx}"))
            custom_name_maps.append(name_map)

    # If absolutely nothing could be built, bail out gracefully
    if not pos_hmm_path and not veto_hmm_path and not custom_dbs:
        log_stderr("# WARNING: HMMER: no HMM DBs available; skipping.\n")
        return proteins

    # Combined name → accession/model-id map (Pfam + custom)
    name_to_acc = {}
    name_to_acc.update(pos_name_map)
    name_to_acc.update(veto_name_map)
    for nm in custom_name_maps:
        name_to_acc.update(nm)

    # ------------------------------------------------------------------
    # 2) Run hmmscan on the FASTA for each available DB
    # ------------------------------------------------------------------
    seq_hits = {}  # seqid → { "pos": [...], "veto": [...], "raw": [...] }

    def _run_and_collect(hmm_path, category, label=None):
        """
        Run hmmscan against a given HMM DB and merge hits into seq_hits.

        category: "pos" or "veto" (used in hit dicts and seq_hits keys)
        label:    suffix for output file naming and logging (default=category)
        """
        if not hmm_path:
            return

        lbl = label or category
        txt_out = os.path.join(out_dir, f"{db_tag}_{lbl}.hmmscan.txt")

        if not skip_cmd:
            cmd = [
                hmmer_bin,
                "--cpu",
                str(hmmer_cpu),
                hmm_path,
                fasta,
            ]
            log_stdout(f"# hmmscan ({lbl}) command: {' '.join(cmd)}")
            try:
                with open(txt_out, "w", encoding="utf-8") as fh_out:
                    subprocess.run(cmd, check=True, stdout=fh_out)
            except FileNotFoundError:
                log_stderr(f"# ERROR: hmmscan binary not found: {hmmer_bin}")
                raise
            except subprocess.CalledProcessError as e:
                log_stderr(
                    f"# ERROR: hmmscan ({lbl}) failed with exit code {e.returncode}"
                )
                raise
        else:
            log_stdout(f"# hmmscan ({lbl}) skipped (hmmer_skip_cmd=True)")

        if os.path.isfile(txt_out):
            hits = _parse_hmmscan_txt(txt_out, name_to_acc, e_cut, category)
            # merge into seq_hits
            for seqid, seq_list in hits.items():
                entry = seq_hits.setdefault(seqid, {"pos": [], "veto": [], "raw": []})
                entry[category].extend(seq_list)
                entry["raw"].extend(seq_list)

    # Pfam positive DB
    if pos_hmm_path:
        _run_and_collect(pos_hmm_path, "pos", label="pos")

    # Custom positive DBs – also counted as "pos"
    for hmm_path, custom_label in custom_dbs:
        _run_and_collect(hmm_path, "pos", label=f"{custom_label}_pos")

    # Pfam veto DB
    if veto_hmm_path:
        _run_and_collect(veto_hmm_path, "veto", label="veto")

    # ------------------------------------------------------------------
    # 3) Attach hits to proteins dict
    # ------------------------------------------------------------------
    n_with_hits = 0
    for seqid, pdata in proteins.items():
        entry = seq_hits.get(seqid)
        if not entry:
            continue

        # Unique model IDs for pos / veto
        pos_pfams = sorted({h["pfam"] for h in entry["pos"] if h.get("pfam")})
        veto_pfams = sorted({h["pfam"] for h in entry["veto"] if h.get("pfam")})

        pdata["hmmer_hits"] = pos_pfams
        pdata["hmmer_veto_hits"] = veto_pfams
        pdata["hmmer_raw_hits"] = entry["raw"]
        n_with_hits += 1

    log_stdout(f"# HMMER annotations complete. Sequences with hits: {n_with_hits}\n")
    return proteins


# ======================================================================
# Helper functions (Pfam DBs)
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


# ======================================================================
# Helper functions (custom HMM archives, e.g. CapsuleFinder)
# ======================================================================

def _ensure_custom_hmm_db(db_root, url, tag):
    """
    Given a URL pointing to an archive (tar/tar.gz/tgz/zip) that contains
    one or more *.hmm files (possibly in nested subdirectories), build a
    custom HMM database under:

        db_root/custom_<archive_basename>/

    Steps:
      - Download archive (if not already present),
      - Extract under db_root/custom_<...>/extracted/,
      - Recursively collect all *.hmm files,
      - Concatenate into <tag>.hmm with NAME-based deduplication,
      - Run hmmpress.

    All custom HMMs are treated as **positive** models. We map:

        HMM NAME → model_id  (here model_id = HMM NAME)

    so that the 'pfam' field in hmmer_hits is simply the HMM NAME.
    """
    db_root = os.path.abspath(os.path.expanduser(db_root))

    # Derive a stable, human-readable directory name from the URL
    base_url = url.split("?", 1)[0]
    base_name = os.path.basename(base_url) or "custom"
    safe_base = re.sub(r"[^A-Za-z0-9_.-]", "_", base_name)
    db_name = f"{tag}_{safe_base}"

    db_dir = os.path.join(db_root, db_name)
    os.makedirs(db_dir, exist_ok=True)

    hmm_path = os.path.join(db_dir, f"{tag}.hmm")
    h3m_path = hmm_path + ".h3m"

    # If we already have a pressed DB, just rebuild the name_map and reuse it
    if os.path.isfile(hmm_path) and os.path.isfile(h3m_path):
        log_stdout(f"# Custom HMM DB already present for {url} → {hmm_path}")
        name_map = _build_name_map_from_hmm(hmm_path)
        if not name_map:
            log_stderr(
                f"# WARNING: existing custom DB {hmm_path} has no NAME entries; "
                "it will be ignored."
            )
            return None, {}
        return hmm_path, name_map

    # Otherwise, download & extract the archive
    archive_path = os.path.join(db_dir, safe_base)
    if not os.path.isfile(archive_path):
        log_stdout(f"# Downloading custom HMM archive from {url}")
        try:
            req = Request(url, headers={"User-Agent": "SerraPHIM-inmembrane"})
            with urlopen(req) as resp:
                data = resp.read()
        except Exception as e:
            log_stderr(f"# ERROR: failed to download custom HMM archive: {e}")
            return None, {}

        try:
            with open(archive_path, "wb") as fh:
                fh.write(data)
        except OSError as e:
            log_stderr(f"# ERROR: cannot write custom archive {archive_path}: {e}")
            return None, {}

    # Clean extraction directory to avoid stale HMMs
    extract_dir = os.path.join(db_dir, "extracted")
    if os.path.isdir(extract_dir):
        shutil.rmtree(extract_dir, ignore_errors=True)
    os.makedirs(extract_dir, exist_ok=True)

    # Extract archive (tar/tar.gz/tgz or zip)
    if zipfile.is_zipfile(archive_path):
        try:
            with zipfile.ZipFile(archive_path, "r") as zf:
                zf.extractall(extract_dir)
        except Exception as e:
            log_stderr(f"# ERROR: failed to extract ZIP archive {archive_path}: {e}")
            return None, {}
    else:
        # tarfile.open(..., 'r:*') auto-detects compression (tar, tar.gz, tgz, etc.)
        try:
            with tarfile.open(archive_path, "r:*") as tf:
                tf.extractall(extract_dir)
        except tarfile.ReadError as e:
            log_stderr(f"# ERROR: failed to extract TAR archive {archive_path}: {e}")
            return None, {}
        except Exception as e:
            log_stderr(f"# ERROR: unexpected error extracting {archive_path}: {e}")
            return None, {}

    # Recursively collect *.hmm files
    hmm_files = []
    for root, dirs, files in os.walk(extract_dir):
        for fname in files:
            if fname.endswith(".hmm"):
                hmm_files.append(os.path.join(root, fname))

    if not hmm_files:
        log_stderr(
            f"# WARNING: no *.hmm files found in custom archive {url} "
            f"after extraction to {extract_dir}"
        )
        return None, {}

    log_stdout(
        f"# Custom HMM archive {url}: found {len(hmm_files)} *.hmm files "
        f"→ building {hmm_path}"
    )

    # Concatenate all *.hmm into a single DB with NAME-based deduplication
    seen_names = set()
    name_map = {}
    first = True

    try:
        with open(hmm_path, "w", encoding="utf-8") as out_fh:
            for hmm_file in sorted(hmm_files):
                try:
                    with open(hmm_file, "r", encoding="utf-8") as in_fh:
                        content = in_fh.read()
                except OSError:
                    log_stderr(f"# WARNING: cannot read custom HMM {hmm_file}; skipping.")
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
                    log_stderr(
                        f"# WARNING: custom HMM file {hmm_file} missing NAME; skipping."
                    )
                    continue

                if hmm_name in seen_names:
                    log_stderr(
                        f"# WARNING: duplicate HMM NAME '{hmm_name}' in {hmm_file}; "
                        f"skipping this entry."
                    )
                    continue

                seen_names.add(hmm_name)
                # For custom models, use NAME itself as identifier
                name_map[hmm_name] = hmm_name

                # Separate models with a blank line
                if not first:
                    if not content.startswith("\n"):
                        out_fh.write("\n")
                else:
                    first = False

                out_fh.write(content)
                if not content.endswith("\n"):
                    out_fh.write("\n")
    except OSError as e:
        log_stderr(f"# ERROR: cannot write concatenated custom HMM DB {hmm_path}: {e}")
        return None, {}

    if not seen_names:
        log_stderr(
            f"# WARNING: no valid HMMs were written for custom DB {url}; "
            "skipping this database."
        )
        return None, {}

    # Press the DB
    try:
        subprocess.run(["hmmpress", hmm_path], check=True)
    except FileNotFoundError:
        log_stderr(
            "# WARNING: hmmpress not found; cannot create indices "
            f"for custom DB {hmm_path}. This DB will be skipped."
        )
        return None, {}
    except subprocess.CalledProcessError as e:
        log_stderr(
            f"# WARNING: hmmpress failed for custom DB {hmm_path} "
            f"with exit code {e.returncode}"
        )
        return None, {}

    return hmm_path, name_map


def _build_name_map_from_hmm(hmm_path):
    """
    Build a {HMM NAME → model_id} map from an existing concatenated HMM file.
    For custom DBs we simply set model_id = NAME.
    """
    name_map = {}
    try:
        with open(hmm_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("NAME"):
                    parts = line.split()
                    if len(parts) >= 2:
                        name = parts[1].strip()
                        # For custom DBs, we don't have Pfam accessions; use NAME
                        name_map.setdefault(name, name)
    except OSError as e:
        log_stderr(f"# WARNING: cannot rebuild NAME map from {hmm_path}: {e}")
    return name_map


# ======================================================================
# hmmscan output parser
# ======================================================================

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

    and map model_name → Pfam accession / custom model ID via name_to_acc if possible.

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
