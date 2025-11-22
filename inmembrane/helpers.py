#
# Common helper functions for inmembrane
#

import os
import re
import sys
import subprocess
import textwrap
from collections import OrderedDict

# BeautifulSoup is legacy-only (html2text); make it optional
try:
    from bs4 import BeautifulSoup  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    BeautifulSoup = None

LOG_SILENT = False


def dict_get(this_dict, prop):
    """
    Convenience helper to get values from dicts.

    Returns the value for 'prop' if present, otherwise returns False
    (legacy sentinel used throughout the original inmembrane codebase).
    """
    if prop not in this_dict:
        return False
    return this_dict[prop]


def run_with_output(cmd):
    """
    Run an external program as a child process and capture stdout as bytes.

    NOTE: Legacy helper, kept for backwards compatibility. New code should
    prefer subprocess.run(..., check=True, capture_output=True).
    """
    p = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, _ = p.communicate()
    return out


def run(cmd, out_file=None):
    """
    Legacy wrapper to run an external program, optionally redirecting stdout.

    - Checks that the binary exists (either as a file or via `which`).
    - Skips execution if out_file already exists.
    """
    if out_file:
        log_stderr(f"# {cmd} > {out_file}")
    else:
        log_stderr(f"# {cmd}")

    # Skip if output already exists
    if out_file and os.path.isfile(out_file):
        log_stderr(f"# -> skipped: {out_file} already exists")
        return

    binary = cmd.split()[0]
    is_binary_there = False

    if os.path.isfile(binary):
        is_binary_there = True
    if run_with_output(f"which {binary}"):
        is_binary_there = True

    if not is_binary_there:
        raise IOError(f"Couldn't find executable binary '{binary}'")

    if out_file:
        os.system(cmd + " > " + out_file)
    else:
        os.system(cmd)


def silence_log(b):
    """
    Turn logging silent mode on/off for log_stderr and log_stdout.
    """
    global LOG_SILENT
    LOG_SILENT = bool(b)


def log_stderr(s, width=76, comment=True):
    """
    Wrapper for all stderr output. Allows future customization.
    """
    if LOG_SILENT:
        return

    if s and not s.endswith("\n"):
        s += "\n"
    if not s.startswith("#"):
        s = "# " + s
    sys.stderr.write(s)


def log_stdout(s, width=76):
    """
    Wrapper for all stdout output. Allows future customization.
    """
    if LOG_SILENT:
        return
    print(s)


def parse_fasta_header(header):
    """
    Parse a FASTA header (with or without leading '>').

    Returns:
        (seqid, name)

    If an NCBI-style header (gi|id|db|accession ...) is detected, uses
    "gi|ginumber" as seqid and the final pipe-separated field as name.
    Otherwise:
        - seqid = first whitespace-separated token
        - name  = full header (without '>')
    """
    if not header:
        raise ValueError("Empty FASTA header")

    if header[0] == ">":
        header = header[1:]

    header = header.rstrip("\r\n")
    tokens = header.split("|")

    # NCBI-style: short prefix + '|' and small token as first element
    if "|" in header and len(tokens[0]) <= 3:
        # "gi|ginumber|gb|accession bla" → seqid "gi|ginumber"
        seqid = f"{tokens[0]}|{tokens[1].split()[0]}"
        name = tokens[-1].strip()
    else:
        parts = header.split()
        seqid = parts[0] if parts else header
        name = header.strip()

    return seqid, name


def seqid_to_filename(seqid):
    """
    Make a sequence id filename-friendly (replace '|' with '_').
    """
    return seqid.replace("|", "_")


def create_proteins_dict(fasta):
    """
    Build the main proteins dictionary used by inmembrane from a FASTA file.

    Returns:
        seqids   : list of sequence IDs in the order encountered in FASTA
        proteins : OrderedDict mapping seqid -> protein dict:
            {
                'seq'            : sequence string,
                'name'           : header-derived name/description,
                'original_header': raw header line without '>',
                'sequence_length': int,
                'safe_seqid'     : safe ID (added by generate_safe_seqids)
            }
    """
    seqids = []
    proteins = OrderedDict()

    current_id = None
    seq_chunks = []

    with open(fasta, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue

            if line.startswith(">"):
                # Flush previous sequence
                if current_id is not None:
                    seq = "".join(seq_chunks)
                    proteins[current_id]["seq"] = seq
                    proteins[current_id]["sequence_length"] = len(seq)

                # Start new entry
                seqid, name = parse_fasta_header(line)
                seqids.append(seqid)
                proteins[seqid] = {
                    "seq": "",
                    "name": name,
                    "original_header": line[1:].strip(),
                }
                current_id = seqid
                seq_chunks = []
            else:
                if current_id is not None:
                    seq_chunks.append(line.strip())

    # Flush final sequence
    if current_id is not None:
        seq = "".join(seq_chunks)
        proteins[current_id]["seq"] = seq
        proteins[current_id]["sequence_length"] = len(seq)

    proteins, id_mapping = generate_safe_seqids(proteins)

    return seqids, proteins


def print_proteins(proteins):
    """
    Pretty-print the proteins dictionary in a Python-eval-compatible way.
    """
    print("{")
    for seqid, pdata in proteins.items():
        print(f"  '{seqid}': {{")
        for key, value in pdata.items():
            print(f"    '{key}': {repr(value)},")
        print("  },")
    print("}")


def write_proteins_fasta(fasta_filename, proteins, seqids, width=50):
    """
    Write a subset of proteins to a FASTA file.

    Args:
        fasta_filename : path to output FASTA
        proteins       : dict as returned by create_proteins_dict
        seqids         : iterable of seqids to include (in that order)
        width          : line wrap width (aa per line)
    """
    out = proteins_to_fasta(proteins, seqids=seqids, width=width)
    with open(fasta_filename, "w", encoding="utf-8") as fh:
        fh.write(out)


def proteins_to_fasta(proteins, seqids=None, use_safe_seqid=False, width=50):
    """
    Convert proteins dict to FASTA string.

    Args:
        proteins      : dict from create_proteins_dict
        seqids        : list of seqids to output (defaults to all in dict order)
        use_safe_seqid: if True, use proteins[seqid]['safe_seqid'] as header
                        instead of proteins[seqid]['name']
        width         : line wrap width

    Returns:
        str with FASTA-formatted sequences.
    """
    if seqids is None:
        idlist = list(proteins.keys())
    else:
        idlist = list(seqids)

    fasta_out = []
    for seqid in idlist:
        seq = proteins[seqid]["seq"]
        seq_wrap = textwrap.fill(seq, width)
        if use_safe_seqid:
            header = proteins[seqid]["safe_seqid"]
        else:
            header = proteins[seqid]["name"]
        fasta_out.append(f">{header}\n{seq_wrap}\n")

    return "".join(fasta_out)


def chop_nterminal_peptide(protein, i_cut):
    """
    Legacy helper: adjust N-terminal coordinates of loop / helix annotations
    after removing i_cut residues (e.g., signal peptides).

    NOTE: This assumes legacy topology structures of the form:
        protein['*_loops']   : list of (start, end)
        protein['*_helices'] : list of (start, end)

    Current TMbed-based metrics use dicts instead; this function is kept
    for backwards compatibility with old protocols and is not used by
    the modern SerraPHIM workflow.
    """
    protein["sequence_length"] -= i_cut

    # Shift coordinates
    for prop in protein:
        if "_loops" in prop or "_helices" in prop:
            sses = protein[prop]
            for i in range(len(sses)):
                j, k = sses[i]
                sses[i] = (j - i_cut, k - i_cut)

    # Remove or fix up elements that are completely/partially cut
    for prop in protein:
        if "_loops" in prop or "_helices" in prop:
            sses = protein[prop]
            for i in reversed(range(len(sses))):
                j, k = sses[i]
                # completely cut out
                if j <= 0 and k <= 0:
                    del sses[i]
                # new N-terminal
                elif j <= 0 < k:
                    sses[i] = (1, k)
                    # If a TM-helix becomes the new N-terminus, convert it to a loop
                    if "_helices" in prop:
                        program = prop.split("_")[0]
                        for x in protein:
                            if x == f"{program}_loops":
                                new_N_loop = protein[x][0]
                                new_N_loop[0] = 1
                        del sses[i]


def generate_safe_seqids(proteins):
    """
    Add 'safe_seqid' to each protein, mapping problematic characters away.

    Returns:
        (proteins, id_mapping)

        proteins   : same dict, with 'safe_seqid' fields added
        id_mapping : dict mapping safe_seqid -> original seqid
    """
    id_mapping = {}
    count = 0
    for seqid in proteins:
        safe_id = re.sub(r"[^\w]", "", seqid) + f"_{count}"
        id_mapping[safe_id] = seqid
        proteins[seqid]["safe_seqid"] = safe_id
        count += 1

    return proteins, id_mapping


def clean_directory(top, excluded_files):
    """
    Recursively delete a directory tree (like `rm -rf <top>`),
    skipping files whose basename is in excluded_files.
    """
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            if name not in excluded_files:
                os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))


def html2text(page, aggressive=False):
    """
    Convert an HTML page to (rough) plain text.

    NOTE: Legacy helper, only available if BeautifulSoup is installed.
    If bs4 is not available, returns the original page unchanged.
    """
    if BeautifulSoup is None:
        # Fallback: no HTML parsing available
        return page

    soup = BeautifulSoup(page, "html.parser")

    # kill all script and style elements
    for script in soup(["script", "style"]):
        script.extract()

    text = soup.get_text()

    # break into lines and remove leading and trailing space on each
    lines = (line.strip() for line in text.splitlines())
    if aggressive:
        # break multi-headlines into a line each
        chunks = (
            phrase.strip()
            for line in lines
            for phrase in line.split("  ")
        )
        # drop blank lines
        text = "\n".join(chunk for chunk in chunks if chunk)
    else:
        text = "\n".join(lines)

    return text
