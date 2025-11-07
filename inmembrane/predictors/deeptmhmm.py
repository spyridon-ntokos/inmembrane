"""
DeepTMHMM wrapper
-----------------
Predicts α-helical and β-barrel transmembrane topologies.

Sets per-protein fields:
    - deeptmhmm_helices (list of tuples)
"""

from inmembrane.helpers import log_stderr, log_stdout


def annotate(params, proteins):
    log_stdout("## Running DeepTMHMM ...")

    # --- Placeholder: local or API call goes here ---
    # e.g., cmd = f"deeptmhmm --fasta input.fasta --format json > deeptmhmm.json"

    # --- Mock data ---
    for seqid, pdata in proteins.items():
        seq = pdata["seq"]
        # Pretend we found 1 TM helix every ~100 residues
        tmhmm = [(i, i + 20) for i in range(1, len(seq), 100) if i + 20 < len(seq)]
        pdata["deeptmhmm_helices"] = tmhmm

    log_stdout("## DeepTMHMM annotations complete.\n")
    return proteins