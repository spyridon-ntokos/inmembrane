"""
MASSP wrapper
-------------
Predicts residue-level membrane association and exposure probability.

Sets per-protein fields:
    - massp_exposed_fraction (float)
    - massp_exposed_prob (float)
"""

from inmembrane.helpers import log_stderr


def annotate(params, proteins):
    log_stderr("# Running MASSP ...")

    # --- Placeholder: local model call ---
    # e.g., massp.predict(fasta="input.fasta", output="massp.json")

    # --- Mock output ---
    for seqid, pdata in proteins.items():
        seq = pdata["seq"]
        exposed_fraction = 0.1 + (len(seq) % 5) * 0.15  # 0.1–0.85 mock range
        pdata["massp_exposed_fraction"] = exposed_fraction
        pdata["massp_exposed_prob"] = exposed_fraction

    log_stderr("# MASSP annotations complete.\n")
    return proteins