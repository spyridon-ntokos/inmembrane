"""
HMMER / pyHMMER wrapper
-----------------------
Searches sequences against Pfam/Superfamily HMM libraries
to identify surface anchor motifs and conserved domains.

Sets per-protein fields:
    - hmmer_hits (list of str)
"""

from inmembrane.helpers import log_stderr


def annotate(params, proteins):
    log_stderr("# Running HMMER (mock mode) ...")

    # --- Placeholder for pyHMMER integration ---
    # import pyhmmer
    # with pyhmmer.plan7.HMMFile("Pfam-A.hmm") as hmms:
    #     ...

    # --- Mock domain hits ---
    mock_hits = ["PF00746", "PF01473", "PF01476", "PF00877"]  # sample Pfam IDs
    for seqid, pdata in proteins.items():
        # Randomly assign one hit for demonstration
        pdata["hmmer_hits"] = [mock_hits[hash(seqid) % len(mock_hits)]]

    log_stderr("# HMMER annotations complete.\n")
    return proteins