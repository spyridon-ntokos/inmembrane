"""
SignalP 6.0 wrapper
-------------------
Predicts Sec/Tat signal peptides and cleavage sites.

Sets per-protein fields:
    - signalp_is_sp (bool)
    - signalp_type (str)
    - signalp_cleave_position (int)
"""

import os
import subprocess
from inmembrane.helpers import log_stderr


def annotate(params, proteins):
    signalp_bin = params.get("signalp_bin", "signalp6")
    fasta = "input.fasta"

    log_stderr("# Running SignalP 6.0 ...")

    # --- Placeholder: local CLI execution example ---
    # cmd = f"{signalp_bin} --fastafile {fasta} --organism {params['signalp_organism']} --format txt > signalp.out"
    # subprocess.run(cmd, shell=True, check=True)

    # --- Mock inference loop ---
    for seqid, pdata in proteins.items():
        seq = pdata["seq"]
        # TODO: Replace this logic with real parsing of SignalP output
        pdata["signalp_is_sp"] = seq.startswith("M")
        pdata["signalp_type"] = "Sec/SPI" if pdata["signalp_is_sp"] else None
        pdata["signalp_cleave_position"] = 25 if pdata["signalp_is_sp"] else 0

    log_stderr("# SignalP 6.0 annotations complete.\n")
    return proteins