"""
DeepLocPro wrapper
------------------
Predicts subcellular localization categories for bacterial proteins.

Sets per-protein fields:
    - deeplocpro_label (str)
"""

from inmembrane.helpers import log_stderr, log_stdout


def annotate(params, proteins):
    log_stdout("## Running DeepLocPro ...")

    # --- Placeholder: model invocation ---
    # e.g., deep_loc_pro.predict(fasta="input.fasta", output="deeplocpro.json")

    # --- Mock output ---
    locs = ["cytoplasmic", "cytoplasmic_membrane", "periplasmic",
            "outer_membrane", "extracellular"]
    for i, (seqid, pdata) in enumerate(proteins.items()):
        pdata["deeplocpro_label"] = locs[i % len(locs)]

    log_stdout("## DeepLocPro annotations complete.\n")
    return proteins