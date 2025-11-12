"""
DeepTMHMM wrapper
-----------------
Predicts α-helical and β-barrel transmembrane topologies.

Sets per-protein fields:
    - deeptmhmm_helices (list of tuples)
"""

from inmembrane.helpers import log_stderr, log_stdout


# Wrapper function with mock predictor
def annotate(params, proteins):
    log_stdout("~ Running DeepTMHMM (mock predictor wrapper)...")

    # --- Placeholder: local or API call goes here ---
    # e.g., cmd = f"deeptmhmm --fasta input.fasta --format json > deeptmhmm.json"

    # --- Mock data ---
    for seqid, pdata in proteins.items():
        seq = pdata["seq"]
        # Pretend we found 1 TM helix every ~100 residues
        tmhmm = [(i, i + 20) for i in range(1, len(seq), 100) if i + 20 < len(seq)]
        pdata["deeptmhmm_helices"] = tmhmm

    log_stdout("# DeepTMHMM annotations complete.\n")
    return proteins

## Testing actual DeepTMHMM (BioLib-wrapped, intended to run locally on GPU)

# import biolib
# from pathlib import Path

# biolib.utils.STREAM_STDOUT = True # Stream progress from app in real time

# # Load app
# deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24')

# # Resolve path (must be absolute, and expand "~")
# test_fasta = Path("~/SerraPHIM_v2/data/bakta_annotations/AP013063.1_Serratia_marcescens_SM39_DNA_+_plasmids/AP013063.1_Serratia_marcescens_SM39_DNA_+_plasmids.faa").expanduser()
# output_path = Path("~/SerraPHIM_v2/data/inmembrane_output/deeptmhmm_test").expanduser()

# # Run job with file upload
# deeptmhmm_job = deeptmhmm.cli(
#     args=f'--fasta {test_fasta}', machine='local'
# )

# # Save results to local directory
# deeptmhmm_job.save_files(output_path)

# print(f"Results saved locally to {output_path}")