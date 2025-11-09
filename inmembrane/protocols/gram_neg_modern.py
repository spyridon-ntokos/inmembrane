"""
Modernized Gram-negative protocol for SerraPHIM–inmembrane
-----------------------------------------------------------

Reimplementation of the original inmembrane Gram-negative surfaceome
analysis pipeline using updated predictors (SignalP 6.0, DeepTMHMM,
DeepLocPro, MASSP, pyHMMER). Incorporates both inner membrane and
outer membrane topologies and modernizes the classification heuristics.

Author: Spyridon Ntokos (SerraPHIM project)
"""

from inmembrane.helpers import dict_get, chop_nterminal_peptide


# ---------------------------------------------------------------------
# Adjustable parameters (default values; can be overridden in config)
# ---------------------------------------------------------------------
DEFAULT_PARAMS = {
    "signalp_organism": "other",
    "massp_threshold_exposed": 0.30,
    "internal_exposed_loop_min": 30,
    "bomp_clearly_cutoff": 3.0,
    "bomp_maybe_cutoff": 1.0,
    "tmbetadisc_positive_required_signal": True,
}


# ---------------------------------------------------------------------
# Protocol core
# ---------------------------------------------------------------------
def get_annotations(params):
    predictors = ["signalp6", "deeptmhmm", "deeplocpro", "massp", "hmmer"]
    params.setdefault("predictors", predictors)
    for k, v in DEFAULT_PARAMS.items():
        params.setdefault(k, v)
    return predictors


def post_process_protein(params, protein):
    """
    Integrate predictor outputs and classify each protein according to
    Gram-negative topology and localization logic.
    """

    # ---- Extract predictor outputs ----
    is_signalp = dict_get(protein, "signalp_is_sp")
    signalp_type = dict_get(protein, "signalp_type") or ""
    signalp_cut = dict_get(protein, "signalp_cleave_position") or 0
    signalp_pr = dict_get(protein, "signalp_cleave_prob") or 0.0

    is_lipop = any(signalp_type.endswith(x) for x in ("Sec/SPII", "Tat/SPII"))
    is_tat = "tat" in signalp_type.lower()

    tm_segments = dict_get(protein, "deeptmhmm_helices") or []
    n_tms = len(tm_segments)

    loc_label = dict_get(protein, "deeplocpro_label") or "unknown"
    exposed_score = dict_get(protein, "massp_exposed_fraction") or 0.0
    exposed_prob = dict_get(protein, "massp_exposed_prob") or 0.0
    hmm_hits = dict_get(protein, "hmmer_hits") or []

    # Outer membrane β-barrel predictions
    bomp_score = dict_get(protein, "bomp_score") or 0.0
    tmbetadisc_pred = dict_get(protein, "tmbetadisc_positive") or False
    num_beta_strands = dict_get(protein, "tmbeta_strands") or 0

    threshold = params.get("massp_threshold_exposed", 0.30)

    details = []
    category = "UNKNOWN"

    # ---- Annotation collation ----
    if is_signalp:
        sp = signalp_type or "SP"
        cs = signalp_cut
        pr = signalp_pr
        if cs:
            details.append(f"signalp({sp};CS={cs};Pr={pr:.2f})")
            chop_nterminal_peptide(protein, cs)
        else:
            details.append(f"signalp({sp})")

    elif is_tat:
        details.append("signalp(Tat-like)")

    if is_lipop:
        details.append("lipoprotein")

    if n_tms > 0:
        details.append(f"deeptmhmm({n_tms}TM)")

    if loc_label != "unknown":
        details.append(f"deeplocpro({loc_label})")

    if exposed_score > 0:
        details.append(f"massp(exp={exposed_score:.2f})")

    if hmm_hits:
        details.append(f"hmm({len(hmm_hits)})")

    if bomp_score > 0:
        details.append(f"bomp({bomp_score:.1f})")

    if tmbetadisc_pred:
        details.append("tmbetadisc(+)")

    if num_beta_strands:
        details.append(f"βstrands({num_beta_strands})")

    # -----------------------------------------------------------------
    # Classification logic
    # -----------------------------------------------------------------
    if (bomp_score >= params["bomp_clearly_cutoff"]) or \
       (tmbetadisc_pred and (not params["tmbetadisc_positive_required_signal"] or is_signalp)):
        category = "OM(barrel)"

    elif is_lipop:
        if "Asp+2" in hmm_hits:
            category = "LIPOPROTEIN(IM)"
        else:
            category = "LIPOPROTEIN(OM)"

    elif n_tms > 0:
        if exposed_score > threshold or exposed_prob > threshold:
            category = "IM(peri+cyto)"
        else:
            category = "IM"

    elif is_signalp and not n_tms:
        category = "PERIPLASMIC/SECRETED"

    elif loc_label in ["outer_membrane"]:
        category = "OM(barrel)"
    elif loc_label in ["periplasmic"]:
        category = "PERIPLASMIC/SECRETED"
    elif loc_label in ["cytoplasmic_membrane"]:
        category = "IM"
    elif loc_label in ["cytoplasmic"]:
        category = "CYTOPLASMIC(non-PSE)"
    else:
        category = "CYTOPLASMIC"

    protein["category"] = category
    protein["details"] = details
    protein["loop_extent"] = exposed_score
    return details, category


# ---------------------------------------------------------------------
# Output formatting helpers
# ---------------------------------------------------------------------
def protein_output_line(seqid, proteins):
    p = proteins[seqid]
    return f"{seqid:<15} {p['category']:<20}  Exp={p.get('loop_extent', 0):.2f}  " \
           f"{';'.join(p.get('details', [])):<60} {p.get('name', '')[:60]}"


def protein_csv_line(seqid, proteins, header=False):
    """
    Returns one CSV line per protein. Adds header if requested.
    """
    lines = []
    if header:
        lines.append("SeqID,Category,LoopExtent,Details,Name\n")
    p = proteins[seqid]
    lines.append(
        f"{seqid},{p['category']},{p.get('loop_extent', 0):.2f},"
        f"{';'.join(p.get('details', []))},\"{p.get('name', '')}\"\n"
    )
    return "".join(lines)


def summary_table(params, proteins):
    counts = {}
    for p in proteins.values():
        cat = p["category"]
        counts[cat] = counts.get(cat, 0) + 1

    out = "\n## Number of proteins in each class:\n"
    for c, n in sorted(counts.items()):
        out += f"~ {c:<25}\t{n}\n"
    return out
