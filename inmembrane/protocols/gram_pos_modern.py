"""
Modernized Gram-positive protocol for SerraPHIM–inmembrane
-----------------------------------------------------------

Reimplementation of the original SurfG+/inmembrane Gram-positive
surfaceome logic using modern predictors (SignalP 6.0, DeepTMHMM,
DeepLocPro, MASSP, pyHMMER).

Author: Spyridon Ntokos (SerraPHIM project)
"""

from inmembrane.helpers import dict_get, chop_nterminal_peptide


# ---------------------------------------------------------------------
# Adjustable parameters
# ---------------------------------------------------------------------
DEFAULT_PARAMS = {
    "signalp_organism": "other",
    "massp_threshold_exposed": 0.35,
    "terminal_exposed_loop_min": 50,
    "internal_exposed_loop_min": 100,
}

GRAM_POS_SURFACE_HMMS = {
    "LPXTG": ["PF00746"],
    "GW_repeat": ["SSF0040855", "SSF0040856", "SSF0040857"],
    "PG_binding_1": ["PF01471", "PF08823", "PF09374"],
    "Choline_binding": ["PF01473"],
    "LysM": ["PF01476"],
    "CW_binding_2": ["PF04122"],
    "SLH": ["PF00395"],
    "NLPC_P60": ["PF00877"],
}


def has_surface_anchor(hmm_hits):
    for family, accs in GRAM_POS_SURFACE_HMMS.items():
        if any(acc in hmm_hits for acc in accs):
            return family
    return None


# ---------------------------------------------------------------------
# Core protocol
# ---------------------------------------------------------------------
def get_annotations(params):
    predictors = ["signalp6", "deeptmhmm", "deeplocpro", "massp", "hmmer"]
    params.setdefault("predictors", predictors)
    for k, v in DEFAULT_PARAMS.items():
        params.setdefault(k, v)
    return predictors


def post_process_protein(params, protein):
    is_signalp = dict_get(protein, "signalp_is_sp")
    signalp_type = dict_get(protein, "signalp_type") or ""
    signalp_cut = dict_get(protein, "signalp_cleave_position") or 0
    signalp_pr = dict_get(protein, "signalp_cleave_prob") or 0.0
    is_lipop = any(signalp_type.endswith(x) for x in ("Sec/SPII", "Tat/SPII"))

    tm_segments = dict_get(protein, "deeptmhmm_helices") or []
    n_tms = len(tm_segments)

    loc_label = dict_get(protein, "deeplocpro_label") or "unknown"
    exposed_score = dict_get(protein, "massp_exposed_fraction") or 0.0
    exposed_prob = dict_get(protein, "massp_exposed_prob") or 0.0
    hmm_hits = dict_get(protein, "hmmer_hits") or []

    threshold = params.get("massp_threshold_exposed", 0.35)

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

    if n_tms > 0:
        details.append(f"deeptmhmm({n_tms}TM)")

    if loc_label != "unknown":
        details.append(f"deeplocpro({loc_label})")

    if exposed_score > 0:
        details.append(f"massp(exp={exposed_score:.2f})")

    if hmm_hits:
        details.append(f"hmm({len(hmm_hits)})")

    # -----------------------------------------------------------------
    # Classification logic
    # -----------------------------------------------------------------
    anchor_family = has_surface_anchor(hmm_hits)
    if anchor_family:
        details.append(f"hmm_anchor({anchor_family})")
        category = "PSE-Cellwall"

    elif is_signalp and not n_tms:
        if is_lipop:
            category = "PSE-Lipoprotein"
        else:
            category = "SECRETED"

    elif n_tms > 0:
        if exposed_prob > threshold or exposed_score > threshold:
            category = "PSE-Membrane"
        else:
            category = "MEMBRANE(non-PSE)"

    elif loc_label in ["extracellular", "cell_wall"]:
        category = "PSE-Cellwall"
    elif loc_label in ["cytoplasmic_membrane"]:
        category = "MEMBRANE(non-PSE)"
    elif loc_label == "cytoplasmic":
        category = "CYTOPLASMIC(non-PSE)"
    else:
        category = "CYTOPLASMIC"

    protein["category"] = category
    protein["details"] = details
    protein["loop_extent"] = exposed_score
    return details, category


# ---------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------
def protein_output_line(seqid, proteins):
    p = proteins[seqid]
    return f"{seqid:<15} {p['category']:<20}  Exp={p.get('loop_extent', 0):.2f}  " \
           f"{';'.join(p.get('details', [])):<60} {p.get('name', '')[:60]}"


def protein_csv_line(seqid, proteins, header=False):
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

    out = "\n\n# Number of proteins in each class:\n"
    for c, n in sorted(counts.items()):
        out += f"# {c:<25}\t{n}\n"
    return out
