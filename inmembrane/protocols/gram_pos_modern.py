"""
Modernized Gram-positive protocol for SerraPHIM–inmembrane
-----------------------------------------------------------

Reimplementation of the original SurfG+/inmembrane Gram-positive
surfaceome logic using modern predictors:

    - SignalP 6.0       (Sec/Tat signal peptides)
    - TMbed             (α-helical / β-strand TM segments, signal-like N-termini)
    - DeepLocPro        (coarse localization; group='positive')
    - HMMER / pyHMMER   (Pfam-based Gram+ receptor anchors and veto Pfams)

Main output classes:
    - PSE-Cellwall
    - PSE-Membrane
    - PSE-Lipoprotein
    - SECRETED
    - MEMBRANE(non-PSE)
    - CYTOPLASMIC

Additionally sets a binary:
    - is_pse (bool): whether the protein is considered “phage surface-exposed”
"""

from inmembrane.helpers import dict_get, chop_nterminal_peptide


# ---------------------------------------------------------------------
# Adjustable parameters
# ---------------------------------------------------------------------
DEFAULT_PARAMS = {
    "signalp_organism": "other",

    # TMbed-derived loop thresholds
    "tmbed_outside_loop_min": 50,          # min contiguous 'o'-segment to call exposed loop
    "tmbed_outside_frac_thr": 0.30,        # min outside_fraction to call PSE-Membrane
    "tmbed_secreted_outside_frac": 0.80,   # fraction 'o' to call SECRETED (if signal)

    # DeepLocPro settings
    "deeplocpro_group": "positive",        # prevent OM/periplasm for Gram+
}


# ---------------------------------------------------------------------
# Core protocol
# ---------------------------------------------------------------------
def get_annotations(params):
    """
    Define the predictor list and inject default parameters.
    """
    predictors = ["signalp6", "tmbed", "deeplocpro", "hmmer"]
    params.setdefault("predictors", predictors)
    for k, v in DEFAULT_PARAMS.items():
        params.setdefault(k, v)
    return predictors


def post_process_protein(params, protein):
    """
    Integrate predictor outputs and classify Gram-positive proteins into
    PSE vs non-PSE categories and finer localization classes, with
    HMMER-based upgrades (anchor Pfams) and vetoes (non-receptor Pfams).
    """

    # ---- SignalP 6.0 outputs ----
    is_signalp = dict_get(protein, "signalp_is_sp")
    signalp_type = dict_get(protein, "signalp_type") or ""
    signalp_cut = dict_get(protein, "signalp_cleave_position") or 0
    signalp_pr = dict_get(protein, "signalp_cleave_prob") or 0.0
    signalp_is_lipop = any(signalp_type.endswith(x) for x in ("Sec/SPII", "Tat/SPII"))

    # ---- TMbed outputs ----
    n_helices = dict_get(protein, "tmbed_n_helices") or 0
    n_strands = dict_get(protein, "tmbed_n_strands") or 0
    longest_outside = dict_get(protein, "tmbed_longest_outside_loop") or 0
    outside_fraction = dict_get(protein, "tmbed_outside_fraction") or 0.0
    has_tmbed_signal = dict_get(protein, "tmbed_has_signal") or False
    tmbed_class = dict_get(protein, "tmbed_class") or None
    has_tm = (n_helices + n_strands) > 0

    # ---- DeepLocPro outputs ----
    loc_label = dict_get(protein, "deeplocpro_label") or "unknown"
    loc_conf = dict_get(protein, "deeplocpro_confidence") or 0.0

    # ---- HMMER outputs ----
    hmm_hits = dict_get(protein, "hmmer_hits") or []
    hmm_hits_set = set(hmm_hits)

    pfam_pos_list = set(params.get("hmmer_pfam_gram_pos", []))
    pfam_veto_list = set(params.get("hmmer_pfam_veto_pos", []))

    pfam_pos_hits = sorted(h for h in hmm_hits_set if h in pfam_pos_list)
    pfam_veto_hits = sorted(h for h in hmm_hits_set if h in pfam_veto_list)

    # ---- Tunable thresholds ----
    outside_loop_min = params.get("tmbed_outside_loop_min", 50)
    outside_frac_thr = params.get("tmbed_outside_frac_thr", 0.30)
    secreted_outside_frac = params.get("tmbed_secreted_outside_frac", 0.80)

    details = []
    category = "UNKNOWN"

    # -----------------------------------------------------------------
    # Annotation collation (build details list)
    # -----------------------------------------------------------------
    # SignalP
    if is_signalp:
        sp = signalp_type or "SP"
        cs = signalp_cut
        pr = signalp_pr
        if cs:
            details.append(f"signalp({sp};CS={cs};Pr={pr:.2f})")
            # If at some point you want to adjust loop coordinates, uncomment:
            # chop_nterminal_peptide(protein, cs)
        else:
            details.append(f"signalp({sp})")

    # TMbed aggregated detail (single "tmbed(...)" field)
    tmbed_items = []
    if n_helices:
        tmbed_items.append(f"H={n_helices}")
    if n_strands:
        tmbed_items.append(f"B={n_strands}")
    if longest_outside:
        tmbed_items.append(f"LoopOut={longest_outside}")
    if outside_fraction:
        tmbed_items.append(f"OutFrac={outside_fraction:.2f}")
    if tmbed_class:
        tmbed_items.append(f"Class={tmbed_class}")
    if tmbed_items:
        details.append("tmbed(" + ";".join(tmbed_items) + ")")

    # DeepLocPro detail
    if loc_label != "unknown":
        if loc_conf:
            details.append(f"deeplocpro({loc_label};Pr={loc_conf:.2f})")
        else:
            details.append(f"deeplocpro({loc_label})")

    # HMMER evidence
    hmm_items = []
    if hmm_hits:
        hmm_items.append(f"#Hits={len(hmm_hits)}")
    if pfam_pos_hits:
        hmm_items.append("PositiveHits=" + "+".join(pfam_pos_hits))
    if pfam_veto_hits:
        hmm_items.append("VetoedHits=" + "+".join(pfam_veto_hits))
        
    if hmm_items:
        details.append("hmmer(" + ";".join(hmm_items) + ")")

    # -----------------------------------------------------------------
    # Primary HMMER-driven classification (Gram+ anchors)
    # -----------------------------------------------------------------
    # If we see Gram+ receptor/anchor Pfams and no veto Pfams, let that
    # drive the main PSE class, then refine with SignalP/TMbed/DeepLocPro.
    if pfam_pos_hits and not pfam_veto_hits:
        # Lipoprotein anchors
        if signalp_is_lipop:
            category = "PSE-Lipoprotein"
        # Clear cell-wall-ish localization
        elif loc_label in ("cell_wall", "cell_wall_&_surface"):
            category = "PSE-Cellwall"
        # Have TM segments → membrane-exposed
        elif has_tm:
            category = "PSE-Membrane"
        # No TM segments → anchored to cell wall or surface
        else:
            category = "PSE-Cellwall"

    else:
        # -----------------------------------------------------------------
        # Structural / localization-based classification (no strong HMMER+)
        # -----------------------------------------------------------------

        # 1) Lipoproteins (SPII / Tat-SPII signals)
        if signalp_is_lipop:
            # Use TMbed + DeepLocPro to catch "SPII but actually secreted enzyme" cases
            secreted_like_lipo = (
                tmbed_class == "signal_only"
                and not has_tm
                and loc_label in ("extracellular", "cell_wall", "cell_wall_&_surface")
            )
            if secreted_like_lipo:
                category = "SECRETED"
            else:
                category = "PSE-Lipoprotein"

        # 2) Signal-containing / secretion-like proteins (non-lipoprotein)
        elif (
            is_signalp
            or has_tmbed_signal
            or loc_label in ["extracellular", "cell_wall", "cell_wall_&_surface"]
        ):
            if loc_label in ["cell_wall", "cell_wall_&_surface"]:
                category = "PSE-Cellwall"
            else:
                # Distinguish SECRETED vs PSE-Membrane by TM / outside evidence
                if not has_tm and outside_fraction >= secreted_outside_frac:
                    # Strongly outside, no TM → classic secreted
                    category = "SECRETED"
                elif has_tm:
                    exposed_like = (
                        longest_outside >= outside_loop_min
                    ) or (outside_fraction >= outside_frac_thr)
                    category = "PSE-Membrane" if exposed_like else "MEMBRANE(non-PSE)"
                else:
                    # Weak / ambiguous TMbed evidence → default to SECRETED
                    category = "SECRETED"

        # 3) Multi-pass membrane proteins (no explicit signal peptide)
        elif has_tm:
            exposed_like = (
                longest_outside >= outside_loop_min
            ) or (outside_fraction >= outside_frac_thr)
            category = "PSE-Membrane" if exposed_like else "MEMBRANE(non-PSE)"

        # 4) Localization-based fallbacks (no clear TMbed / signal info)
        elif loc_label in ["extracellular"]:
            category = "SECRETED"
        elif loc_label in ["cell_wall", "cell_wall_&_surface"]:
            category = "PSE-Cellwall"
        elif loc_label == "cytoplasmic_membrane":
            category = "MEMBRANE(non-PSE)"
        elif loc_label == "cytoplasmic":
            category = "CYTOPLASMIC"

        # 5) Default
        else:
            category = "CYTOPLASMIC"

    # -----------------------------------------------------------------
    # HMMER-based veto (demote obvious non-receptor Pfams)
    # -----------------------------------------------------------------
    # Only demote when *no* positive Pfams are present; if both appear,
    # we let positive Pfams dominate (you can tighten this later if needed).
    if pfam_veto_hits and not pfam_pos_hits:
        if category in {"PSE-Cellwall", "PSE-Membrane", "PSE-Lipoprotein"}:
            # Push into non-PSE versions based on structure / localization
            if has_tm:
                category = "MEMBRANE(non-PSE)"
            elif loc_label in ("extracellular", "cell_wall", "cell_wall_&_surface"):
                category = "SECRETED"
            elif loc_label == "cytoplasmic":
                category = "CYTOPLASMIC"
            else:
                category = "CYTOPLASMIC"

    # -----------------------------------------------------------------
    # PSE flag (phage receptor candidates)
    # -----------------------------------------------------------------
    # Only anchored / surface-accessible classes are PSE:
    pse_positive = {
        "PSE-Cellwall",
        "PSE-Membrane",
        "PSE-Lipoprotein",
    }

    # Extra condition: require a reasonably long outside loop
    outside_loop_min_for_pse = params.get("tmbed_outside_loop_min", 30)

    protein["category"] = category
    protein["is_pse"] = (category in pse_positive) and (longest_outside >= outside_loop_min_for_pse)
    details.append(f"PSE={protein['is_pse']}")

    # For backwards compatibility: store loop_extent as TMbed outside_fraction
    protein["details"] = details
    protein["loop_extent"] = outside_fraction

    return details, category


# ---------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------
def protein_output_line(seqid, proteins):
    p = proteins[seqid]
    return (
        f"{seqid:<15} {p['category']:<20}  "
        f"Exp={p.get('loop_extent', 0):.2f}  "
        f"{';'.join(p.get('details', [])):<60} "
        f"{p.get('name', '')[:60]}"
    )


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
