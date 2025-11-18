"""
Modernized Gram-negative protocol for SerraPHIM–inmembrane
-----------------------------------------------------------

Reimplementation of the original inmembrane Gram-negative surfaceome
analysis pipeline using updated predictors:

    - SignalP 6.0  (signal peptides + cleavage sites)
    - TMbed        (α-helical & β-barrel TM segments, signal segments, i/o loops)
    - DeepLocPro   (subcellular localization, Gram-negative mode)
    - HMMER/pyHMMER (domain-level hints via Pfam lists)

Author: Spyridon Ntokos (SerraPHIM project)
"""

from inmembrane.helpers import dict_get, chop_nterminal_peptide


# ---------------------------------------------------------------------
# Adjustable parameters (can be overridden in inmembrane.config)
# ---------------------------------------------------------------------
DEFAULT_PARAMS = {
    # SignalP 6
    "signalp_organism": "other",

    # TMbed heuristics
    "tmbed_barrel_min_strands": 8,      # minimal β-strands to call barrel
    "tmbed_outside_loop_min": 30,       # min outside loop length → PSE-capable
    "tmbed_outside_fraction_exp": 0.25, # outside fraction to support IM(peri+cyto)

    # DeepLocPro
    "deeplocpro_group": "negative",
    "deeplocpro_min_conf": 0.60,        # ignore DeepLocPro labels below this
}


# ---------------------------------------------------------------------
# Protocol core
# ---------------------------------------------------------------------
def get_annotations(params):
    """
    Return predictor IDs required for the Gram-negative workflow and
    populate default tunable parameters.
    """
    predictors = ["signalp6", "tmbed", "deeplocpro", "hmmer"]
    params.setdefault("predictors", predictors)

    for k, v in DEFAULT_PARAMS.items():
        params.setdefault(k, v)

    return predictors


def post_process_protein(params, protein):
    """
    Integrate predictor outputs and classify each protein according to
    Gram-negative topology + localization logic, with HMMER-based
    upgrades (receptor Pfams) and vetoes (non-receptor Pfams).
    """

    # ----------------------------
    # SignalP 6.0 features
    # ----------------------------
    is_signalp = bool(dict_get(protein, "signalp_is_sp"))
    signalp_type = dict_get(protein, "signalp_type") or ""
    signalp_cut = dict_get(protein, "signalp_cleave_position") or 0
    signalp_pr = dict_get(protein, "signalp_cleave_prob") or 0.0

    # Lipoprotein if Sec/SPII or Tat/SPII
    signalp_is_lipop = any(signalp_type.endswith(x) for x in ("Sec/SPII", "Tat/SPII"))
    is_tat = "tat" in signalp_type.lower()

    # ----------------------------
    # TMbed features
    # ----------------------------
    n_helices = dict_get(protein, "tmbed_n_helices") or 0
    n_strands = dict_get(protein, "tmbed_n_strands") or 0
    is_barrel_candidate = bool(dict_get(protein, "tmbed_is_barrel_candidate"))
    outside_fraction = dict_get(protein, "tmbed_outside_fraction") or 0.0
    longest_outside = dict_get(protein, "tmbed_longest_outside_loop") or 0
    has_tmbed_signal = bool(dict_get(protein, "tmbed_has_signal"))
    tmbed_class = dict_get(protein, "tmbed_class") or None

    # ----------------------------
    # DeepLocPro features
    # ----------------------------
    raw_loc_label = dict_get(protein, "deeplocpro_label") or "unknown"
    loc_conf = dict_get(protein, "deeplocpro_confidence") or 0.0
    min_loc_conf = params.get("deeplocpro_min_conf", 0.6)

    # Only trust localization if confidence is high enough
    loc_label = raw_loc_label if loc_conf >= min_loc_conf else "unknown"

    # ----------------------------
    # HMMER / pyHMMER features
    # ----------------------------
    hmm_hits = dict_get(protein, "hmmer_hits") or []
    hmm_hits_set = set(hmm_hits)

    pfam_pos_list = set(params.get("hmmer_pfam_gram_neg", []))
    pfam_veto_list = set(params.get("hmmer_pfam_veto_neg", []))

    pfam_pos_hits = sorted(h for h in hmm_hits_set if h in pfam_pos_list)
    pfam_veto_hits = sorted(h for h in hmm_hits_set if h in pfam_veto_list)

    # ----------------------------
    # Heuristic parameters
    # ----------------------------
    barrel_min_strands = params.get("tmbed_barrel_min_strands", 8)
    outside_loop_min = params.get("tmbed_outside_loop_min", 30)
    outside_frac_thr = params.get("tmbed_outside_fraction_exp", 0.25)

    details = []
    category = "UNKNOWN"

    # ============================================================
    # Annotation collation → human-readable "Details" column
    # ============================================================

    # Signal peptides
    if is_signalp:
        sp = signalp_type or "SP"
        # Lipoprotein
        if signalp_is_lipop:
            sp += "~lipoprotein"
        cs = signalp_cut
        pr = signalp_pr
        if cs:
            details.append(f"signalp({sp};CS={cs};Pr={pr:.2f})")
            # chop_nterminal_peptide(protein, cs)
        else:
            details.append(f"signalp({sp})")
    elif is_tat:
        details.append("signalp(Tat-like)")

    # TMbed signal-only (when SignalP did not trigger)
    if (not is_signalp) and has_tmbed_signal:
        details.append("tmbed(SP)")

    # --- TMbed evidence aggregation ---
    tmbed_has_tm = False
    tmbed_items = []
    if n_helices:
        tmbed_items.append(f"H={n_helices}")
        tmbed_has_tm = True
    if n_strands:
        tmbed_items.append(f"B={n_strands}")
        tmbed_has_tm = True
    if longest_outside:
        tmbed_items.append(f"LoopOut={longest_outside}")
    if outside_fraction:
        tmbed_items.append(f"OutFrac={outside_fraction:.2f}")
    if tmbed_class:
        tmbed_items.append(f"Class={tmbed_class}")

    if tmbed_items:
        details.append("tmbed(" + ";".join(tmbed_items) + ")")

    # DeepLocPro
    if raw_loc_label != "unknown":
        if loc_conf:
            details.append(f"deeplocpro({raw_loc_label};Pr={loc_conf:.2f})")
        else:
            details.append(f"deeplocpro({raw_loc_label})")

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
    # Base classification logic (TMbed + SignalP + DeepLocPro)
    # -----------------------------------------------------------------

    # 1) Outer-membrane β-barrel proteins (TMbed + DeepLocPro)
    is_barrel = (
        is_barrel_candidate
        and n_strands >= barrel_min_strands
    ) or (
        loc_label == "outer_membrane"
        and n_strands >= (barrel_min_strands // 2)
    )

    if is_barrel:
        category = "OM(barrel)"

    # 2) Lipoproteins (SPII / Tat-SPII)
    elif signalp_is_lipop:
        # Likely periplasmic lipoprotein / secreted enzyme:
        # - DeepLocPro says periplasmic
        # - TMbed says "signal_only" and sees no TM segments
        periplip_like = (
            loc_label == "periplasmic"
            or (tmbed_class == "signal_only" and not tmbed_has_tm)
        )

        if periplip_like:
            # Treat as periplasmic, not surface-exposed
            category = "PERIPLASMIC"
        elif loc_label == "outer_membrane":
            category = "LIPOPROTEIN(OM)"
        elif loc_label == "cytoplasmic_membrane":
            category = "LIPOPROTEIN(IM)"
        else:
            # Fallback: still call OM, but only when there is no
            # strong periplasmic / signal_only evidence
            category = "LIPOPROTEIN(OM)"

    # 3) α-helical IM proteins
    elif n_helices > 0:
        exposed_like = (
            longest_outside >= outside_loop_min
            or outside_fraction >= outside_frac_thr
        )
        category = "IM(peri+cyto)" if exposed_like else "IM"

    # 4) Signal-containing proteins → PERIPLASMIC or SECRETED
    elif is_signalp or has_tmbed_signal:
        if n_helices == 0 and n_strands == 0:
            if loc_label == "extracellular" or outside_fraction >= 0.80:
                category = "SECRETED"
            else:
                category = "PERIPLASMIC"
        else:
            category = "PERIPLASMIC"

    # 5) DeepLocPro-only OM fallback:
    #    require at least one TMbed β-strand to accept OM(barrel)
    elif loc_label == "outer_membrane":
        if n_strands > 0:
            # Some β-strand evidence + OM label → accept as OM(barrel)
            category = "OM(barrel)"
        else:
            # No TMbed TM segments at all → treat as non-membrane
            category = "CYTOPLASMIC"

    # 6) Other DeepLocPro-based fallbacks
    elif loc_label == "periplasmic":
        category = "PERIPLASMIC"
    elif loc_label == "cytoplasmic_membrane":
        category = "IM"
    elif loc_label == "cytoplasmic":
        category = "CYTOPLASMIC"

    # 7) Default
    else:
        category = "CYTOPLASMIC"

    # -----------------------------------------------------------------
    # HMMER-based refinements (Pfam + / veto)
    # -----------------------------------------------------------------

    # (A) Upgrade ambiguous non-OM classes if we see strong OM receptor Pfams
    if pfam_pos_hits and not pfam_veto_hits:
        if category not in {"OM(barrel)", "LIPOPROTEIN(OM)"}:
            # Require at least some structural or localization support
            if (n_strands > 0) or tmbed_has_tm or (loc_label == "outer_membrane"):
                category = "OM"

    # (B) Veto obvious non-receptor cases: veto Pfam, no positive Pfam
    if pfam_veto_hits and not pfam_pos_hits:
        if category in {"OM(barrel)", "LIPOPROTEIN(OM)"}:
            # Collapse to non-PSE; default to cytosolic unless DeepLoc suggests otherwise
            if loc_label == "periplasmic":
                category = "PERIPLASMIC"
            elif loc_label == "cytoplasmic":
                category = "CYTOPLASMIC"
            else:
                category = "CYTOPLASMIC"

    # -----------------------------------------------------------------
    # PSE flag (phage receptor candidates)
    # -----------------------------------------------------------------
    pse_positive = {
        "OM(barrel)",
        "LIPOPROTEIN(OM)",
        "OM"
    }

    protein["category"] = category
    protein["is_pse"] = (category in pse_positive) and (longest_outside >= outside_loop_min)
    details.append(f"PSE={protein['is_pse']}")

    # ----------------------------
    # Store results
    # ----------------------------
    # For historical reasons, loop_extent was used for "exposed fraction".
    # Here we reuse it as the *longest outside loop length* from TMbed.
    protein["loop_extent"] = longest_outside
    protein["details"] = details

    return details, category


# ---------------------------------------------------------------------
# Output formatting helpers
# ---------------------------------------------------------------------
def protein_output_line(seqid, proteins):
    p = proteins[seqid]
    return (
        f"{seqid:<15} {p['category']:<20}  "
        f"LoopOut={p.get('loop_extent', 0):>3}  "
        f"{';'.join(p.get('details', [])):<60} {p.get('name', '')[:60]}"
    )


def protein_csv_line(seqid, proteins, header: bool = False):
    """
    Returns one CSV line per protein. Adds header if requested.
    NOTE: __init__.py currently calls this *without* header=True.
    """
    lines = []
    if header:
        lines.append("SeqID,Category,LoopExtent,Details,Name\n")
    p = proteins[seqid]
    lines.append(
        f"{seqid},{p['category']},{p.get('loop_extent', 0):.0f},"
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
