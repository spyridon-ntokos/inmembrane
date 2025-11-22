"""
Modernized Gram-negative protocol for SerraPHIM–inmembrane
-----------------------------------------------------------

Reimplementation of the original inmembrane Gram-negative surfaceome
analysis pipeline using updated predictors:

    - SignalP 6.0   (signal peptides + cleavage sites)
    - TMbed         (α-helical & β-barrel TM segments, signal segments, i/o loops)
    - DeepLocPro    (subcellular localization, Gram-negative mode)
    - HMMER/pyHMMER (Pfam-based OM receptor Pfams + custom HMM panels + veto Pfams)

Main output classes (category):

    - OM(barrel)        : classical OM β-barrel proteins
    - OM                : OM proteins inferred from HMMER when TMbed/DeepLoc are ambiguous
    - LIPOPROTEIN(OM)   : OM-anchored lipoproteins
    - LIPOPROTEIN(IM)   : IM-anchored lipoproteins
    - IM(peri+cyto)     : IM proteins with sizeable periplasmic/cytoplasmic loops
    - IM                : generic inner-membrane proteins
    - PERIPLASMIC
    - SECRETED
    - CYTOPLASMIC

Additionally sets:

    - is_pse (bool): whether the protein is considered “phage surface-exposed”
"""

from inmembrane.helpers import dict_get, chop_nterminal_peptide


# ---------------------------------------------------------------------
# Adjustable parameters (can be overridden in inmembrane.config)
# ---------------------------------------------------------------------
DEFAULT_PARAMS = {
    # SignalP 6
    "signalp_organism": "other",

    # TMbed heuristics
    "tmbed_barrel_min_strands": 8,        # minimal β-strands to call barrel
    "tmbed_outside_loop_min_neg": 14,     # min outside loop length → PSE-capable
    "tmbed_outside_fraction_exp": 0.25,   # outside fraction to support IM(peri+cyto)

    # DeepLocPro
    "deeplocpro_group": "negative",
    "deeplocpro_min_conf": 0.60,          # ignore DeepLocPro labels below this level

    # Optional: stricter HMM-rescue behaviour
    "tmbed_custom_anchor_min_loop": 14,   # min loop length to rescue via custom HMM
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
    upgrades (receptor Pfams/custom panels) and vetoes (non-receptor Pfams).

    Phage receptor candidacy is handled separately by
    _assign_phage_receptor_candidate_gram_neg(), which sets:

        protein["phage_receptor_candidate"]  (canonical)
        protein["is_pse"]                    (backwards-compatible alias)
        and appends 'PRC=<bool>' plus any PRC_via_* tags to details.
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
    # hmmer_hits      : positive DBs (Pfam + custom)
    # hmmer_veto_hits : veto DB (Pfam veto panel)
    hmm_hits = dict_get(protein, "hmmer_hits") or []
    hmm_veto_hits = dict_get(protein, "hmmer_veto_hits") or []

    hmm_hits_set = set(hmm_hits)
    hmm_veto_set = set(hmm_veto_hits)

    pfam_pos_list = set(params.get("hmmer_pfam_gram_neg", []))
    pfam_veto_list = set(params.get("hmmer_pfam_veto_neg", []))

    # Pfam positives
    pfam_pos_hits = sorted(h for h in hmm_hits_set if h in pfam_pos_list)
    # Pfam veto hits
    pfam_veto_hits = sorted(h for h in hmm_veto_set if h in pfam_veto_list)
    # Custom positive models (from custom HMM DBs)
    custom_hits = sorted(
        h
        for h in hmm_hits_set
        if (h not in pfam_pos_list and h not in pfam_veto_list)
    )

    # ----------------------------
    # Heuristic parameters
    # ----------------------------
    barrel_min_strands = params.get("tmbed_barrel_min_strands", 8)
    outside_loop_min = params.get("tmbed_outside_loop_min_neg", 14)
    outside_frac_thr = params.get("tmbed_outside_fraction_exp", 0.25)

    details = []
    category = "UNKNOWN"

    # ============================================================
    # Annotation collation → human-readable "Details" column
    # ============================================================

    # Signal peptides
    if is_signalp:
        sp = signalp_type or "SP"
        if signalp_is_lipop:
            sp += "~lipoprotein"
        cs = signalp_cut
        pr = signalp_pr
        if cs:
            details.append(f"signalp({sp};CS={cs};Pr={pr:.2f})")
            # If you ever want to adjust topology explicitly:
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

    # DeepLocPro (always show raw label, but category logic uses loc_label)
    if raw_loc_label != "unknown":
        if loc_conf:
            details.append(f"deeplocpro({raw_loc_label};Pr={loc_conf:.2f})")
        else:
            details.append(f"deeplocpro({raw_loc_label})")

    # HMMER evidence (Pfam + custom)
    hmm_items = []
    if hmm_hits:
        hmm_items.append(f"#Hits={len(hmm_hits)}")
    if pfam_pos_hits:
        hmm_items.append("PositiveHits=" + "+".join(pfam_pos_hits))
    if custom_hits:
        # models from custom positive DBs (e.g. CapsuleFinder HMMs)
        hmm_items.append("CustomHits=" + "+".join(custom_hits))
    if pfam_veto_hits:
        hmm_items.append("VetoedHits=" + "+".join(pfam_veto_hits))

    if hmm_items:
        details.append("hmmer(" + ";".join(hmm_items) + ")")

    # -----------------------------------------------------------------
    # Base classification logic (TMbed + SignalP + DeepLocPro)
    # -----------------------------------------------------------------

    # 1) Outer-membrane β-barrel proteins (TMbed + DeepLocPro)
    is_barrel = (
        is_barrel_candidate and n_strands >= barrel_min_strands
    ) or (
        loc_label == "outer_membrane" and n_strands >= (barrel_min_strands // 2)
    )

    if is_barrel:
        category = "OM(barrel)"

    # 2) Lipoproteins (SPII / Tat-SPII)
    elif signalp_is_lipop:
        # Likely periplasmic lipoprotein / secreted enzyme:
        periplip_like = (
            loc_label == "periplasmic"
            or (tmbed_class == "signal_only" and not tmbed_has_tm)
        )

        if periplip_like:
            category = "PERIPLASMIC"
        elif loc_label == "outer_membrane":
            category = "LIPOPROTEIN(OM)"
        elif loc_label == "cytoplasmic_membrane":
            category = "LIPOPROTEIN(IM)"
        else:
            # Fallback: assume OM-anchored if no strong periplasmic evidence
            category = "LIPOPROTEIN(OM)"

    # 3) α-helical IM proteins
    elif n_helices > 0:
        exposed_like = (
            longest_outside >= outside_loop_min
            # or outside_fraction >= outside_frac_thr
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
            category = "OM(barrel)"
        else:
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
    # HMMER-based refinements (Pfam + custom) on the *category*
    # -----------------------------------------------------------------

    # (A) Upgrade ambiguous non-OM classes if we see strong OM receptor Pfams
    #     or trusted custom HMMs, but only when structure/localization supports OM.
    if (pfam_pos_hits or custom_hits) and not pfam_veto_hits:
        if category not in {"OM(barrel)", "LIPOPROTEIN(OM)"}:
            upgrade_to_om = False

            # Pfam OM anchors (porins, TonB-dependent β-barrels, etc.)
            if pfam_pos_hits:
                if (n_strands > 0) or (loc_label == "outer_membrane"):
                    upgrade_to_om = True

            # Custom anchors (e.g. CapsuleFinder) – stricter:
            # require β-strand barrel-like evidence or strong OM localization.
            if custom_hits and not upgrade_to_om:
                if (n_strands >= barrel_min_strands // 2) or (loc_label == "outer_membrane"):
                    upgrade_to_om = True

            if upgrade_to_om:
                category = "OM"

    # (B) Veto obvious non-receptor cases: veto Pfam, no positive Pfam/custom
    if pfam_veto_hits and not pfam_pos_hits and not custom_hits:
        if category in {"OM(barrel)", "LIPOPROTEIN(OM)", "OM"}:
            if loc_label == "periplasmic":
                category = "PERIPLASMIC"
            elif loc_label == "cytoplasmic":
                category = "CYTOPLASMIC"
            else:
                category = "CYTOPLASMIC"

    # -----------------------------------------------------------------
    # PhReD / Bakta annotation evidence (tags only; PRC handled below)
    # -----------------------------------------------------------------
    name_lower = (dict_get(protein, "name") or "").lower()

    phred_hits = _phred_keyword_hits(
        name_lower,
        gram_is_negative=True,
    )

    for tag in phred_hits:
        details.append(f"PhReD={tag}")

    # -----------------------------------------------------------------
    # Store category + details and then assign phage receptor candidacy
    # -----------------------------------------------------------------
    protein["category"] = category
    protein["details"] = details
    protein["loop_extent"] = longest_outside

    _assign_phage_receptor_candidate_gram_neg(
        params=params,
        protein=protein,
        category=category,
        longest_outside=longest_outside,
        loc_label=loc_label,
        loc_conf=loc_conf,
        pfam_pos_hits=pfam_pos_hits,
        pfam_veto_hits=pfam_veto_hits,
        custom_hits=custom_hits,
        phred_hits=phred_hits,
    )

    return protein["details"], category


def _assign_phage_receptor_candidate_gram_neg(
    params,
    protein,
    category,
    longest_outside,
    loc_label,
    loc_conf,
    pfam_pos_hits,
    pfam_veto_hits,
    custom_hits,
    phred_hits,
):
    """
    Decide whether a Gram-negative protein is a phage receptor candidate.

    Sets:
        protein["phage_receptor_candidate"]
        protein["is_pse"]
        and appends PRC-related tags to protein["details"].
    """
    details = protein.get("details", [])

    # Categories that can in principle act as phage receptors
    prc_positive_categories = {
        "OM(barrel)",
        "LIPOPROTEIN(OM)",
        "OM",
    }

    # 1) Default structural requirement: sufficiently long outside loop
    outside_loop_min_for_prc = params.get(
        "tmbed_outside_loop_min_neg",
        params.get("tmbed_outside_loop_min", 14),
    )
    has_structural_exposure = (longest_outside >= outside_loop_min_for_prc)

    # 2) HMM-based overrides (Pfam or custom OM anchors)
    hmm_anchor_pfam = bool(pfam_pos_hits)
    hmm_anchor_custom = bool(custom_hits)

    strongly_cytoplasmic = (loc_label == "cytoplasmic" and loc_conf >= 0.80)

    # OM-like:
    om_like_category = category in {"OM(barrel)", "LIPOPROTEIN(OM)", "OM"}
    om_like_loc = (loc_label == "outer_membrane")

    custom_min_loop = params.get("tmbed_custom_anchor_min_loop", 14)

    # Custom HMM override:
    custom_anchor_exposure = (
        hmm_anchor_custom
        and (om_like_category or om_like_loc)
        and (longest_outside >= custom_min_loop)
        and not strongly_cytoplasmic
    )

    # Pfam anchor override (trusted OM Pfams)
    pfam_anchor_exposure = (
        hmm_anchor_pfam
        and om_like_category
        and not strongly_cytoplasmic
    )

    has_hmm_based_exposure = pfam_anchor_exposure or custom_anchor_exposure

    prc = (
        (category in prc_positive_categories)
        and (has_structural_exposure or has_hmm_based_exposure)
    )

    # annotate why it was accepted as PRC via HMM anchor
    if prc and has_hmm_based_exposure and not has_structural_exposure:
        details.append("PRC_via_HMM_anchor")

    # -----------------------------------------------------------------
    # Optional PhReD-driven upgrades (very conservative)
    # -----------------------------------------------------------------
    apply_phred_appendage = bool(params.get("phred_appendage_pse", False))

    def _has_tag(prefix: str) -> bool:
        return any(t.startswith(prefix) for t in phred_hits)

    name_lower = (dict_get(protein, "name") or "").lower()

    if apply_phred_appendage and not prc:
        # ---------- PILI / FIMBRIAE ----------
        if _has_tag("GN_appendages:pili"):
            looks_like_structural_pilin = (
                "pilin" in name_lower
                or "fimbrial" in name_lower
                or "fimbrillin" in name_lower
            )

            # Exclude assembly/stability/accessory proteins
            if looks_like_structural_pilin and not any(
                bad in name_lower for bad in ("assembly", "biogenesis", "stability")
            ):
                if longest_outside >= 50:
                    prc = True
                    details.append("PRC_via_PhReD_pilus")

        # ---------- FLAGELLA (very strict) ----------
        if _has_tag("GN_appendages:flagella") and not prc:
            # Only accept true flagellin
            if "flagellin" in name_lower or " flic" in name_lower:
                if longest_outside >= 100:
                    prc = True
                    details.append("PRC_via_PhReD_flagellin")

    # -----------------------------------------------------------------
    # Finalize fields
    # -----------------------------------------------------------------
    protein["phage_receptor_candidate"] = prc
    # backwards-compatible alias
    protein["is_pse"] = prc

    details.append(f"PRC={prc}")
    protein["details"] = details

def _phred_keyword_hits(name_lower, gram_is_negative=True):
    hits = []
    PHRED_KEYWORDS_GRAM_NEG_PROTEIN = {
        # Classical OM porins
        "OM_porins": [
            " ompa",          # "outer membrane protein A"
            " omp c", " ompc",
            " omp f", " ompf",
            " omp u", " ompu",
            " omp w", " ompw",
            " omplc", " omp lc",
            " ompk",
            " rop a0", " rop a1", " rop a2", " rop a3", " rop a4",
            " yaet", " bam a",   # YaeT/BamA
            " yncd",
            " ail ",            # Ail/Lom
        ],

        # LamB / Tsx
        "LamB_Tsx": [
            " lamb",            # LamB
            " tsx ",            # Tsx
        ],

        # TonB-dependent nutrient receptors
        "TonB_dep": [
            " fhu a", " ton a",     # FhuA/TonA
            " btu b",               # BtuB
            " fep a",               # FepA
            " tonb-dependent",      # general
            " ton b",               # the TonB protein
        ],

        # OM efflux factors and friends
        "OM_efflux": [
            " tolc",
            " oprm",
        ],

        # Other named OM receptor proteins
        "Other_OM": [
            " btu b",        # also in TonB_dep; tag both if matched
            " vcpq",
            " vp0980",
            " pertactin",
        ],
    }
    PHRED_KEYWORDS_GRAM_NEG_POLYSACCH = {
        "LPS_LipidA_core": [
            " lps ", " lipopolysaccharide", " lipooligosaccharide",
            " heptose ", " hepi", " hep ii", " hep iii",
            " kdo ", " ko residues",
            " lps core", " inner core", " outer core",
            " glucosamine residues of lps",
        ],

        "LPS_O_antigen": [
            " o antigen", " o-antigen", " o specific antigen",
            " o-specific antigen", " o polysaccharide",
            " o side chains", " o chain ", "osa lps",
        ],

        "capsule_cps_exopolysaccharide": [
            " capsule ", " capsular polysaccharide",
            " cps ", " exopolysaccharide",
            " extracellular polysaccharide",
            " colanic acid",
            " vi exopolysaccharide",
        ],

        # Generic “cell surface sugar” class
        "surface_glycans_generic": [
            " polysaccharide", " cell wall saccharides",
            " sugar moieties", " glycan", " glycans",
        ],
    }
    PHRED_KEYWORDS_GRAM_NEG_APPENDAGES = {
        "flagella": [
            " flagella", " flagellum", " flagellin", " flic", " fljb", " flaa",
        ],

        "pili": [
            " pilus", " pili", " pil a", " pilA",
            " type iv pilus", " type iv pili",
            " msha pilus", " msha type iv pilus",
            " conjugative pili", " rp4 plasmid",  # RP4-encoded pili
            " f pilus", " i pilus", " h pili", " j pili",
            " m pili", " n pili", " r pili", " t pili", " x pili",
        ],

        "fimbriae": [
            " fimbriae", " type 4 fimbriae",
        ],
    }
    PHRED_KEYWORDS_GRAM_POS_POLYSACCH = {
        "WTA": [
            " wall teichoic acid", " wta ", " teichoic acid",
            " glcna c residues of wall teichoic",
            " galactosylated glcnac", " glucosylated wall teichoic acid",
        ],

        "LTA": [
            " lipoteichoic acid", " lta ",
            " polyglycerophosphate teichoic acid",
            " poly(glycerophosphate)",  # lta backbone
        ],

        "CW_polysaccharide_capsule": [
            " cell wall polysaccharide",
            " capsular polysaccharide",
            " enterococcal polysaccharide antigen",
            " cell wall saccharides",
        ],
    }
    PHRED_KEYWORDS_GRAM_POS_PROTEIN = {
        "S_layer": [
            " s-layer", " s layer ", " slayer ",
            " sap ",          # Sap S-layer protein
        ],

        "phage_infection_proteins": [
            " pip ", " pipef", " yueb", " gamr", " yjae",
            " phage infection protein", " pip protein",
        ],

        # generic flagella
        "flagella": [
            " flagella", " flagellum", " flagellin", " flic",
        ],
    }

    if gram_is_negative:
        kw_sets = {
            "GN_OM_protein": PHRED_KEYWORDS_GRAM_NEG_PROTEIN,
            "GN_polysacch": PHRED_KEYWORDS_GRAM_NEG_POLYSACCH,
            "GN_appendages": PHRED_KEYWORDS_GRAM_NEG_APPENDAGES,
        }
    else:
        kw_sets = {
            "GP_polysacch": PHRED_KEYWORDS_GRAM_POS_POLYSACCH,
            "GP_protein": PHRED_KEYWORDS_GRAM_POS_PROTEIN,
        }

    for super_class, groups in kw_sets.items():
        for label, terms in groups.items():
            if any(t in name_lower for t in terms):
                hits.append(f"{super_class}:{label}")

    return hits

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
