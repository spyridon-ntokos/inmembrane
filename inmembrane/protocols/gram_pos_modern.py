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

import re
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
    HMMER-based upgrades (anchor Pfams), vetoes (non-receptor Pfams),
    and optional annotation-based veto for phage structural / core IM.

    Phage receptor candidacy is handled separately by
    _assign_phage_receptor_candidate_gram_pos(), which sets:

        protein["phage_receptor_candidate"]  (canonical)
        protein["is_pse"]                    (backwards-compatible alias)
        and appends 'PRC=<bool>' plus any PRC_via_* tags to details.
    """

    # -----------------------------------------------------------------
    # Predictor outputs
    # -----------------------------------------------------------------
    # ---- SignalP 6.0 ----
    is_signalp = dict_get(protein, "signalp_is_sp")
    signalp_type = dict_get(protein, "signalp_type") or ""
    signalp_cut = dict_get(protein, "signalp_cleave_position") or 0
    signalp_pr = dict_get(protein, "signalp_cleave_prob") or 0.0
    signalp_is_lipop = any(
        signalp_type.endswith(x) for x in ("Sec/SPII", "Tat/SPII")
    )

    # ---- TMbed ----
    n_helices = dict_get(protein, "tmbed_n_helices") or 0
    n_strands = dict_get(protein, "tmbed_n_strands") or 0
    longest_outside = dict_get(protein, "tmbed_longest_outside_loop") or 0
    outside_fraction = dict_get(protein, "tmbed_outside_fraction") or 0.0
    has_tmbed_signal = dict_get(protein, "tmbed_has_signal") or False
    tmbed_class = dict_get(protein, "tmbed_class") or None
    # Only count beta strands as TM if there are enough for a stable element
    has_tm = (n_helices > 0) or (n_strands >= 4)

    # ---- DeepLocPro ----
    loc_label = dict_get(protein, "deeplocpro_label") or "unknown"
    loc_conf = dict_get(protein, "deeplocpro_confidence") or 0.0

    # ---- HMMER (Pfam + custom hits) ----
    hmm_hits = dict_get(protein, "hmmer_hits") or []           # positive DBs (Pfam + custom)
    hmm_veto_hits = dict_get(protein, "hmmer_veto_hits") or [] # veto DB (Pfam veto panel)
    hmm_hits_set = set(hmm_hits)
    hmm_veto_set = set(hmm_veto_hits)

    pfam_pos_list = set(params.get("hmmer_pfam_gram_pos", []))
    pfam_veto_list = set(params.get("hmmer_pfam_veto_pos", []))

    # Pfam positives (from positive DB)
    pfam_pos_hits = sorted(h for h in hmm_hits_set if h in pfam_pos_list)
    # Pfam veto (from veto DB)
    pfam_veto_hits = sorted(h for h in hmm_veto_set if h in pfam_veto_list)
    # Custom positives (present in hmmer_hits but not in Pfam panels)
    custom_hits = sorted(
        h
        for h in hmm_hits_set
        if (h not in pfam_pos_list and h not in pfam_veto_list)
    )

    # -----------------------------------------------------------------
    # Tunable thresholds
    # -----------------------------------------------------------------
    outside_loop_min = params.get("tmbed_outside_loop_min_pos", 50)
    outside_frac_thr = params.get("tmbed_outside_frac_thr", 0.30)
    secreted_outside_frac = params.get("tmbed_secreted_outside_frac", 0.80)

    # Whether to apply annotation-based veto to PRC flag (project-specific)
    apply_annotation_veto = bool(params.get("phage_receptor_veto", False))

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

    # TMbed
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

    # DeepLocPro
    if loc_label != "unknown":
        if loc_conf:
            details.append(f"deeplocpro({loc_label};Pr={loc_conf:.2f})")
        else:
            details.append(f"deeplocpro({loc_label})")

    # HMMER summary (Pfam + custom)
    hmm_items = []
    if hmm_hits:
        hmm_items.append(f"#Hits={len(hmm_hits)}")
    if pfam_pos_hits:
        hmm_items.append("PositiveHits=" + "+".join(pfam_pos_hits))
    if custom_hits:
        # These are models coming from custom positive DBs (e.g. CapsuleFinder)
        # 'h' is whatever we stored as hit["pfam"] → Pfam accession or model NAME.
        hmm_items.append("CustomHits=" + "+".join(custom_hits))
    if pfam_veto_hits:
        hmm_items.append("VetoedHits=" + "+".join(pfam_veto_hits))

    if hmm_items:
        details.append("hmmer(" + ";".join(hmm_items) + ")")

    # -----------------------------------------------------------------
    # Primary HMMER-driven classification (Gram+ anchors)
    # -----------------------------------------------------------------
    if pfam_pos_hits or custom_hits and not pfam_veto_hits:
        # Lipoprotein anchors
        if signalp_is_lipop:
            category = "PSE-Lipoprotein"
        # Clear cell-wall-ish localization
        elif loc_label in ("cell_wall", "cell_wall_&_surface"):
            category = "PSE-Cellwall"
        # TM segments → exposed membrane protein
        elif has_tm:
            category = "PSE-Membrane"
        # No TM segments → likely cell-wall/surface anchored
        else:
            category = "PSE-Cellwall"

    else:
        # -----------------------------------------------------------------
        # Structural / localization-based classification (no strong HMM+)
        # -----------------------------------------------------------------

        # 1) Lipoproteins (SPII / Tat-SPII)
        if signalp_is_lipop:
            # Catch "lipoprotein but effectively secreted" cases
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
            or loc_label in ("extracellular", "cell_wall", "cell_wall_&_surface")
        ):
            if loc_label in ("cell_wall", "cell_wall_&_surface"):
                category = "PSE-Cellwall"
            else:
                # Distinguish SECRETED vs PSE-Membrane by TM / outside evidence
                if not has_tm and outside_fraction >= secreted_outside_frac:
                    # Strongly outside, no TM → classic secreted
                    category = "SECRETED"
                elif has_tm:
                    exposed_like = (longest_outside >= outside_loop_min)
                    # If outside_fraction should also count, uncomment:
                    # exposed_like = exposed_like or (outside_fraction >= outside_frac_thr)
                    category = "PSE-Membrane" if exposed_like else "MEMBRANE(non-PSE)"
                else:
                    # Weak / ambiguous TMbed evidence → default to SECRETED
                    category = "SECRETED"

        # 3) Multi-pass membrane proteins (no explicit signal peptide)
        elif has_tm:
            exposed_like = (longest_outside >= outside_loop_min)
            # If outside_fraction should also count, uncomment:
            # exposed_like = exposed_like or (outside_fraction >= outside_frac_thr)
            category = "PSE-Membrane" if exposed_like else "MEMBRANE(non-PSE)"

        # 4) Localization-based fallbacks (no clear TMbed / signal info)
        elif loc_label == "extracellular":
            category = "SECRETED"
        elif loc_label in ("cell_wall", "cell_wall_&_surface"):
            category = "PSE-Cellwall"
        elif loc_label == "cytoplasmic_membrane":
            category = "MEMBRANE(non-PSE)"
        elif loc_label == "cytoplasmic":
            category = "CYTOPLASMIC"

        # 5) Default
        else:
            category = "CYTOPLASMIC"

    # -----------------------------------------------------------------
    # HMMER-based veto (non-receptor Pfams) on *category*
    # -----------------------------------------------------------------
    if pfam_veto_hits and not pfam_pos_hits and not custom_hits:
        if category in {"PSE-Cellwall", "PSE-Membrane", "PSE-Lipoprotein"}:
            if has_tm:
                category = "MEMBRANE(non-PSE)"
            elif loc_label in ("extracellular", "cell_wall", "cell_wall_&_surface"):
                category = "SECRETED"
            elif loc_label == "cytoplasmic":
                category = "CYTOPLASMIC"
            else:
                category = "CYTOPLASMIC"

    # -----------------------------------------------------------------
    # Annotation-based veto for phage structural and housekeeping IM
    # (only affects phage receptor candidacy, not the main category)
    # -----------------------------------------------------------------
    name_lower = (dict_get(protein, "name") or "").lower()
    force_non_prc = False
    veto_tags = []

    if apply_annotation_veto and name_lower:
        # phage structural / lysis proteins
        phage_structural_keywords = [
            "tape measure protein",
            "tail protein",
            "tail fibre",
            "tail fiber",
            "tail component",
            "baseplate",
            "capsid",
            "portal protein",
            "major tail protein",
            "minor tail protein",
        ]

        if (
            "holin" in name_lower
            or any(kw in name_lower for kw in phage_structural_keywords)
            or ("phage" in name_lower and any(
                kw in name_lower for kw in ["tail", "capsid", "portal", "structural"]
            ))
        ):
            force_non_prc = True
            veto_tags.append("veto_phage_structural")

        # housekeeping IM / divisome / generic multi-pass transporters
        if not force_non_prc and category == "PSE-Membrane":
            housekeeping_im_phrases = [
                "cell division protein ftsw",
                "cell division protein ftsl",
                "cell division protein ftsi",
                "ftsw",
                "ftsl",
                "ftsi",
                "ftsq",
                "divib",
                "mrec",
                "mred",
                "rod shape-determining protein",
                "penicillin-binding protein",
                "secdf",
                "lipoteichoic acid synthase",
                "lipoteichoic acid-specific glycosyltransferase",
            ]

            if any(phrase in name_lower for phrase in housekeeping_im_phrases):
                force_non_prc = True
                veto_tags.append("veto_housekeeping_im")
            else:
                transporter_words = [
                    "permease",
                    "symporter",
                    "antiporter",
                    "transporter",
                    "efflux",
                ]
                if (
                    n_helices >= 7
                    and any(w in name_lower for w in transporter_words)
                    and "substrate-binding" not in name_lower
                    and "binding protein" not in name_lower
                    and "binding lipoprotein" not in name_lower
                    and "lipoprotein" not in name_lower
                ):
                    force_non_prc = True
                    veto_tags.append("veto_generic_im_transporter")

    for tag in veto_tags:
        details.append(tag)

    # -----------------------------------------------------------------
    # PhReD / Bakta annotation evidence (tags only; PRC handled below)
    # -----------------------------------------------------------------
    phred_hits = _phred_keyword_hits(
        name_lower,
        gram_is_negative=False,
    )

    for tag in phred_hits:
        details.append(f"PhReD={tag}")

    # -----------------------------------------------------------------
    # Store category + details and then assign phage receptor candidacy
    # -----------------------------------------------------------------
    protein["category"] = category
    protein["details"] = details  # helper will extend this list
    protein["loop_extent"] = longest_outside

    _assign_phage_receptor_candidate_gram_pos(
        params=params,
        protein=protein,
        category=category,
        longest_outside=longest_outside,
        loc_label=loc_label,
        loc_conf=loc_conf,
        pfam_pos_hits=pfam_pos_hits,
        custom_hits=custom_hits,
        force_non_prc=force_non_prc,
        n_helices=n_helices,
        n_strands=n_strands,
        has_tm=has_tm,
        name_lower=name_lower,
    )

    return protein["details"], category

def _assign_phage_receptor_candidate_gram_pos(
    params,
    protein,
    category,
    longest_outside,
    loc_label,
    loc_conf,
    pfam_pos_hits,
    custom_hits,
    force_non_prc,
    n_helices,
    n_strands,
    has_tm,
    name_lower,
):
    """
    Decide whether a Gram-positive protein is a phage receptor candidate.

    Sets:
        protein["phage_receptor_candidate"]
        protein["is_pse"]              (backwards-compatible alias)
        and appends PRC-related tags to protein["details"].
    """
    details = protein.get("details", [])

    prc_positive_categories = {
        "PSE-Cellwall",
        "PSE-Membrane",
        "PSE-Lipoprotein",
    }

    outside_loop_min_for_prc = params.get("tmbed_outside_loop_min", 30)

    # 1) Default structural requirement: sufficiently long outside loop
    has_structural_exposure = (longest_outside >= outside_loop_min_for_prc)

    # 2) HMM-based override – conservative
    hmm_anchor_pfam = bool(pfam_pos_hits)
    hmm_anchor_custom = bool(custom_hits)

    # Strongly cytoplasmic according to DeepLocPro → do NOT override with HMMs
    strongly_cytoplasmic = (loc_label == "cytoplasmic" and loc_conf >= 0.80)

    # “Sensor-like” IM architecture: few helices, no beta-barrel
    few_helices = (n_helices > 0 and n_helices <= 3 and n_strands == 0)

    # words that scream “core transporter / pump”
    transporter_words = ["transporter", "permease", "efflux"]

    is_likely_core_transporter = (
        has_tm
        and n_helices >= 5
        and any(w in name_lower for w in transporter_words)
    )

    # Minimum loop length required to rescue PRC via custom HMMs
    custom_min_loop = params.get("tmbed_custom_anchor_min_loop", 20)

    # Custom HMM override is allowed only if:
    #   - we have custom hits,
    #   - protein is clearly an IM protein (TM + CM loc),
    #   - it is NOT a generic core transporter,
    #   - it has at least a modest outside loop,
    #   - and is not strongly cytoplasmic.
    custom_anchor_exposure = (
        hmm_anchor_custom
        and has_tm
        and loc_label == "cytoplasmic_membrane"
        and loc_conf >= 0.70
        and not strongly_cytoplasmic
        and not is_likely_core_transporter
        and longest_outside >= custom_min_loop
    )

    # Pfam anchor override (trusted Gram+ Pfams; avoid big transporters)
    pfam_anchor_exposure = (
        hmm_anchor_pfam
        and few_helices
        and not strongly_cytoplasmic
    )

    has_hmm_based_exposure = pfam_anchor_exposure or custom_anchor_exposure

    prc = (
        (category in prc_positive_categories)
        and (has_structural_exposure or has_hmm_based_exposure)
        and not force_non_prc
    )

    # annotate why it was accepted as PRC via HMM anchor
    if prc and has_hmm_based_exposure and not has_structural_exposure:
        details.append("PRC_via_HMM_anchor")

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
# Output formatting
# ---------------------------------------------------------------------
def protein_output_line(seqid, proteins):
    p = proteins[seqid]
    return (
        f"{seqid:<15} {p['category']:<20}  "
        f"Exp={p.get('loop_extent', 0):>3}  "
        f"{';'.join(p.get('details', [])):<60} "
        f"{p.get('name', '')[:60]}"
    )


def protein_csv_line(seqid, proteins, header=False):
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

    out = "\n\n# Number of proteins in each class:\n"
    for c, n in sorted(counts.items()):
        out += f"# {c:<25}\t{n}\n"
    return out
