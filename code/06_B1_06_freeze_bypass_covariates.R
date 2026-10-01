# ============================================================
# A7 / B-line
# 06_B1_06_freeze_bypass_covariates.R
#
# PURPOSE
#   Outcome-blind freezing of the genomic exposure definition and
#   candidate covariates for the anti-EGFR biomarker study.
#
# PRIMARY STUDY COHORT
#   - OLD_ELIGIBLE_803_R == 1
#   - FIRST_BIOLOGIC_R == "anti-EGFR"
#   - sequencing completed on/before first anti-EGFR exposure
#
# IMPORTANT
#   This script does NOT read progression / rwPFS / OS / death.
#   It does NOT run KM, Cox, HR, p-values, or outcome comparisons.
#   It does NOT modify the two formal *_R.csv files.
#
# FROZEN EXPOSURE TIERS
#
# Tier 1: PRESSING-like strict core
#   - ERBB2 amplification
#   - high-confidence activating ERBB2 mutation
#   - MET amplification
#   - PIK3CA exon-20 / kinase-domain hotspot
#   - clear PTEN loss-of-function
#   - AKT1 E17K
#   - canonical ALK/ROS1/NTRK1-3/RET fusion
#
# Tier 2: PRIMARY high-confidence EGFR-bypass panel
#   Tier 1 PLUS:
#   - KRAS amplification
#   - MAP2K1 resistance/activating hotspot
#   - clear NF1 loss-of-function
#
# Tier 3: Original broad composite
#   Existing BYPASS_ORIGINAL_R; sensitivity only.
#
# NON-V600 BRAF
#   Class III BRAF variants are NOT classified as bypass.
#   They are flagged for a pre-specified exclusion sensitivity.
#
# OUTPUTS
#   06_logs_and_audit/B1_06_bypass_covariate_freeze_audit.txt
#   06_logs_and_audit/B1_06_variant_inventory.csv
#   06_logs_and_audit/B1_06_patient_genomic_covariate_audit.csv
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory
raw_dir <- file.path(project_root, "01_raw_data")
audit_dir <- file.path(project_root, "06_logs_and_audit")

main_file <- file.path(
    raw_dir,
    "A7_MSK_CHORD_CRC_clean_R.csv"
)

tar_file <- file.path(
    raw_dir,
    "msk_chord_2024.tar"
)

if (!file.exists(main_file)) {
    stop("Missing main R-ready CSV:\n", main_file)
}

if (!file.exists(tar_file)) {
    stop("Missing MSK-CHORD tar:\n", tar_file)
}

if (!dir.exists(audit_dir)) {
    dir.create(
        audit_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
}

if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}
library(data.table)

cat("\n============================================================\n")
cat("B1-06 OUTCOME-BLIND BYPASS / COVARIATE FREEZE\n")
cat("============================================================\n\n")


# ============================================================
# Helpers
# ============================================================

min_or_na <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
        return(NA_real_)
    }

    min(x)
}

pct <- function(n, d) {

    if (is.na(d) || d == 0L) {
        return("NA")
    }

    sprintf(
        "%d/%d (%.1f%%)",
        n,
        d,
        100 * n / d
    )
}

median_iqr <- function(x) {

    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
        return("NA")
    }

    q <- quantile(
        x,
        probs = c(0.25, 0.50, 0.75),
        na.rm = TRUE,
        names = FALSE
    )

    sprintf(
        "%.1f [%.1f, %.1f]",
        q[2],
        q[1],
        q[3]
    )
}

safe_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
}

collapse_unique <- function(x) {

    x <- safe_chr(x)
    x <- trimws(x)
    x <- x[nzchar(x)]
    x <- unique(x)

    if (length(x) == 0L) {
        return(NA_character_)
    }

    paste(
        sort(x),
        collapse = "|"
    )
}

extract_site_code <- function(x) {

    x <- toupper(
        safe_chr(x)
    )

    out <- rep(
        NA_character_,
        length(x)
    )

    hit <- grepl(
        "C[0-9]{3}",
        x,
        perl = TRUE
    )

    if (any(hit)) {

        temp <- regmatches(
            x[hit],
            gregexpr(
                "C[0-9]{3}",
                x[hit],
                perl = TRUE
            )
        )

        out[hit] <- vapply(
            temp,
            function(z) {

                if (length(z) == 0L) {
                    NA_character_
                } else {
                    tail(z, 1L)
                }
            },
            character(1)
        )
    }

    out
}


# ============================================================
# 1. Primary outcome-blind cohort from the formal R-ready CSV
# ============================================================

main <- fread(
    main_file,
    check.names = FALSE,
    showProgress = TRUE
)

required_main <- c(
    "PATIENT_ID",
    "SAMPLE_ID",
    "OLD_ELIGIBLE_803_R",
    "FIRST_BIOLOGIC_R",
    "FIRST_BIOLOGIC_DAY_R",
    "BYPASS_ORIGINAL_R",
    "SEQUENCING_FIRST_DAY_R",
    "SPECIMEN_ACQUISITION_FIRST_DAY_R",
    "GENDER",
    "PRIMARY_SITE"
)

missing_main <- setdiff(
    required_main,
    names(main)
)

if (length(missing_main) > 0L) {

    stop(
        "Main CSV missing required columns:\n",
        paste(
            missing_main,
            collapse = "\n"
        )
    )
}

cohort <- main[
    OLD_ELIGIBLE_803_R == 1L &
        FIRST_BIOLOGIC_R == "anti-EGFR" &
        is.finite(SEQUENCING_FIRST_DAY_R) &
        is.finite(FIRST_BIOLOGIC_DAY_R) &
        SEQUENCING_FIRST_DAY_R <= FIRST_BIOLOGIC_DAY_R
]

if (nrow(cohort) != 191L) {

    stop(
        "Primary clean anti-EGFR anchor failed. ",
        "Expected 191; got ",
        nrow(cohort)
    )
}

cohort <- cohort[
    ,
    .(
        PATIENT_ID,
        SAMPLE_ID,
        FIRST_BIOLOGIC_DAY_R,
        BYPASS_ORIGINAL_R,
        SEQUENCING_FIRST_DAY_R,
        SPECIMEN_ACQUISITION_FIRST_DAY_R,
        GENDER,
        PRIMARY_SITE
    )
]

cat(
    "Primary clean anti-EGFR cohort:",
    nrow(cohort),
    "\n"
)


# ============================================================
# 2. Temporary extraction of original MSK-CHORD files
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_06_"
)

dir.create(
    tmp_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

on.exit(
    unlink(
        tmp_dir,
        recursive = TRUE,
        force = TRUE
    ),
    add = TRUE
)

untar(
    tar_file,
    exdir = tmp_dir
)

all_files <- list.files(
    tmp_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
)

find_raw <- function(filename) {

    hit <- all_files[
        basename(all_files) == filename
    ]

    if (length(hit) == 0L) {
        stop(
            "Missing raw MSK file: ",
            filename
        )
    }

    hit[1]
}


# ============================================================
# 3. Sample-level assay / specimen metadata
# ============================================================

sample <- fread(
    find_raw(
        "data_clinical_sample.txt"
    ),
    sep = "\t",
    skip = 4,
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

sample_primary <- sample[
    SAMPLE_ID %chin% cohort$SAMPLE_ID,
    .(
        SAMPLE_ID,
        PATIENT_ID,
        GENE_PANEL,
        SAMPLE_TYPE,
        METASTATIC_SITE,
        PRIMARY_SITE_RAW = PRIMARY_SITE,
        TMB_NONSYNONYMOUS
    )
]

if (nrow(sample_primary) != 191L) {

    stop(
        "Sample map failed. Expected 191 sample rows; got ",
        nrow(sample_primary)
    )
}

sample_map <- sample_primary[
    ,
    .(
        SAMPLE_ID,
        PATIENT_ID
    )
]


# ============================================================
# 4. Mutation inventory
# ============================================================

mutation_cols <- c(
    "Hugo_Symbol",
    "Tumor_Sample_Barcode",
    "HGVSp",
    "HGVSp_Short",
    "Variant_Classification"
)

mut <- fread(
    find_raw(
        "data_mutations.txt"
    ),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    select = mutation_cols,
    showProgress = TRUE
)

mut <- mut[
    Tumor_Sample_Barcode %chin%
        cohort$SAMPLE_ID
]

mut[
    ,
    PATIENT_ID := sample_map$PATIENT_ID[
        match(
            Tumor_Sample_Barcode,
            sample_map$SAMPLE_ID
        )
    ]
]

genes_inventory <- c(
    "ERBB2",
    "MET",
    "PIK3CA",
    "PTEN",
    "AKT1",
    "KRAS",
    "MAP2K1",
    "NF1",
    "ERBB3",
    "FGFR2",
    "IGF1R",
    "ARAF",
    "AKT2",
    "MAP2K2",
    "MAP2K4",
    "EGFR",
    "BRAF",
    "GNAS"
)

mut_goi <- mut[
    Hugo_Symbol %chin%
        genes_inventory
]


# ============================================================
# 5. Pre-specified high-confidence mutation definitions
# ============================================================

# ERBB2 activating mutations:
# literature-supported recurrent activating CRC variants.
erbb2_activating <- c(
    "p.S310F",
    "p.S310Y",
    "p.L755S",
    "p.V777L",
    "p.V842I",
    "p.L866M"
)

# PIK3CA exon-20 / kinase-domain hotspot:
# exon 20 has stronger anti-EGFR negative-predictive evidence
# than exon 9.
pik3ca_ex20 <- c(
    "p.M1043I",
    "p.H1047R",
    "p.H1047L",
    "p.H1047Y",
    "p.G1049R"
)

# PIK3CA exon-9 hotspots are descriptive only.
pik3ca_ex9 <- c(
    "p.E542K",
    "p.E545K",
    "p.E545G",
    "p.Q546K",
    "p.Q546R"
)

lof_classes <- c(
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Nonsense_Mutation",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation"
)

# MAP2K1:
# pre-specified activating/resistance hotspot codons.
# The observed D67E is retained because codon 67 is a recurrent
# MEK1 activation/resistance region; it is separately flagged so
# a sensitivity excluding D67E can be performed later.
map2k1_hotspot <- c(
    "p.Q56P",
    "p.K57E",
    "p.K57N",
    "p.K57T",
    "p.D67E",
    "p.D67N",
    "p.C121S",
    "p.P124S",
    "p.P124L",
    "p.E203K"
)

# Class III / kinase-impaired BRAF observed in this cohort.
# NOT considered bypass-positive.
braf_class3 <- c(
    "p.D594G",
    "p.D594N",
    "p.N581S"
)

mut_goi[
    ,
    ERBB2_ACTIVATING_R := as.integer(
        Hugo_Symbol == "ERBB2" &
            HGVSp_Short %chin%
                erbb2_activating
    )
]

mut_goi[
    ,
    PIK3CA_EX20_R := as.integer(
        Hugo_Symbol == "PIK3CA" &
            HGVSp_Short %chin%
                pik3ca_ex20
    )
]

mut_goi[
    ,
    PIK3CA_EX9_R := as.integer(
        Hugo_Symbol == "PIK3CA" &
            HGVSp_Short %chin%
                pik3ca_ex9
    )
]

mut_goi[
    ,
    PTEN_LOF_R := as.integer(
        Hugo_Symbol == "PTEN" &
            Variant_Classification %chin%
                lof_classes
    )
]

mut_goi[
    ,
    AKT1_E17K_R := as.integer(
        Hugo_Symbol == "AKT1" &
            HGVSp_Short == "p.E17K"
    )
]

mut_goi[
    ,
    MAP2K1_HOTSPOT_R := as.integer(
        Hugo_Symbol == "MAP2K1" &
            HGVSp_Short %chin%
                map2k1_hotspot
    )
]

mut_goi[
    ,
    MAP2K1_D67E_R := as.integer(
        Hugo_Symbol == "MAP2K1" &
            HGVSp_Short == "p.D67E"
    )
]

mut_goi[
    ,
    NF1_LOF_R := as.integer(
        Hugo_Symbol == "NF1" &
            Variant_Classification %chin%
                lof_classes
    )
]

mut_goi[
    ,
    BRAF_CLASS3_R := as.integer(
        Hugo_Symbol == "BRAF" &
            HGVSp_Short %chin%
                braf_class3
    )
]

mut_goi[
    ,
    EGFR_MUTATION_DESCRIPTIVE_R := as.integer(
        Hugo_Symbol == "EGFR"
    )
]

mut_goi[
    ,
    GNAS_MUTATION_DESCRIPTIVE_R := as.integer(
        Hugo_Symbol == "GNAS"
    )
]


# ============================================================
# 6. Patient mutation flags
# ============================================================

mutation_patient <- mut_goi[
    ,
    .(
        ERBB2_ACTIVATING_R =
            max(
                ERBB2_ACTIVATING_R,
                na.rm = TRUE
            ),

        PIK3CA_EX20_R =
            max(
                PIK3CA_EX20_R,
                na.rm = TRUE
            ),

        PIK3CA_EX9_R =
            max(
                PIK3CA_EX9_R,
                na.rm = TRUE
            ),

        PTEN_LOF_R =
            max(
                PTEN_LOF_R,
                na.rm = TRUE
            ),

        AKT1_E17K_R =
            max(
                AKT1_E17K_R,
                na.rm = TRUE
            ),

        MAP2K1_HOTSPOT_R =
            max(
                MAP2K1_HOTSPOT_R,
                na.rm = TRUE
            ),

        MAP2K1_D67E_R =
            max(
                MAP2K1_D67E_R,
                na.rm = TRUE
            ),

        NF1_LOF_R =
            max(
                NF1_LOF_R,
                na.rm = TRUE
            ),

        BRAF_CLASS3_R =
            max(
                BRAF_CLASS3_R,
                na.rm = TRUE
            ),

        EGFR_MUTATION_DESCRIPTIVE_R =
            max(
                EGFR_MUTATION_DESCRIPTIVE_R,
                na.rm = TRUE
            ),

        GNAS_MUTATION_DESCRIPTIVE_R =
            max(
                GNAS_MUTATION_DESCRIPTIVE_R,
                na.rm = TRUE
            )
    ),
    by = PATIENT_ID
]


# ============================================================
# 7. CNA: ERBB2 / MET / KRAS amplification
# ============================================================

cna_file <- find_raw(
    "data_cna.txt"
)

cna_header <- names(
    fread(
        cna_file,
        sep = "\t",
        header = TRUE,
        nrows = 0,
        check.names = FALSE,
        showProgress = FALSE
    )
)

usable_samples <- intersect(
    cohort$SAMPLE_ID,
    cna_header
)

if (length(usable_samples) != 191L) {

    stop(
        "Expected all 191 cohort samples in CNA matrix; found ",
        length(usable_samples)
    )
}

cna_small <- fread(
    cna_file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    select = c(
        "Hugo_Symbol",
        usable_samples
    ),
    showProgress = TRUE
)

cna_small <- cna_small[
    Hugo_Symbol %chin%
        c(
            "ERBB2",
            "MET",
            "KRAS"
        )
]

build_amp <- function(gene) {

    row_g <- cna_small[
        Hugo_Symbol == gene
    ]

    if (nrow(row_g) == 0L) {

        return(
            data.table(
                SAMPLE_ID = usable_samples,
                GENE = gene,
                AMP = 0L
            )
        )
    }

    vals <- suppressWarnings(
        as.numeric(
            unlist(
                row_g[
                    1,
                    ..usable_samples
                ],
                use.names = FALSE
            )
        )
    )

    data.table(
        SAMPLE_ID = usable_samples,
        GENE = gene,
        AMP = as.integer(
            vals == 2
        )
    )
}

amp_long <- rbindlist(
    lapply(
        c(
            "ERBB2",
            "MET",
            "KRAS"
        ),
        build_amp
    )
)

amp_wide <- dcast(
    amp_long,
    SAMPLE_ID ~ GENE,
    value.var = "AMP",
    fill = 0L
)

setnames(
    amp_wide,
    c(
        "ERBB2",
        "MET",
        "KRAS"
    ),
    c(
        "ERBB2_AMP_R",
        "MET_AMP_R",
        "KRAS_AMP_R"
    )
)

amp_wide[
    ,
    PATIENT_ID := sample_map$PATIENT_ID[
        match(
            SAMPLE_ID,
            sample_map$SAMPLE_ID
        )
    ]
]

amp_patient <- amp_wide[
    ,
    .(
        ERBB2_AMP_R =
            max(
                ERBB2_AMP_R,
                na.rm = TRUE
            ),

        MET_AMP_R =
            max(
                MET_AMP_R,
                na.rm = TRUE
            ),

        KRAS_AMP_R =
            max(
                KRAS_AMP_R,
                na.rm = TRUE
            )
    ),
    by = PATIENT_ID
]


# ============================================================
# 8. Canonical resistance fusions
# ============================================================

sv <- fread(
    find_raw(
        "data_sv.txt"
    ),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

sv <- sv[
    Sample_Id %chin%
        cohort$SAMPLE_ID
]

sv[
    ,
    PATIENT_ID := sample_map$PATIENT_ID[
        match(
            Sample_Id,
            sample_map$SAMPLE_ID
        )
    ]
]

canonical_fusion_genes <- c(
    "ALK",
    "ROS1",
    "NTRK1",
    "NTRK2",
    "NTRK3",
    "RET"
)

sv[
    ,
    CANONICAL_PRESSING_FUSION_R := as.integer(
        (
            Site1_Hugo_Symbol %chin%
                canonical_fusion_genes |
            Site2_Hugo_Symbol %chin%
                canonical_fusion_genes
        ) &
        grepl(
            "fusion",
            paste(
                safe_chr(Event_Info),
                safe_chr(Annotation)
            ),
            ignore.case = TRUE
        )
    )
]

fusion_patient <- sv[
    ,
    .(
        CANONICAL_PRESSING_FUSION_R =
            max(
                CANONICAL_PRESSING_FUSION_R,
                na.rm = TRUE
            )
    ),
    by = PATIENT_ID
]


# ============================================================
# 9. Treatment context / anti-EGFR agent
# ============================================================

treatment <- fread(
    find_raw(
        "data_timeline_treatment.txt"
    ),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

treatment <- treatment[
    PATIENT_ID %chin%
        cohort$PATIENT_ID
]

treatment[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(
            START_DATE
        )
    )
]

treatment[
    ,
    STOP_DATE_R := suppressWarnings(
        as.numeric(
            STOP_DATE
        )
    )
]

treatment[
    ,
    AGENT_UPPER_R := toupper(
        trimws(
            safe_chr(
                AGENT
            )
        )
    )
]

treatment[
    ,
    IS_FP_R :=
        AGENT_UPPER_R %chin%
        c(
            "FLUOROURACIL",
            "5-FLUOROURACIL",
            "5-FU",
            "CAPECITABINE"
        )
]

treatment[
    ,
    IS_OX_R :=
        AGENT_UPPER_R ==
        "OXALIPLATIN"
]

treatment[
    ,
    IS_IRI_R :=
        AGENT_UPPER_R ==
        "IRINOTECAN"
]

treatment[
    ,
    IS_EGFR_R :=
        AGENT_UPPER_R %chin%
        c(
            "CETUXIMAB",
            "PANITUMUMAB"
        )
]

build_tx_context <- function(pid, index_day) {

    p <- treatment[
        PATIENT_ID == pid &
            is.finite(
                START_DATE_R
            )
    ]

    if (nrow(p) == 0L) {

        return(
            data.table(
                PATIENT_ID = pid,
                ANTI_EGFR_AGENT_R = NA_character_,
                PRIOR_HEAVY_R = NA_integer_,
                PRIOR_CYTOTOXIC_PATTERN_R = NA_character_,
                ACTIVE_BACKBONE_AT_INDEX_R = NA_character_
            )
        )
    }

    index_egfr <- p[
        IS_EGFR_R == TRUE &
            START_DATE_R ==
                index_day
    ]

    anti_agent <- collapse_unique(
        index_egfr$AGENT
    )

    prior <- p[
        START_DATE_R <
            index_day
    ]

    prior_fp <- any(
        prior$IS_FP_R
    )

    prior_ox <- any(
        prior$IS_OX_R
    )

    prior_iri <- any(
        prior$IS_IRI_R
    )

    prior_pattern <- fcase(
        prior_ox & prior_iri,
        "prior oxaliplatin + irinotecan",

        xor(
            prior_ox,
            prior_iri
        ),
        "prior one major backbone",

        prior_fp &
            !prior_ox &
            !prior_iri,
        "prior fluoropyrimidine only",

        default =
            "no prior major cytotoxic"
    )

    # More parsimonious binary covariate for later primary model
    prior_heavy <- as.integer(
        prior_ox &
            prior_iri
    )

    # Active regimen around index:
    # overlapping medication or starts within 7 days after index.
    active <- p[
        (
            START_DATE_R <=
                index_day &
            (
                is.na(
                    STOP_DATE_R
                ) |
                STOP_DATE_R >=
                    index_day
            )
        ) |
        (
            START_DATE_R >
                index_day &
            START_DATE_R <=
                index_day + 7
        )
    ]

    has_fp <- any(
        active$IS_FP_R
    )

    has_ox <- any(
        active$IS_OX_R
    )

    has_iri <- any(
        active$IS_IRI_R
    )

    active_backbone <- fcase(
        has_fp & has_ox & has_iri,
        "triplet/FOLFOXIRI-like",

        has_fp & has_ox & !has_iri,
        "FOLFOX/CAPOX-like",

        has_fp & has_iri & !has_ox,
        "FOLFIRI/CAPIRI-like",

        has_fp & !has_ox & !has_iri,
        "fluoropyrimidine-only",

        !has_fp & !has_ox & has_iri,
        "irinotecan-only",

        !has_fp & has_ox & !has_iri,
        "oxaliplatin-only",

        default =
            "no major cytotoxic active"
    )

    data.table(
        PATIENT_ID = pid,
        ANTI_EGFR_AGENT_R =
            anti_agent,
        PRIOR_HEAVY_R =
            prior_heavy,
        PRIOR_CYTOTOXIC_PATTERN_R =
            prior_pattern,
        ACTIVE_BACKBONE_AT_INDEX_R =
            active_backbone
    )
}

tx_context <- rbindlist(
    Map(
        build_tx_context,
        cohort$PATIENT_ID,
        cohort$FIRST_BIOLOGIC_DAY_R
    ),
    fill = TRUE
)


# ============================================================
# 10. Explicit CRC Stage IV date
# ============================================================

diagnosis <- fread(
    find_raw(
        "data_timeline_diagnosis.txt"
    ),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

diagnosis <- diagnosis[
    PATIENT_ID %chin%
        cohort$PATIENT_ID
]

diagnosis[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(
            START_DATE
        )
    )
]

diagnosis[
    ,
    SITE_CODE_R :=
        extract_site_code(
            DX_DESCRIPTION
        )
]

diagnosis[
    ,
    IS_CRC_DX_R := grepl(
        "^C18[0-9]$|^C19[0-9]$|^C20[0-9]$|^C218$",
        SITE_CODE_R,
        perl = TRUE
    )
]

diagnosis[
    ,
    IS_STAGE4_CRC_R :=
        IS_CRC_DX_R &
        (
            STAGE_CDM_DERIVED ==
                "Stage 4" |
            grepl(
                "Distant",
                SUMMARY,
                ignore.case = TRUE
            )
        )
]

stage4 <- diagnosis[
    IS_STAGE4_CRC_R == TRUE,
    .(
        EXPLICIT_STAGE4_CRC_DAY_R =
            min_or_na(
                START_DATE_R
            )
    ),
    by = PATIENT_ID
]


# ============================================================
# 11. ECOG performance status
# ============================================================

performance <- fread(
    find_raw(
        "data_timeline_performance_status.txt"
    ),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

required_ps <- c(
    "PATIENT_ID",
    "START_DATE",
    "ECOG"
)

if (!all(required_ps %in% names(performance))) {

    stop(
        "Performance-status file header mismatch."
    )
}

performance <- performance[
    PATIENT_ID %chin%
        cohort$PATIENT_ID
]

performance[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(
            START_DATE
        )
    )
]

performance[
    ,
    ECOG_R := suppressWarnings(
        as.numeric(
            ECOG
        )
    )
]

build_ecog <- function(pid, index_day) {

    p <- performance[
        PATIENT_ID == pid &
            is.finite(
                START_DATE_R
            ) &
            START_DATE_R <=
                index_day &
            is.finite(
                ECOG_R
            )
    ]

    if (nrow(p) == 0L) {

        return(
            data.table(
                PATIENT_ID = pid,
                ECOG_LATEST_PRE_R = NA_real_,
                ECOG_LAG_DAYS_R = NA_real_,
                ECOG_WITHIN_90D_R = NA_real_,
                ECOG_WITHIN_180D_R = NA_real_
            )
        )
    }

    setorder(
        p,
        -START_DATE_R
    )

    r <- p[1]

    lag <- index_day -
        r$START_DATE_R

    data.table(
        PATIENT_ID = pid,
        ECOG_LATEST_PRE_R =
            r$ECOG_R,
        ECOG_LAG_DAYS_R =
            lag,
        ECOG_WITHIN_90D_R =
            fifelse(
                lag <= 90,
                r$ECOG_R,
                NA_real_
            ),
        ECOG_WITHIN_180D_R =
            fifelse(
                lag <= 180,
                r$ECOG_R,
                NA_real_
            )
    )
}

ecog_context <- rbindlist(
    Map(
        build_ecog,
        cohort$PATIENT_ID,
        cohort$FIRST_BIOLOGIC_DAY_R
    ),
    fill = TRUE
)


# ============================================================
# 12. Metastatic-site burden within 90 days before index
# ============================================================

tumor_sites <- fread(
    find_raw(
        "data_timeline_tumor_sites.txt"
    ),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

tumor_sites <- tumor_sites[
    PATIENT_ID %chin%
        cohort$PATIENT_ID
]

tumor_sites[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(
            START_DATE
        )
    )
]

build_site_burden <- function(pid, index_day) {

    p <- tumor_sites[
        PATIENT_ID == pid &
            is.finite(
                START_DATE_R
            ) &
            START_DATE_R >=
                index_day - 90 &
            START_DATE_R <=
                index_day
    ]

    if (nrow(p) == 0L) {

        return(
            data.table(
                PATIENT_ID = pid,
                MET_SITE_COUNT_90D_R = NA_integer_,
                MET_SITE_TYPES_90D_R = NA_character_
            )
        )
    }

    sites <- trimws(
        safe_chr(
            p$TUMOR_SITE
        )
    )

    sites <- sort(
        unique(
            sites[
                nzchar(sites)
            ]
        )
    )

    data.table(
        PATIENT_ID = pid,
        MET_SITE_COUNT_90D_R =
            length(sites),
        MET_SITE_TYPES_90D_R =
            if (length(sites) == 0L) {
                NA_character_
            } else {
                paste(
                    sites,
                    collapse = "|"
                )
            }
    )
}

site_context <- rbindlist(
    Map(
        build_site_burden,
        cohort$PATIENT_ID,
        cohort$FIRST_BIOLOGIC_DAY_R
    ),
    fill = TRUE
)


# ============================================================
# 13. Assemble patient-level outcome-blind audit
# ============================================================

audit <- merge(
    cohort,
    sample_primary,
    by = c(
        "PATIENT_ID",
        "SAMPLE_ID"
    ),
    all.x = TRUE,
    sort = FALSE
)

for (obj in list(
    mutation_patient,
    amp_patient,
    fusion_patient,
    tx_context,
    stage4,
    ecog_context,
    site_context
)) {

    audit <- merge(
        audit,
        obj,
        by = "PATIENT_ID",
        all.x = TRUE,
        sort = FALSE
    )
}

binary_flags <- c(
    "ERBB2_ACTIVATING_R",
    "PIK3CA_EX20_R",
    "PIK3CA_EX9_R",
    "PTEN_LOF_R",
    "AKT1_E17K_R",
    "MAP2K1_HOTSPOT_R",
    "MAP2K1_D67E_R",
    "NF1_LOF_R",
    "BRAF_CLASS3_R",
    "EGFR_MUTATION_DESCRIPTIVE_R",
    "GNAS_MUTATION_DESCRIPTIVE_R",
    "ERBB2_AMP_R",
    "MET_AMP_R",
    "KRAS_AMP_R",
    "CANONICAL_PRESSING_FUSION_R"
)

for (v in binary_flags) {

    if (!v %in% names(audit)) {
        audit[
            ,
            (v) := 0L
        ]
    }

    audit[
        is.na(
            get(v)
        ),
        (v) := 0L
    ]
}


# ============================================================
# 14. FROZEN exposure tiers
# ============================================================

# Tier 1: strict PRESSING-like core
audit[
    ,
    BYPASS_PRESSING_STRICT_R := as.integer(
        ERBB2_AMP_R == 1L |
        ERBB2_ACTIVATING_R == 1L |
        MET_AMP_R == 1L |
        PIK3CA_EX20_R == 1L |
        PTEN_LOF_R == 1L |
        AKT1_E17K_R == 1L |
        CANONICAL_PRESSING_FUSION_R == 1L
    )
]

# Tier 2: PRIMARY high-confidence panel
audit[
    ,
    BYPASS_HIGH_CONFIDENCE_R := as.integer(
        BYPASS_PRESSING_STRICT_R == 1L |
        KRAS_AMP_R == 1L |
        MAP2K1_HOTSPOT_R == 1L |
        NF1_LOF_R == 1L
    )
]

# Sensitivity excluding the single MAP2K1 D67E if present
audit[
    ,
    BYPASS_HIGH_CONFIDENCE_NO_D67E_R := as.integer(
        BYPASS_PRESSING_STRICT_R == 1L |
        KRAS_AMP_R == 1L |
        (
            MAP2K1_HOTSPOT_R == 1L &
            MAP2K1_D67E_R == 0L
        ) |
        NF1_LOF_R == 1L
    )
]

# Tier 3: old broad composite retained as sensitivity only
audit[
    ,
    BYPASS_ORIGINAL_BROAD_R :=
        BYPASS_ORIGINAL_R
]


# ============================================================
# 15. Time / specimen / assay covariates
# ============================================================

audit[
    ,
    STAGE4_TO_ANTI_EGFR_DAYS_R :=
        FIRST_BIOLOGIC_DAY_R -
        EXPLICIT_STAGE4_CRC_DAY_R
]

audit[
    ,
    LOG_STAGE4_TO_ANTI_EGFR_R :=
        log1p(
            pmax(
                STAGE4_TO_ANTI_EGFR_DAYS_R,
                0
            )
        )
]

audit[
    ,
    SPECIMEN_TO_ANTI_EGFR_DAYS_R :=
        FIRST_BIOLOGIC_DAY_R -
        SPECIMEN_ACQUISITION_FIRST_DAY_R
]

audit[
    ,
    SEQUENCING_TO_ANTI_EGFR_DAYS_R :=
        FIRST_BIOLOGIC_DAY_R -
        SEQUENCING_FIRST_DAY_R
]

audit[
    ,
    SPECIMEN_WITHIN_365D_R := as.integer(
        is.finite(
            SPECIMEN_TO_ANTI_EGFR_DAYS_R
        ) &
        SPECIMEN_TO_ANTI_EGFR_DAYS_R <=
            365
    )
]

audit[
    ,
    SPECIMEN_WITHIN_730D_R := as.integer(
        is.finite(
            SPECIMEN_TO_ANTI_EGFR_DAYS_R
        ) &
        SPECIMEN_TO_ANTI_EGFR_DAYS_R <=
            730
    )
]

audit[
    ,
    PANEL_OLD_R := as.integer(
        GENE_PANEL %chin%
            c(
                "IMPACT341",
                "IMPACT410"
            )
    )
]

primary_site_lower <- tolower(
    trimws(
        safe_chr(
            audit$PRIMARY_SITE
        )
    )
)

audit[
    ,
    LEFT_SUBSITE_R := fcase(
        grepl(
            "rect",
            primary_site_lower
        ),
        "rectum/rectosigmoid",

        grepl(
            "sigmoid|descending|left",
            primary_site_lower
        ),
        "left colon",

        default =
            "other/unknown"
    )
]


# ============================================================
# 16. Hard outcome-blind anchors from independent raw precheck
# ============================================================

n_strict <- audit[
    BYPASS_PRESSING_STRICT_R == 1L,
    .N
]

n_primary <- audit[
    BYPASS_HIGH_CONFIDENCE_R == 1L,
    .N
]

n_primary_no_d67e <- audit[
    BYPASS_HIGH_CONFIDENCE_NO_D67E_R == 1L,
    .N
]

n_broad <- audit[
    BYPASS_ORIGINAL_BROAD_R == 1L,
    .N
]

n_braf3 <- audit[
    BRAF_CLASS3_R == 1L,
    .N
]

n_pik3ca_ex9 <- audit[
    PIK3CA_EX9_R == 1L,
    .N
]

expected_ok <- (
    nrow(audit) == 191L &&
        n_strict == 15L &&
        n_primary == 23L &&
        n_primary_no_d67e == 22L &&
        n_broad == 35L &&
        n_braf3 == 7L &&
        n_pik3ca_ex9 == 6L
)

if (!expected_ok) {

    stop(
        paste0(
            "B1-06 hard audit failed.\n",
            "Expected cohort/strict/primary/noD67E/broad/BRAF3/PIK3CAex9 = ",
            "191/15/23/22/35/7/6.\n",
            "Observed = ",
            nrow(audit), "/",
            n_strict, "/",
            n_primary, "/",
            n_primary_no_d67e, "/",
            n_broad, "/",
            n_braf3, "/",
            n_pik3ca_ex9
        )
    )
}

cat(
    "Hard anchors PASS:",
    "191 / 15 / 23 / 22 / 35 / 7 / 6\n"
)


# ============================================================
# 17. Component counts
# ============================================================

component_counts <- data.table(
    COMPONENT = c(
        "ERBB2 amplification",
        "ERBB2 activating mutation",
        "MET amplification",
        "PIK3CA exon20 hotspot",
        "PTEN loss-of-function",
        "AKT1 E17K",
        "Canonical ALK/ROS1/NTRK/RET fusion",
        "KRAS amplification",
        "MAP2K1 hotspot",
        "NF1 loss-of-function",
        "PIK3CA exon9 hotspot (not primary)",
        "BRAF class III (not bypass)",
        "PRESSING-like strict",
        "PRIMARY high-confidence bypass",
        "PRIMARY high-confidence excluding D67E",
        "Original broad sensitivity"
    ),

    N = c(
        audit[
            ERBB2_AMP_R == 1L,
            .N
        ],
        audit[
            ERBB2_ACTIVATING_R == 1L,
            .N
        ],
        audit[
            MET_AMP_R == 1L,
            .N
        ],
        audit[
            PIK3CA_EX20_R == 1L,
            .N
        ],
        audit[
            PTEN_LOF_R == 1L,
            .N
        ],
        audit[
            AKT1_E17K_R == 1L,
            .N
        ],
        audit[
            CANONICAL_PRESSING_FUSION_R == 1L,
            .N
        ],
        audit[
            KRAS_AMP_R == 1L,
            .N
        ],
        audit[
            MAP2K1_HOTSPOT_R == 1L,
            .N
        ],
        audit[
            NF1_LOF_R == 1L,
            .N
        ],
        n_pik3ca_ex9,
        n_braf3,
        n_strict,
        n_primary,
        n_primary_no_d67e,
        n_broad
    )
)

component_counts[
    ,
    PERCENT :=
        round(
            100 *
                N /
                nrow(audit),
            1
        )
]


# ============================================================
# 18. Variant inventory
# ============================================================

variant_inventory <- mut_goi[
    ,
    .(
        N_PATIENTS =
            uniqueN(
                PATIENT_ID
            ),

        ERBB2_ACTIVATING =
            max(
                ERBB2_ACTIVATING_R,
                na.rm = TRUE
            ),

        PIK3CA_EX20 =
            max(
                PIK3CA_EX20_R,
                na.rm = TRUE
            ),

        PIK3CA_EX9 =
            max(
                PIK3CA_EX9_R,
                na.rm = TRUE
            ),

        PTEN_LOF =
            max(
                PTEN_LOF_R,
                na.rm = TRUE
            ),

        AKT1_E17K =
            max(
                AKT1_E17K_R,
                na.rm = TRUE
            ),

        MAP2K1_HOTSPOT =
            max(
                MAP2K1_HOTSPOT_R,
                na.rm = TRUE
            ),

        NF1_LOF =
            max(
                NF1_LOF_R,
                na.rm = TRUE
            ),

        BRAF_CLASS3 =
            max(
                BRAF_CLASS3_R,
                na.rm = TRUE
            )
    ),
    by = .(
        Hugo_Symbol,
        HGVSp_Short,
        Variant_Classification
    )
][order(
    Hugo_Symbol,
    -N_PATIENTS,
    HGVSp_Short
)]


# ============================================================
# 19. Outcome-blind covariate availability
# ============================================================

covariate_availability <- data.table(
    VARIABLE = c(
        "Sex",
        "Left subsite",
        "Anti-EGFR agent",
        "Stage IV -> anti-EGFR interval",
        "Prior heavy exposure: both oxaliplatin + irinotecan",
        "Detailed prior cytotoxic pattern",
        "Active backbone at anti-EGFR index",
        "Sample type",
        "Gene panel",
        "Specimen -> anti-EGFR interval",
        "Sequencing -> anti-EGFR interval",
        "Latest ECOG on/before index",
        "ECOG within 90 days",
        "ECOG within 180 days",
        "Metastatic-site burden within 90 days"
    ),

    NONMISSING_N = c(
        audit[
            nzchar(
                safe_chr(
                    GENDER
                )
            ),
            .N
        ],
        audit[
            nzchar(
                safe_chr(
                    LEFT_SUBSITE_R
                )
            ),
            .N
        ],
        audit[
            nzchar(
                safe_chr(
                    ANTI_EGFR_AGENT_R
                )
            ),
            .N
        ],
        audit[
            is.finite(
                STAGE4_TO_ANTI_EGFR_DAYS_R
            ),
            .N
        ],
        audit[
            !is.na(
                PRIOR_HEAVY_R
            ),
            .N
        ],
        audit[
            nzchar(
                safe_chr(
                    PRIOR_CYTOTOXIC_PATTERN_R
                )
            ),
            .N
        ],
        audit[
            nzchar(
                safe_chr(
                    ACTIVE_BACKBONE_AT_INDEX_R
                )
            ),
            .N
        ],
        audit[
            nzchar(
                safe_chr(
                    SAMPLE_TYPE
                )
            ),
            .N
        ],
        audit[
            nzchar(
                safe_chr(
                    GENE_PANEL
                )
            ),
            .N
        ],
        audit[
            is.finite(
                SPECIMEN_TO_ANTI_EGFR_DAYS_R
            ),
            .N
        ],
        audit[
            is.finite(
                SEQUENCING_TO_ANTI_EGFR_DAYS_R
            ),
            .N
        ],
        audit[
            is.finite(
                ECOG_LATEST_PRE_R
            ),
            .N
        ],
        audit[
            is.finite(
                ECOG_WITHIN_90D_R
            ),
            .N
        ],
        audit[
            is.finite(
                ECOG_WITHIN_180D_R
            ),
            .N
        ],
        audit[
            !is.na(
                MET_SITE_COUNT_90D_R
            ),
            .N
        ]
    )
)

covariate_availability[
    ,
    NONMISSING_PERCENT :=
        round(
            100 *
                NONMISSING_N /
                nrow(audit),
            1
        )
]


# ============================================================
# 20. Outcome-blind distributions by PRIMARY bypass
# ============================================================

primary_group_n <- audit[
    ,
    .N,
    by = BYPASS_HIGH_CONFIDENCE_R
][order(BYPASS_HIGH_CONFIDENCE_R)]

sample_type_tab <- audit[
    ,
    .N,
    by = .(
        BYPASS_HIGH_CONFIDENCE_R,
        SAMPLE_TYPE
    )
][order(
    BYPASS_HIGH_CONFIDENCE_R,
    -N
)]

panel_tab <- audit[
    ,
    .N,
    by = .(
        BYPASS_HIGH_CONFIDENCE_R,
        GENE_PANEL
    )
][order(
    BYPASS_HIGH_CONFIDENCE_R,
    GENE_PANEL
)]

agent_tab <- audit[
    ,
    .N,
    by = .(
        BYPASS_HIGH_CONFIDENCE_R,
        ANTI_EGFR_AGENT_R
    )
][order(
    BYPASS_HIGH_CONFIDENCE_R,
    ANTI_EGFR_AGENT_R
)]

prior_tab <- audit[
    ,
    .N,
    by = .(
        BYPASS_HIGH_CONFIDENCE_R,
        PRIOR_CYTOTOXIC_PATTERN_R
    )
][order(
    BYPASS_HIGH_CONFIDENCE_R,
    -N
)]

subsite_tab <- audit[
    ,
    .N,
    by = .(
        BYPASS_HIGH_CONFIDENCE_R,
        LEFT_SUBSITE_R
    )
][order(
    BYPASS_HIGH_CONFIDENCE_R,
    LEFT_SUBSITE_R
)]

timing_summary <- audit[
    ,
    .(
        STAGE4_TO_INDEX =
            median_iqr(
                STAGE4_TO_ANTI_EGFR_DAYS_R
            ),

        SPECIMEN_TO_INDEX =
            median_iqr(
                SPECIMEN_TO_ANTI_EGFR_DAYS_R
            ),

        SEQUENCING_TO_INDEX =
            median_iqr(
                SEQUENCING_TO_ANTI_EGFR_DAYS_R
            ),

        SPECIMEN_WITHIN_365D =
            sum(
                SPECIMEN_WITHIN_365D_R,
                na.rm = TRUE
            ),

        SPECIMEN_WITHIN_730D =
            sum(
                SPECIMEN_WITHIN_730D_R,
                na.rm = TRUE
            ),

        N = .N
    ),
    by = BYPASS_HIGH_CONFIDENCE_R
][order(BYPASS_HIGH_CONFIDENCE_R)]


# ============================================================
# 21. Frozen covariate strategy
# ============================================================

# Keep the PRIMARY model parsimonious because only 23 patients
# are high-confidence bypass-positive.
#
# PRIMARY adjustment:
#   sex
#   left subsite
#   anti-EGFR agent
#   log(Stage IV -> anti-EGFR interval)
#   prior heavy cytotoxic exposure (both oxaliplatin + irinotecan)
#
# Sensitivity additions / restrictions:
#   sample type
#   specimen-to-index interval
#   old vs newer IMPACT panel
#   active backbone at index
#   ECOG complete-case within 180 days
#   metastatic-site burden within 90 days
#   specimen acquired <=730 days before anti-EGFR
#   exclude BRAF class III
#   exclude MAP2K1 D67E
#
primary_covariates <- c(
    "GENDER",
    "LEFT_SUBSITE_R",
    "ANTI_EGFR_AGENT_R",
    "LOG_STAGE4_TO_ANTI_EGFR_R",
    "PRIOR_HEAVY_R"
)

sensitivity_covariates <- c(
    "SAMPLE_TYPE",
    "SPECIMEN_TO_ANTI_EGFR_DAYS_R",
    "PANEL_OLD_R",
    "ACTIVE_BACKBONE_AT_INDEX_R",
    "ECOG_WITHIN_180D_R",
    "MET_SITE_COUNT_90D_R"
)


# ============================================================
# 22. Write audit outputs
# ============================================================

variant_file <- file.path(
    audit_dir,
    "B1_06_variant_inventory.csv"
)

patient_file <- file.path(
    audit_dir,
    "B1_06_patient_genomic_covariate_audit.csv"
)

audit_file <- file.path(
    audit_dir,
    "B1_06_bypass_covariate_freeze_audit.txt"
)

fwrite(
    variant_inventory,
    variant_file,
    na = ""
)

fwrite(
    audit,
    patient_file,
    na = ""
)

lines <- c(
    "B1-06 OUTCOME-BLIND BYPASS / COVARIATE FREEZE",
    "================================================",
    "",
    "IMPORTANT",
    "---------",
    "No progression, rwPFS, death, OS, HR, p-value, KM or Cox result was read or calculated.",
    "",
    "1. Primary cohort",
    "-----------------",
    paste0(
        "Sequencing-pre-first-anti-EGFR cohort: ",
        nrow(audit)
    ),
    "",
    "2. Frozen genomic exposure tiers",
    "--------------------------------",
    "",
    "Tier 1: PRESSING-like strict core",
    "  ERBB2 amplification",
    "  OR high-confidence activating ERBB2 mutation",
    "  OR MET amplification",
    "  OR PIK3CA exon20/kinase-domain hotspot",
    "  OR clear PTEN loss-of-function",
    "  OR AKT1 E17K",
    "  OR canonical ALK/ROS1/NTRK1-3/RET fusion",
    paste0(
        "  N positive = ",
        n_strict,
        "/191"
    ),
    "",
    "Tier 2: PRIMARY high-confidence EGFR-bypass panel",
    "  Tier 1 PLUS KRAS amplification, MAP2K1 hotspot, or clear NF1 loss-of-function",
    paste0(
        "  N positive = ",
        n_primary,
        "/191"
    ),
    "",
    "Pre-specified MAP2K1 sensitivity:",
    paste0(
        "  Excluding MAP2K1 D67E -> N positive = ",
        n_primary_no_d67e,
        "/191"
    ),
    "",
    "Tier 3: Original broad composite",
    "  Existing BYPASS_ORIGINAL_R retained for sensitivity only.",
    paste0(
        "  N positive = ",
        n_broad,
        "/191"
    ),
    "",
    "Non-V600 BRAF:",
    paste0(
        "  Class III BRAF (D594G/D594N/N581S) = ",
        n_braf3,
        "/191"
    ),
    "  These patients are NOT called bypass-positive solely because of BRAF.",
    "  A sensitivity analysis will exclude them.",
    "",
    "PIK3CA exon 9:",
    paste0(
        "  Exon-9 hotspot = ",
        n_pik3ca_ex9,
        "/191"
    ),
    "  Not included in the primary bypass definition.",
    "",
    "Component counts:",
    capture.output(
        print(
            component_counts
        )
    ),
    "",
    "3. Outcome-blind covariate availability",
    "---------------------------------------",
    capture.output(
        print(
            covariate_availability
        )
    ),
    "",
    "4. PRIMARY high-confidence group sizes",
    "--------------------------------------",
    capture.output(
        print(
            primary_group_n
        )
    ),
    "",
    "Sample type:",
    capture.output(
        print(
            sample_type_tab
        )
    ),
    "",
    "Gene panel:",
    capture.output(
        print(
            panel_tab
        )
    ),
    "",
    "Anti-EGFR agent:",
    capture.output(
        print(
            agent_tab
        )
    ),
    "",
    "Prior cytotoxic exposure:",
    capture.output(
        print(
            prior_tab
        )
    ),
    "",
    "Left-sided subsite:",
    capture.output(
        print(
            subsite_tab
        )
    ),
    "",
    "Timing:",
    capture.output(
        print(
            timing_summary
        )
    ),
    "",
    "5. Frozen covariate strategy",
    "----------------------------",
    "PRIMARY adjusted model:",
    paste(
        primary_covariates,
        collapse = " + "
    ),
    "",
    "Sensitivity covariates / restrictions:",
    paste(
        sensitivity_covariates,
        collapse = " + "
    ),
    "",
    "Pre-specified restrictions:",
    "  specimen acquisition <=730 days before anti-EGFR",
    "  exclude BRAF class III",
    "  exclude MAP2K1 D67E from bypass classification",
    "",
    "ECOG:",
    "  Not used in the full-cohort primary model because recent ECOG is incompletely observed.",
    "  Use complete-case sensitivity with ECOG within 180 days.",
    "",
    "Age:",
    "  Do NOT use CURRENT_AGE_DEID as age at anti-EGFR index.",
    "  MSK-CHORD does not provide a valid age-at-index variable here.",
    "",
    "6. Manuscript terminology freeze",
    "--------------------------------",
    "Do NOT call this cohort 'BRAF wild-type'.",
    "Use: left-sided, RAS wild-type, BRAF V600E-negative, MSS mCRC.",
    "",
    "Do NOT call the tumor specimens 'immediate pretreatment biopsies'.",
    "Use: pre-treatment sequenced tumor / archival baseline genomic profile.",
    "",
    "Do NOT claim predictive treatment-effect modification from this single anti-EGFR cohort.",
    "Primary claim: baseline high-confidence EGFR-bypass alterations are associated with outcomes after first anti-EGFR exposure.",
    "",
    "7. Next-step gate",
    "-----------------",
    "If these outcome-blind definitions are accepted, the next script will LOCK the analysis cohort and then run the pre-specified rwPFS/OS models."
)

writeLines(
    lines,
    audit_file
)

unlink(
    tmp_dir,
    recursive = TRUE,
    force = TRUE
)

cat("\n============================================================\n")
cat("B1-06 FREEZE AUDIT COMPLETE\n")
cat("============================================================\n")
cat(
    "Cohort / strict / primary / primary-noD67E / broad / BRAF3 = ",
    nrow(audit), " / ",
    n_strict, " / ",
    n_primary, " / ",
    n_primary_no_d67e, " / ",
    n_broad, " / ",
    n_braf3,
    "\n",
    sep = ""
)
cat("\nMain audit:\n", audit_file, "\n")
cat("\nVariant inventory:\n", variant_file, "\n")
cat("\nPatient audit:\n", patient_file, "\n")
