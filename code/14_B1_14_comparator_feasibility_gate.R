# ============================================================
# A7 / B-line
# 14_B1_14_bevacizumab_comparator_feasibility_gate.R
#
# OUTCOME-BLIND COMPARATOR FEASIBILITY GATE ONLY
#
# PURPOSE
#   Determine whether the sequencing-pre-biologic bevacizumab cohort
#   has adequate high-confidence EGFR-bypass positivity and clinical
#   overlap with the locked anti-EGFR cohort to justify a later,
#   explicitly exploratory treatment x bypass analysis.
#
# ABSOLUTE RED LINES
#   - NO rwPFS / OS / death / progression is read or modeled.
#   - NO treatment-effect estimate or interaction is calculated.
#   - The frozen high-confidence genomic definition is NOT changed.
#   - The anti-EGFR primary analysis is NOT changed.
#   - A GO result means only "an exploratory comparator analysis may
#     be estimable"; it does NOT establish causal comparability.
#
# COHORTS
#   anti-EGFR reference:
#     OLD_ELIGIBLE_803_R == 1
#     FIRST_BIOLOGIC_R == "anti-EGFR"
#     sequencing completed on/before first biologic
#     expected N = 191
#
#   bevacizumab comparator feasibility cohort:
#     OLD_ELIGIBLE_803_R == 1
#     FIRST_BIOLOGIC_R == "bevacizumab"
#     sequencing completed on/before first biologic
#     expected N = 137
#
# GENOMIC ALGORITHM
#   Reconstructed from raw mutation/CNA/SV files using the exact
#   B1-06 high-confidence rules.
#
# HARD ANTI-EGFR VALIDATION BEFORE BEV INTERPRETATION
#   N / strict / high-confidence / no-D67E / broad / BRAF III /
#   PIK3CA exon9 =
#     191 / 15 / 23 / 22 / 35 / 7 / 6
#
#   component counts =
#     ERBB2 amp 8
#     ERBB2 activating 3
#     MET amp 1
#     PIK3CA exon20 3
#     PTEN LoF 2
#     AKT1 E17K 0
#     canonical fusion 0
#     KRAS amp 3
#     MAP2K1 hotspot 2
#     NF1 LoF 3
#
# FEASIBILITY / OVERLAP GATE DECLARED BEFORE RESULTS
#
#   HARD COUNT REQUIREMENTS
#     - high-confidence bypass+ >=10 in EACH treatment arm
#     - bypass- >=50 in EACH treatment arm
#
#   PROPENSITY-OVERLAP REQUIREMENTS
#     Outcome-blind logistic PS:
#       anti-EGFR vs bevacizumab ~
#       sex + frozen left subsite +
#       log(Stage IV -> first biologic) +
#       prior major cytotoxic classes (0/1/2) +
#       active irinotecan + active oxaliplatin +
#       old panel + metastatic specimen +
#       log(specimen -> biologic interval) +
#       recent metastatic-site count
#
#     GO requires:
#       - PS complete-case coverage >=75% in EACH arm
#       - >=80% of EACH arm inside PS common support
#       - >=70% of bypass+ patients in EACH arm inside common support
#       - overlap-weight ESS >=50 in EACH arm
#       - overlap-weight ESS among bypass+ >=8 in EACH arm
#
#   DECISION
#     GO:
#       hard count requirements + all overlap requirements pass.
#
#     CONDITIONAL GO:
#       hard counts pass and PS model is estimable, but >=1 overlap
#       threshold misses. DO NOT automatically run interaction;
#       review overlap diagnostics first.
#
#     NO-GO:
#       hard count requirement fails, PS model is not estimable, or
#       there is no PS common-support interval.
#
# OUTPUTS
#   04_results/B1_14_bevacizumab_comparator_feasibility.rds
#   06_logs_and_audit/B1_14_bevacizumab_comparator_feasibility.txt
#   07_tables/TableS7_bevacizumab_comparator_feasibility.docx
#
# If NO-GO:
#   comparator feasibility must still be reported in the manuscript
#   Discussion rather than hidden only in the Supplement.
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

raw_dir <- file.path(
    project_root,
    "01_raw_data"
)

results_dir <- file.path(
    project_root,
    "04_results"
)

audit_dir <- file.path(
    project_root,
    "06_logs_and_audit"
)

table_dir <- file.path(
    project_root,
    "07_tables"
)

main_file <- file.path(
    raw_dir,
    "A7_MSK_CHORD_CRC_clean_R.csv"
)

tar_file <- file.path(
    raw_dir,
    "msk_chord_2024.tar"
)

result_file <- file.path(
    results_dir,
    "B1_14_bevacizumab_comparator_feasibility.rds"
)

audit_file <- file.path(
    audit_dir,
    "B1_14_bevacizumab_comparator_feasibility.txt"
)

table_file <- file.path(
    table_dir,
    "TableS7_bevacizumab_comparator_feasibility.docx"
)

for (d0 in c(
    results_dir,
    audit_dir,
    table_dir
)) {
    if (!dir.exists(d0)) {
        dir.create(
            d0,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }
}

for (f0 in c(
    main_file,
    tar_file
)) {
    if (!file.exists(f0)) {
        stop(
            "Missing required input:\n",
            f0
        )
    }
}

packages <- c(
    "data.table",
    "officer",
    "flextable"
)

for (pkg in packages) {
    if (!requireNamespace(
        pkg,
        quietly = TRUE
    )) {
        install.packages(pkg)
    }
}

library(data.table)
library(officer)
library(flextable)

cat("\n============================================================\n")
cat("B1-14 BEVACIZUMAB COMPARATOR FEASIBILITY GATE\n")
cat("OUTCOME-BLIND ONLY\n")
cat("============================================================\n\n")


# ============================================================
# 1. Helpers
# ============================================================

safe_chr <- function(x) {
    z <- as.character(x)
    z[is.na(z)] <- ""
    z
}

min_or_na <- function(x) {

    z <- suppressWarnings(
        as.numeric(x)
    )

    z <- z[
        is.finite(z)
    ]

    if (
        length(z) == 0L
    ) {
        return(NA_real_)
    }

    min(z)
}

collapse_unique <- function(x) {

    z <- trimws(
        safe_chr(x)
    )

    z <- unique(
        z[
            nzchar(z)
        ]
    )

    if (
        length(z) == 0L
    ) {
        return(NA_character_)
    }

    paste(
        sort(z),
        collapse = "|"
    )
}

median_iqr <- function(x) {

    z <- suppressWarnings(
        as.numeric(x)
    )

    z <- z[
        is.finite(z)
    ]

    if (
        length(z) == 0L
    ) {
        return("NA")
    }

    q <- quantile(
        z,
        probs = c(
            0.25,
            0.50,
            0.75
        ),
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

pct_chr <- function(n, den) {

    if (
        is.na(den) ||
        den <= 0L
    ) {
        return("NA")
    }

    sprintf(
        "%d/%d (%.1f%%)",
        n,
        den,
        100 *
            n /
            den
    )
}

ess <- function(w) {

    z <- suppressWarnings(
        as.numeric(w)
    )

    z <- z[
        is.finite(z) &
            z >= 0
    ]

    if (
        length(z) == 0L ||
        sum(z^2) <= 0
    ) {
        return(NA_real_)
    }

    (
        sum(z)^2
    ) /
        sum(z^2)
}

extract_site_code <- function(x) {

    z <- toupper(
        safe_chr(x)
    )

    out <- rep(
        NA_character_,
        length(z)
    )

    hit <- grepl(
        "C[0-9]{3}",
        z,
        perl = TRUE
    )

    if (any(hit)) {

        temp <- regmatches(
            z[hit],
            gregexpr(
                "C[0-9]{3}",
                z[hit],
                perl = TRUE
            )
        )

        out[hit] <- vapply(
            temp,
            function(v) {

                if (
                    length(v) == 0L
                ) {
                    NA_character_
                } else {
                    tail(
                        v,
                        1L
                    )
                }
            },
            character(1)
        )
    }

    out
}


# ============================================================
# 2. Read formal clean CSV and reproduce cohort anchors
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

if (
    length(
        missing_main
    ) > 0L
) {
    stop(
        "Formal MSK CSV missing required columns:\n",
        paste(
            missing_main,
            collapse = "\n"
        )
    )
}

n_crc <- nrow(main)

# B1-14 anchor definition:
# the 458-patient first-relevant-biologic cohort is specifically
# anti-EGFR or bevacizumab, NOT every eligible patient with a nonmissing
# FIRST_BIOLOGIC_R value.

n_eligible <- main[
    OLD_ELIGIBLE_803_R == 1L,
    .N
]

n_biologic <- main[
    OLD_ELIGIBLE_803_R == 1L &
        FIRST_BIOLOGIC_R %chin%
            c(
                "anti-EGFR",
                "bevacizumab"
            ),
    .N
]

n_anti_all <- main[
    OLD_ELIGIBLE_803_R == 1L &
        FIRST_BIOLOGIC_R ==
            "anti-EGFR",
    .N
]

n_bev_all <- main[
    OLD_ELIGIBLE_803_R == 1L &
        FIRST_BIOLOGIC_R ==
            "bevacizumab",
    .N
]

if (
    n_crc != 5543L ||
    n_eligible != 803L ||
    n_biologic != 458L ||
    n_anti_all != 230L ||
    n_bev_all != 228L
) {
    stop(
        paste0(
            "Formal clean CSV cohort anchors failed.\n",
            "Observed CRC/eligible/biologic/anti/BEV = ",
            n_crc, "/",
            n_eligible, "/",
            n_biologic, "/",
            n_anti_all, "/",
            n_bev_all,
            "\nExpected 5543/803/458/230/228."
        )
    )
}

cohort <- main[
    OLD_ELIGIBLE_803_R == 1L &
        FIRST_BIOLOGIC_R %chin%
            c(
                "anti-EGFR",
                "bevacizumab"
            ) &
        is.finite(
            SEQUENCING_FIRST_DAY_R
        ) &
        is.finite(
            FIRST_BIOLOGIC_DAY_R
        ) &
        SEQUENCING_FIRST_DAY_R <=
            FIRST_BIOLOGIC_DAY_R
]

n_anti_seqpre <- cohort[
    FIRST_BIOLOGIC_R ==
        "anti-EGFR",
    .N
]

n_bev_seqpre <- cohort[
    FIRST_BIOLOGIC_R ==
        "bevacizumab",
    .N
]

if (
    n_anti_seqpre != 191L ||
    n_bev_seqpre != 137L ||
    nrow(cohort) != 328L
) {
    stop(
        paste0(
            "Sequencing-pre-biologic anchor failed.\n",
            "Observed anti/BEV/combined = ",
            n_anti_seqpre, "/",
            n_bev_seqpre, "/",
            nrow(cohort),
            "\nExpected 191/137/328."
        )
    )
}

cohort <- cohort[
    ,
    .(
        PATIENT_ID,
        SAMPLE_ID,
        FIRST_BIOLOGIC_R,
        FIRST_BIOLOGIC_DAY_R,
        BYPASS_ORIGINAL_R,
        SEQUENCING_FIRST_DAY_R,
        SPECIMEN_ACQUISITION_FIRST_DAY_R,
        GENDER,
        PRIMARY_SITE
    )
]

cohort[
    ,
    TREATMENT_R :=
        factor(
            FIRST_BIOLOGIC_R,
            levels = c(
                "bevacizumab",
                "anti-EGFR"
            )
        )
]

cat(
    "Sequencing-pre-biologic cohorts: anti-EGFR=",
    n_anti_seqpre,
    "; bevacizumab=",
    n_bev_seqpre,
    "\n",
    sep = ""
)


# ============================================================
# 3. Extract raw MSK files
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_14_"
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
        basename(
            all_files
        ) == filename
    ]

    if (
        length(hit) == 0L
    ) {
        stop(
            "Missing raw MSK file: ",
            filename
        )
    }

    hit[1]
}


# ============================================================
# 4. Sample metadata
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
    SAMPLE_ID %chin%
        cohort$SAMPLE_ID,
    .(
        SAMPLE_ID,
        PATIENT_ID,
        GENE_PANEL,
        SAMPLE_TYPE,
        METASTATIC_SITE,
        PRIMARY_SITE_RAW =
            PRIMARY_SITE,
        TMB_NONSYNONYMOUS
    )
]

if (
    nrow(
        sample_primary
    ) != 328L
) {
    stop(
        "Sample metadata map failed. Expected 328 rows; got ",
        nrow(
            sample_primary
        )
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
# 5. Exact B1-06 mutation definitions
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
    PATIENT_ID :=
        sample_map$PATIENT_ID[
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

erbb2_activating <- c(
    "p.S310F",
    "p.S310Y",
    "p.L755S",
    "p.V777L",
    "p.V842I",
    "p.L866M"
)

pik3ca_ex20 <- c(
    "p.M1043I",
    "p.H1047R",
    "p.H1047L",
    "p.H1047Y",
    "p.G1049R"
)

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

braf_class3 <- c(
    "p.D594G",
    "p.D594N",
    "p.N581S"
)

mut_goi[
    ,
    ERBB2_ACTIVATING_R :=
        as.integer(
            Hugo_Symbol ==
                "ERBB2" &
            HGVSp_Short %chin%
                erbb2_activating
        )
]

mut_goi[
    ,
    PIK3CA_EX20_R :=
        as.integer(
            Hugo_Symbol ==
                "PIK3CA" &
            HGVSp_Short %chin%
                pik3ca_ex20
        )
]

mut_goi[
    ,
    PIK3CA_EX9_R :=
        as.integer(
            Hugo_Symbol ==
                "PIK3CA" &
            HGVSp_Short %chin%
                pik3ca_ex9
        )
]

mut_goi[
    ,
    PTEN_LOF_R :=
        as.integer(
            Hugo_Symbol ==
                "PTEN" &
            Variant_Classification %chin%
                lof_classes
        )
]

mut_goi[
    ,
    AKT1_E17K_R :=
        as.integer(
            Hugo_Symbol ==
                "AKT1" &
            HGVSp_Short ==
                "p.E17K"
        )
]

mut_goi[
    ,
    MAP2K1_HOTSPOT_R :=
        as.integer(
            Hugo_Symbol ==
                "MAP2K1" &
            HGVSp_Short %chin%
                map2k1_hotspot
        )
]

mut_goi[
    ,
    MAP2K1_D67E_R :=
        as.integer(
            Hugo_Symbol ==
                "MAP2K1" &
            HGVSp_Short ==
                "p.D67E"
        )
]

mut_goi[
    ,
    NF1_LOF_R :=
        as.integer(
            Hugo_Symbol ==
                "NF1" &
            Variant_Classification %chin%
                lof_classes
        )
]

mut_goi[
    ,
    BRAF_CLASS3_R :=
        as.integer(
            Hugo_Symbol ==
                "BRAF" &
            HGVSp_Short %chin%
                braf_class3
        )
]

mut_goi[
    ,
    EGFR_MUTATION_DESCRIPTIVE_R :=
        as.integer(
            Hugo_Symbol ==
                "EGFR"
        )
]

mut_goi[
    ,
    GNAS_MUTATION_DESCRIPTIVE_R :=
        as.integer(
            Hugo_Symbol ==
                "GNAS"
        )
]

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
# 6. Exact B1-06 CNA definitions
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

if (
    length(
        usable_samples
    ) != 328L
) {
    stop(
        paste0(
            "CNA ascertainment is not complete across the combined cohort.\n",
            "Expected 328 samples in the CNA matrix; found ",
            length(
                usable_samples
            ),
            ".\nDo not interpret comparator genomic negativity until reviewed."
        )
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
        Hugo_Symbol ==
            gene
    ]

    if (
        nrow(row_g) == 0L
    ) {
        return(
            data.table(
                SAMPLE_ID =
                    usable_samples,
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
        SAMPLE_ID =
            usable_samples,
        GENE = gene,
        AMP =
            as.integer(
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
    PATIENT_ID :=
        sample_map$PATIENT_ID[
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
# 7. Exact B1-06 canonical fusion definition
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
    PATIENT_ID :=
        sample_map$PATIENT_ID[
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

if (
    nrow(sv) > 0L
) {

    sv[
        ,
        CANONICAL_PRESSING_FUSION_R :=
            as.integer(
                (
                    Site1_Hugo_Symbol %chin%
                        canonical_fusion_genes |
                    Site2_Hugo_Symbol %chin%
                        canonical_fusion_genes
                ) &
                grepl(
                    "fusion",
                    paste(
                        safe_chr(
                            Event_Info
                        ),
                        safe_chr(
                            Annotation
                        )
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

} else {

    fusion_patient <- data.table(
        PATIENT_ID =
            character(),
        CANONICAL_PRESSING_FUSION_R =
            integer()
    )
}


# ============================================================
# 8. Outcome-blind treatment context, exact B1-06 architecture
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
    START_DATE_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

treatment[
    ,
    STOP_DATE_R :=
        suppressWarnings(
            as.numeric(
                STOP_DATE
            )
        )
]

treatment[
    ,
    AGENT_UPPER_R :=
        toupper(
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

treatment[
    ,
    IS_BEV_R :=
        AGENT_UPPER_R ==
        "BEVACIZUMAB"
]

build_tx_context <- function(
    pid,
    index_day,
    arm
) {

    p <- treatment[
        PATIENT_ID == pid &
            is.finite(
                START_DATE_R
            )
    ]

    if (
        nrow(p) == 0L
    ) {
        return(
            data.table(
                PATIENT_ID = pid,
                INDEX_BIOLOGIC_AGENT_R =
                    NA_character_,
                PRIOR_HEAVY_R =
                    NA_integer_,
                PRIOR_CYTOTOXIC_PATTERN_R =
                    NA_character_,
                PRIOR_MAJOR_N_R =
                    NA_integer_,
                ACTIVE_BACKBONE_AT_INDEX_R =
                    NA_character_,
                ACTIVE_IRI_R =
                    NA_integer_,
                ACTIVE_OX_R =
                    NA_integer_
            )
        )
    }

    index_rows <- if (
        arm ==
            "anti-EGFR"
    ) {
        p[
            IS_EGFR_R == TRUE &
                START_DATE_R ==
                    index_day
        ]
    } else {
        p[
            IS_BEV_R == TRUE &
                START_DATE_R ==
                    index_day
        ]
    }

    index_agent <- collapse_unique(
        index_rows$AGENT
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
        prior_ox &
            prior_iri,
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

    prior_major_n <- as.integer(
        prior_ox
    ) +
        as.integer(
            prior_iri
        )

    prior_heavy <- as.integer(
        prior_major_n ==
            2L
    )

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
        has_fp &
            has_ox &
            has_iri,
        "triplet/FOLFOXIRI-like",

        has_fp &
            has_ox &
            !has_iri,
        "FOLFOX/CAPOX-like",

        has_fp &
            has_iri &
            !has_ox,
        "FOLFIRI/CAPIRI-like",

        has_fp &
            !has_ox &
            !has_iri,
        "fluoropyrimidine-only",

        !has_fp &
            !has_ox &
            has_iri,
        "irinotecan-only",

        !has_fp &
            has_ox &
            !has_iri,
        "oxaliplatin-only",

        default =
            "no major cytotoxic active"
    )

    data.table(
        PATIENT_ID = pid,
        INDEX_BIOLOGIC_AGENT_R =
            index_agent,
        PRIOR_HEAVY_R =
            prior_heavy,
        PRIOR_CYTOTOXIC_PATTERN_R =
            prior_pattern,
        PRIOR_MAJOR_N_R =
            prior_major_n,
        ACTIVE_BACKBONE_AT_INDEX_R =
            active_backbone,
        ACTIVE_IRI_R =
            as.integer(
                has_iri
            ),
        ACTIVE_OX_R =
            as.integer(
                has_ox
            )
    )
}

tx_context <- rbindlist(
    Map(
        build_tx_context,
        cohort$PATIENT_ID,
        cohort$FIRST_BIOLOGIC_DAY_R,
        cohort$FIRST_BIOLOGIC_R
    ),
    fill = TRUE
)


# ============================================================
# 9. Diagnosis: explicit Stage IV + frozen left subsite
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
    START_DATE_NUM_R :=
        suppressWarnings(
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
    IS_CRC_DX_R :=
        grepl(
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
                safe_chr(
                    SUMMARY
                ),
                ignore.case = TRUE
            )
        )
]

stage4 <- diagnosis[
    IS_STAGE4_CRC_R ==
        TRUE,
    .(
        EXPLICIT_STAGE4_CRC_DAY_R =
            min_or_na(
                START_DATE_NUM_R
            )
    ),
    by = PATIENT_ID
]

left_colon_codes <- c(
    "C185",
    "C186",
    "C187"
)

rectal_codes <- c(
    "C199",
    "C209",
    "C218"
)

left_all_codes <- c(
    left_colon_codes,
    rectal_codes
)

right_codes <- c(
    "C180",
    "C181",
    "C182",
    "C183",
    "C184"
)

derive_dx_patient <- function(pid) {

    p <- diagnosis[
        PATIENT_ID == pid
    ]

    codes <- unique(
        p$SITE_CODE_R[
            !is.na(
                p$SITE_CODE_R
            )
        ]
    )

    has_left <- any(
        codes %chin%
            left_all_codes
    )

    has_right <- any(
        codes %chin%
            right_codes
    )

    dx_side <- if (
        has_left &&
        !has_right
    ) {
        "left"
    } else if (
        has_right &&
        !has_left
    ) {
        "right"
    } else if (
        has_left &&
        has_right
    ) {
        "conflict"
    } else {
        "unknown"
    }

    p_crc_left <- p[
        SITE_CODE_R %chin%
            left_all_codes &
        is.finite(
            START_DATE_NUM_R
        )
    ]

    if (
        nrow(
            p_crc_left
        ) > 0L
    ) {

        setorder(
            p_crc_left,
            START_DATE_NUM_R
        )

        first_code <-
            p_crc_left$SITE_CODE_R[1]

        dx_subsite <- if (
            first_code %chin%
                rectal_codes
        ) {
            "rectum/rectosigmoid"
        } else {
            "left colon"
        }

    } else {

        dx_subsite <-
            NA_character_
    }

    data.table(
        PATIENT_ID = pid,
        DX_SIDE_R =
            dx_side,
        DX_SUBSITE_R =
            dx_subsite
    )
}

dx_patient <- rbindlist(
    lapply(
        cohort$PATIENT_ID,
        derive_dx_patient
    )
)


# ============================================================
# 10. ECOG and recent metastatic-site burden
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

performance <- performance[
    PATIENT_ID %chin%
        cohort$PATIENT_ID
]

performance[
    ,
    START_DATE_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

performance[
    ,
    ECOG_R :=
        suppressWarnings(
            as.numeric(
                ECOG
            )
        )
]

build_ecog <- function(
    pid,
    index_day
) {

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

    if (
        nrow(p) == 0L
    ) {
        return(
            data.table(
                PATIENT_ID = pid,
                ECOG_LATEST_PRE_R =
                    NA_real_,
                ECOG_LAG_DAYS_R =
                    NA_real_,
                ECOG_WITHIN_90D_R =
                    NA_real_
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
                lag <=
                    90,
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
    START_DATE_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

build_site_burden <- function(
    pid,
    index_day
) {

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

    if (
        nrow(p) == 0L
    ) {
        return(
            data.table(
                PATIENT_ID = pid,
                MET_SITE_COUNT_90D_R =
                    NA_integer_,
                MET_SITE_TYPES_90D_R =
                    NA_character_
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
            if (
                length(sites) ==
                    0L
            ) {
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
# 11. Assemble outcome-blind comparator audit
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

for (obj0 in list(
    mutation_patient,
    amp_patient,
    fusion_patient,
    tx_context,
    stage4,
    dx_patient,
    ecog_context,
    site_context
)) {

    audit <- merge(
        audit,
        obj0,
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

    if (
        !v %in%
            names(audit)
    ) {
        audit[
            ,
            (v) :=
                0L
        ]
    }

    audit[
        is.na(
            get(v)
        ),
        (v) :=
            0L
    ]
}

audit[
    ,
    BYPASS_PRESSING_STRICT_R :=
        as.integer(
            ERBB2_AMP_R == 1L |
            ERBB2_ACTIVATING_R == 1L |
            MET_AMP_R == 1L |
            PIK3CA_EX20_R == 1L |
            PTEN_LOF_R == 1L |
            AKT1_E17K_R == 1L |
            CANONICAL_PRESSING_FUSION_R == 1L
        )
]

audit[
    ,
    BYPASS_HIGH_CONFIDENCE_R :=
        as.integer(
            BYPASS_PRESSING_STRICT_R == 1L |
            KRAS_AMP_R == 1L |
            MAP2K1_HOTSPOT_R == 1L |
            NF1_LOF_R == 1L
        )
]

audit[
    ,
    BYPASS_HIGH_CONFIDENCE_NO_D67E_R :=
        as.integer(
            BYPASS_PRESSING_STRICT_R == 1L |
            KRAS_AMP_R == 1L |
            (
                MAP2K1_HOTSPOT_R == 1L &
                MAP2K1_D67E_R == 0L
            ) |
            NF1_LOF_R == 1L
        )
]

audit[
    ,
    BYPASS_ORIGINAL_BROAD_R :=
        BYPASS_ORIGINAL_R
]

audit[
    ,
    STAGE4_TO_BIOLOGIC_DAYS_R :=
        FIRST_BIOLOGIC_DAY_R -
        EXPLICIT_STAGE4_CRC_DAY_R
]

audit[
    ,
    LOG_STAGE4_TO_BIOLOGIC_R :=
        fifelse(
            is.finite(
                STAGE4_TO_BIOLOGIC_DAYS_R
            ) &
            STAGE4_TO_BIOLOGIC_DAYS_R >=
                0,
            log1p(
                STAGE4_TO_BIOLOGIC_DAYS_R
            ),
            NA_real_
        )
]

audit[
    ,
    SPECIMEN_TO_BIOLOGIC_DAYS_R :=
        FIRST_BIOLOGIC_DAY_R -
        SPECIMEN_ACQUISITION_FIRST_DAY_R
]

audit[
    ,
    SEQUENCING_TO_BIOLOGIC_DAYS_R :=
        FIRST_BIOLOGIC_DAY_R -
        SEQUENCING_FIRST_DAY_R
]

audit[
    ,
    LOG_SPECIMEN_TO_BIOLOGIC_R :=
        fifelse(
            is.finite(
                SPECIMEN_TO_BIOLOGIC_DAYS_R
            ) &
            SPECIMEN_TO_BIOLOGIC_DAYS_R >=
                0,
            log1p(
                SPECIMEN_TO_BIOLOGIC_DAYS_R
            ),
            NA_real_
        )
]

audit[
    ,
    PANEL_OLD_R :=
        as.integer(
            GENE_PANEL %chin%
                c(
                    "IMPACT341",
                    "IMPACT410"
                )
        )
]

sample_type_upper <- toupper(
    trimws(
        safe_chr(
            audit$SAMPLE_TYPE
        )
    )
)

audit[
    ,
    SAMPLE_METASTATIC_R :=
        fifelse(
            nzchar(
                sample_type_upper
            ),
            as.integer(
                grepl(
                    "METAST",
                    sample_type_upper
                )
            ),
            NA_integer_
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
    STRUCTURED_SIDE_R :=
        fcase(
            primary_site_lower %chin%
                c(
                    "rectum",
                    "sigmoid colon",
                    "rectosigmoid colon",
                    "descending colon",
                    "left colon"
                ),
            "left",

            primary_site_lower %chin%
                c(
                    "cecum",
                    "ascending colon",
                    "transverse colon",
                    "hepatic flexure"
                ),
            "right",

            default =
                "unknown"
        )
]

audit[
    ,
    STRUCTURED_SUBSITE_R :=
        fcase(
            primary_site_lower %chin%
                c(
                    "rectum",
                    "rectosigmoid colon"
                ),
            "rectum/rectosigmoid",

            primary_site_lower %chin%
                c(
                    "sigmoid colon",
                    "descending colon",
                    "left colon"
                ),
            "left colon",

            default =
                NA_character_
        )
]

audit[
    ,
    SUBSITE_FROZEN_R :=
        fifelse(
            !is.na(
                STRUCTURED_SUBSITE_R
            ),
            STRUCTURED_SUBSITE_R,
            DX_SUBSITE_R
        )
]

audit[
    ,
    SIDE_CONFLICT_R :=
        as.integer(
            STRUCTURED_SIDE_R %chin%
                c(
                    "left",
                    "right"
                ) &
            DX_SIDE_R %chin%
                c(
                    "left",
                    "right"
                ) &
            STRUCTURED_SIDE_R !=
                DX_SIDE_R
        )
]

# Validate subsite algorithm against locked anti-EGFR reference.
anti_subsite_check <- audit[
    FIRST_BIOLOGIC_R ==
        "anti-EGFR",
    .(
        LEFT_COLON =
            sum(
                SUBSITE_FROZEN_R ==
                    "left colon",
                na.rm = TRUE
            ),
        RECTUM =
            sum(
                SUBSITE_FROZEN_R ==
                    "rectum/rectosigmoid",
                na.rm = TRUE
            ),
        MISSING =
            sum(
                is.na(
                    SUBSITE_FROZEN_R
                )
            ),
        CONFLICT =
            sum(
                SIDE_CONFLICT_R ==
                    1L,
                na.rm = TRUE
            )
    )
]

if (
    anti_subsite_check$LEFT_COLON !=
        96L ||
    anti_subsite_check$RECTUM !=
        95L ||
    anti_subsite_check$MISSING !=
        0L ||
    anti_subsite_check$CONFLICT !=
        1L
) {
    stop(
        paste0(
            "Frozen subsite algorithm failed anti-EGFR validation.\n",
            "Observed left/rectum/missing/conflict = ",
            anti_subsite_check$LEFT_COLON, "/",
            anti_subsite_check$RECTUM, "/",
            anti_subsite_check$MISSING, "/",
            anti_subsite_check$CONFLICT,
            "\nExpected 96/95/0/1."
        )
    )
}


# ============================================================
# 12. HARD anti-EGFR genomic validation BEFORE interpreting BEV
# ============================================================

anti <- audit[
    FIRST_BIOLOGIC_R ==
        "anti-EGFR"
]

bev <- audit[
    FIRST_BIOLOGIC_R ==
        "bevacizumab"
]

component_vars <- c(
    "ERBB2_AMP_R",
    "ERBB2_ACTIVATING_R",
    "MET_AMP_R",
    "PIK3CA_EX20_R",
    "PTEN_LOF_R",
    "AKT1_E17K_R",
    "CANONICAL_PRESSING_FUSION_R",
    "KRAS_AMP_R",
    "MAP2K1_HOTSPOT_R",
    "NF1_LOF_R"
)

expected_component_counts <- c(
    ERBB2_AMP_R = 8L,
    ERBB2_ACTIVATING_R = 3L,
    MET_AMP_R = 1L,
    PIK3CA_EX20_R = 3L,
    PTEN_LOF_R = 2L,
    AKT1_E17K_R = 0L,
    CANONICAL_PRESSING_FUSION_R = 0L,
    KRAS_AMP_R = 3L,
    MAP2K1_HOTSPOT_R = 2L,
    NF1_LOF_R = 3L
)

observed_component_counts <- vapply(
    component_vars,
    function(v) {
        anti[
            get(v) ==
                1L,
            .N
        ]
    },
    integer(1)
)

anti_anchor_vector <- c(
    N =
        nrow(
            anti
        ),
    STRICT =
        anti[
            BYPASS_PRESSING_STRICT_R ==
                1L,
            .N
        ],
    HIGH_CONF =
        anti[
            BYPASS_HIGH_CONFIDENCE_R ==
                1L,
            .N
        ],
    NO_D67E =
        anti[
            BYPASS_HIGH_CONFIDENCE_NO_D67E_R ==
                1L,
            .N
        ],
    BROAD =
        anti[
            BYPASS_ORIGINAL_BROAD_R ==
                1L,
            .N
        ],
    BRAF3 =
        anti[
            BRAF_CLASS3_R ==
                1L,
            .N
        ],
    PIK3CA_EX9 =
        anti[
            PIK3CA_EX9_R ==
                1L,
            .N
        ]
)

expected_anchor_vector <- c(
    N = 191L,
    STRICT = 15L,
    HIGH_CONF = 23L,
    NO_D67E = 22L,
    BROAD = 35L,
    BRAF3 = 7L,
    PIK3CA_EX9 = 6L
)

if (
    !identical(
        as.integer(
            anti_anchor_vector
        ),
        as.integer(
            expected_anchor_vector
        )
    ) ||
    !identical(
        as.integer(
            observed_component_counts
        ),
        as.integer(
            expected_component_counts[
                component_vars
            ]
        )
    )
) {
    stop(
        paste0(
            "B1-14 STOP: exact B1-06 anti-EGFR genomic validation failed.\n",
            "Do NOT interpret bevacizumab classification.\n",
            "Observed anchors:\n",
            paste(
                names(
                    anti_anchor_vector
                ),
                anti_anchor_vector,
                sep = "=",
                collapse = ", "
            ),
            "\nObserved component counts:\n",
            paste(
                component_vars,
                observed_component_counts,
                sep = "=",
                collapse = ", "
            )
        )
    )
}

cat(
    "Exact B1-06 anti-EGFR genomic validation: PASS\n"
)


# ============================================================
# 13. Genomic prevalence by treatment arm
# ============================================================

component_labels <- c(
    ERBB2_AMP_R =
        "ERBB2 amplification",
    ERBB2_ACTIVATING_R =
        "ERBB2 activating mutation",
    MET_AMP_R =
        "MET amplification",
    PIK3CA_EX20_R =
        "PIK3CA exon 20",
    PTEN_LOF_R =
        "PTEN loss-of-function",
    AKT1_E17K_R =
        "AKT1 E17K",
    CANONICAL_PRESSING_FUSION_R =
        "Canonical ALK/ROS1/NTRK/RET fusion",
    KRAS_AMP_R =
        "KRAS amplification",
    MAP2K1_HOTSPOT_R =
        "MAP2K1 hotspot",
    NF1_LOF_R =
        "NF1 loss-of-function"
)

genomic_summary <- rbindlist(
    lapply(
        c(
            "anti-EGFR",
            "bevacizumab"
        ),
        function(arm) {

            ds <- audit[
                FIRST_BIOLOGIC_R ==
                    arm
            ]

            data.table(
                Treatment = arm,
                N = nrow(ds),
                `High-confidence+, n` =
                    ds[
                        BYPASS_HIGH_CONFIDENCE_R ==
                            1L,
                        .N
                    ],
                `High-confidence+, %` =
                    100 *
                    mean(
                        ds$BYPASS_HIGH_CONFIDENCE_R ==
                            1L
                    ),
                `Strict+, n` =
                    ds[
                        BYPASS_PRESSING_STRICT_R ==
                            1L,
                        .N
                    ],
                `No-D67E+, n` =
                    ds[
                        BYPASS_HIGH_CONFIDENCE_NO_D67E_R ==
                            1L,
                        .N
                    ],
                `Broad+, n` =
                    ds[
                        BYPASS_ORIGINAL_BROAD_R ==
                            1L,
                        .N
                    ],
                `BRAF class III, n` =
                    ds[
                        BRAF_CLASS3_R ==
                            1L,
                        .N
                    ],
                `PIK3CA exon9, n` =
                    ds[
                        PIK3CA_EX9_R ==
                            1L,
                        .N
                    ]
            )
        }
    ),
    fill = TRUE
)

component_summary <- rbindlist(
    lapply(
        c(
            "anti-EGFR",
            "bevacizumab"
        ),
        function(arm) {

            ds <- audit[
                FIRST_BIOLOGIC_R ==
                    arm
            ]

            rbindlist(
                lapply(
                    component_vars,
                    function(v) {

                        data.table(
                            Treatment =
                                arm,
                            Component =
                                unname(
                                    component_labels[
                                        v
                                    ]
                                ),
                            N =
                                ds[
                                    get(v) ==
                                        1L,
                                    .N
                                ],
                            Percent =
                                100 *
                                mean(
                                    ds[[v]] ==
                                        1L
                                )
                        )
                    }
                )
            )
        }
    )
)


# ============================================================
# 14. Treatment-context / positivity overlap audits
# ============================================================

audit[
    ,
    PRIOR_MAJOR_N_F_R :=
        factor(
            PRIOR_MAJOR_N_R,
            levels = 0:2,
            labels = c(
                "0",
                "1",
                "2"
            )
        )
]

audit[
    ,
    ACTIVE_IRI_F_R :=
        factor(
            ACTIVE_IRI_R,
            levels = c(
                0,
                1
            ),
            labels = c(
                "No",
                "Yes"
            )
        )
]

audit[
    ,
    ACTIVE_OX_F_R :=
        factor(
            ACTIVE_OX_R,
            levels = c(
                0,
                1
            ),
            labels = c(
                "No",
                "Yes"
            )
        )
]

audit[
    ,
    PANEL_OLD_F_R :=
        factor(
            PANEL_OLD_R,
            levels = c(
                0,
                1
            ),
            labels = c(
                "Newer panel",
                "IMPACT341/410"
            )
        )
]

audit[
    ,
    SAMPLE_METASTATIC_F_R :=
        factor(
            SAMPLE_METASTATIC_R,
            levels = c(
                0,
                1
            ),
            labels = c(
                "Other/nonmetastatic",
                "Metastatic"
            )
        )
]

audit[
    ,
    SUBSITE_F_R :=
        factor(
            SUBSITE_FROZEN_R,
            levels = c(
                "left colon",
                "rectum/rectosigmoid"
            )
        )
]

audit[
    ,
    GENDER_F_R :=
        factor(
            GENDER,
            levels = c(
                "Female",
                "Male"
            )
        )
]

audit[
    ,
    STAGE4_TO_BIOLOGIC_BAND_R :=
        fcase(
            !is.finite(
                STAGE4_TO_BIOLOGIC_DAYS_R
            ),
            "missing",

            STAGE4_TO_BIOLOGIC_DAYS_R <=
                60,
            "<=60 d",

            STAGE4_TO_BIOLOGIC_DAYS_R <=
                180,
            "61-180 d",

            STAGE4_TO_BIOLOGIC_DAYS_R <=
                365,
            "181-365 d",

            default =
                ">365 d"
        )
]

prior_cells <- audit[
    ,
    .N,
    by = .(
        PRIOR_MAJOR_N =
            PRIOR_MAJOR_N_F_R,
        Treatment =
            FIRST_BIOLOGIC_R,
        Bypass =
            BYPASS_HIGH_CONFIDENCE_R
    )
][
    order(
        PRIOR_MAJOR_N,
        Treatment,
        Bypass
    )
]

backbone_cells <- audit[
    ,
    .N,
    by = .(
        ACTIVE_BACKBONE_AT_INDEX_R,
        Treatment =
            FIRST_BIOLOGIC_R,
        Bypass =
            BYPASS_HIGH_CONFIDENCE_R
    )
][
    order(
        ACTIVE_BACKBONE_AT_INDEX_R,
        Treatment,
        Bypass
    )
]

timing_cells <- audit[
    ,
    .N,
    by = .(
        STAGE4_TO_BIOLOGIC_BAND_R,
        Treatment =
            FIRST_BIOLOGIC_R,
        Bypass =
            BYPASS_HIGH_CONFIDENCE_R
    )
][
    order(
        STAGE4_TO_BIOLOGIC_BAND_R,
        Treatment,
        Bypass
    )
]

context_summary <- audit[
    ,
    .(
        N = .N,
        BYPASS_POS =
            sum(
                BYPASS_HIGH_CONFIDENCE_R ==
                    1L
            ),
        PRIOR_0 =
            sum(
                PRIOR_MAJOR_N_R ==
                    0L,
                na.rm = TRUE
            ),
        PRIOR_1 =
            sum(
                PRIOR_MAJOR_N_R ==
                    1L,
                na.rm = TRUE
            ),
        PRIOR_2 =
            sum(
                PRIOR_MAJOR_N_R ==
                    2L,
                na.rm = TRUE
            ),
        ACTIVE_IRI =
            sum(
                ACTIVE_IRI_R ==
                    1L,
                na.rm = TRUE
            ),
        ACTIVE_OX =
            sum(
                ACTIVE_OX_R ==
                    1L,
                na.rm = TRUE
            ),
        STAGE4_TO_INDEX =
            median_iqr(
                STAGE4_TO_BIOLOGIC_DAYS_R
            ),
        SPECIMEN_TO_INDEX =
            median_iqr(
                SPECIMEN_TO_BIOLOGIC_DAYS_R
            ),
        SEQUENCING_TO_INDEX =
            median_iqr(
                SEQUENCING_TO_BIOLOGIC_DAYS_R
            ),
        MET_SITE_COUNT =
            median_iqr(
                MET_SITE_COUNT_90D_R
            ),
        ECOG_90D_AVAILABLE =
            sum(
                is.finite(
                    ECOG_WITHIN_90D_R
                )
            )
    ),
    by =
        FIRST_BIOLOGIC_R
]


# ============================================================
# 15. Outcome-blind propensity-overlap diagnostic
# ============================================================

ps_formula <- as.formula(
    paste0(
        "I(FIRST_BIOLOGIC_R == 'anti-EGFR') ~ ",
        "GENDER_F_R + SUBSITE_F_R + ",
        "LOG_STAGE4_TO_BIOLOGIC_R + ",
        "PRIOR_MAJOR_N_F_R + ",
        "ACTIVE_IRI_F_R + ACTIVE_OX_F_R + ",
        "PANEL_OLD_F_R + SAMPLE_METASTATIC_F_R + ",
        "LOG_SPECIMEN_TO_BIOLOGIC_R + ",
        "MET_SITE_COUNT_90D_R"
    )
)

ps_vars <- all.vars(
    ps_formula
)

# all.vars() includes FIRST_BIOLOGIC_R via I()
ps_covars <- setdiff(
    ps_vars,
    "FIRST_BIOLOGIC_R"
)

ps_complete <- complete.cases(
    audit[
        ,
        ..ps_covars
    ]
)

audit[
    ,
    PS_COMPLETE_R :=
        ps_complete
]

ps_coverage <- audit[
    ,
    .(
        N = .N,
        PS_COMPLETE_N =
            sum(
                PS_COMPLETE_R
            ),
        PS_COMPLETE_PCT =
            100 *
            mean(
                PS_COMPLETE_R
            ),
        BYPASS_POS =
            sum(
                BYPASS_HIGH_CONFIDENCE_R ==
                    1L
            ),
        BYPASS_POS_PS_COMPLETE =
            sum(
                BYPASS_HIGH_CONFIDENCE_R ==
                    1L &
                PS_COMPLETE_R
            )
    ),
    by =
        FIRST_BIOLOGIC_R
]

dps <- audit[
    PS_COMPLETE_R ==
        TRUE
]

ps_fit <- tryCatch(
    suppressWarnings(
        glm(
            ps_formula,
            data = dps,
            family = binomial()
        )
    ),
    error = function(e) NULL
)

ps_estimable <- (
    !is.null(
        ps_fit
    ) &&
    isTRUE(
        ps_fit$converged
    ) &&
    all(
        is.finite(
            coef(
                ps_fit
            )
        )
    )
)

common_support_exists <- FALSE
common_lower <- NA_real_
common_upper <- NA_real_

if (
    ps_estimable
) {

    dps[
        ,
        PS_ANTI_R :=
            as.numeric(
                predict(
                    ps_fit,
                    type = "response"
                )
            )
    ]

    if (
        any(
            !is.finite(
                dps$PS_ANTI_R
            )
        )
    ) {
        ps_estimable <- FALSE
    }
}

if (
    ps_estimable
) {

    anti_ps <- dps[
        FIRST_BIOLOGIC_R ==
            "anti-EGFR",
        PS_ANTI_R
    ]

    bev_ps <- dps[
        FIRST_BIOLOGIC_R ==
            "bevacizumab",
        PS_ANTI_R
    ]

    common_lower <- max(
        min(
            anti_ps
        ),
        min(
            bev_ps
        )
    )

    common_upper <- min(
        max(
            anti_ps
        ),
        max(
            bev_ps
        )
    )

    common_support_exists <-
        is.finite(
            common_lower
        ) &&
        is.finite(
            common_upper
        ) &&
        common_lower <
            common_upper

    if (
        common_support_exists
    ) {

        dps[
            ,
            IN_COMMON_SUPPORT_R :=
                PS_ANTI_R >=
                    common_lower &
                PS_ANTI_R <=
                    common_upper
        ]

        dps[
            ,
            OVERLAP_WEIGHT_R :=
                fifelse(
                    FIRST_BIOLOGIC_R ==
                        "anti-EGFR",
                    1 -
                        PS_ANTI_R,
                    PS_ANTI_R
                )
        ]
    }
}

if (
    !ps_estimable ||
    !common_support_exists
) {

    overlap_metrics <- data.table(
        Treatment =
            c(
                "anti-EGFR",
                "bevacizumab"
            ),
        PS_COMPLETE_N =
            NA_integer_,
        COMMON_SUPPORT_N =
            NA_integer_,
        COMMON_SUPPORT_PCT =
            NA_real_,
        BYPASS_POS_PS_N =
            NA_integer_,
        BYPASS_POS_IN_SUPPORT =
            NA_integer_,
        BYPASS_POS_SUPPORT_PCT =
            NA_real_,
        OVERLAP_ESS =
            NA_real_,
        BYPASS_POS_OVERLAP_ESS =
            NA_real_
    )

    ps_quantiles <- data.table()

} else {

    overlap_metrics <- dps[
        ,
        .(
            PS_COMPLETE_N =
                .N,
            COMMON_SUPPORT_N =
                sum(
                    IN_COMMON_SUPPORT_R
                ),
            COMMON_SUPPORT_PCT =
                100 *
                mean(
                    IN_COMMON_SUPPORT_R
                ),
            BYPASS_POS_PS_N =
                sum(
                    BYPASS_HIGH_CONFIDENCE_R ==
                        1L
                ),
            BYPASS_POS_IN_SUPPORT =
                sum(
                    BYPASS_HIGH_CONFIDENCE_R ==
                        1L &
                    IN_COMMON_SUPPORT_R
                ),
            BYPASS_POS_SUPPORT_PCT =
                100 *
                mean(
                    IN_COMMON_SUPPORT_R[
                        BYPASS_HIGH_CONFIDENCE_R ==
                            1L
                    ]
                ),
            OVERLAP_ESS =
                ess(
                    OVERLAP_WEIGHT_R
                ),
            BYPASS_POS_OVERLAP_ESS =
                ess(
                    OVERLAP_WEIGHT_R[
                        BYPASS_HIGH_CONFIDENCE_R ==
                            1L
                    ]
                )
        ),
        by =
            FIRST_BIOLOGIC_R
    ]

    setnames(
        overlap_metrics,
        "FIRST_BIOLOGIC_R",
        "Treatment"
    )

    ps_quantiles <- dps[
        ,
        {
            q <- quantile(
                PS_ANTI_R,
                probs = c(
                    0,
                    0.05,
                    0.25,
                    0.50,
                    0.75,
                    0.95,
                    1
                ),
                na.rm = TRUE,
                names = FALSE
            )

            .(
                Min = q[1],
                P05 = q[2],
                P25 = q[3],
                Median = q[4],
                P75 = q[5],
                P95 = q[6],
                Max = q[7]
            )
        },
        by =
            FIRST_BIOLOGIC_R
    ]

    setnames(
        ps_quantiles,
        "FIRST_BIOLOGIC_R",
        "Treatment"
    )
}


# ============================================================
# 16. Predeclared feasibility gate
# ============================================================

anti_pos <- anti[
    BYPASS_HIGH_CONFIDENCE_R ==
        1L,
    .N
]

anti_neg <- anti[
    BYPASS_HIGH_CONFIDENCE_R ==
        0L,
    .N
]

bev_pos <- bev[
    BYPASS_HIGH_CONFIDENCE_R ==
        1L,
    .N
]

bev_neg <- bev[
    BYPASS_HIGH_CONFIDENCE_R ==
        0L,
    .N
]

hard_counts_pass <- (
    anti_pos >= 10L &&
    bev_pos >= 10L &&
    anti_neg >= 50L &&
    bev_neg >= 50L
)

ps_coverage_pass <- FALSE
support_pass <- FALSE
positive_support_pass <- FALSE
ess_pass <- FALSE
positive_ess_pass <- FALSE

if (
    ps_estimable &&
    common_support_exists
) {

    cov_anti <- ps_coverage[
        FIRST_BIOLOGIC_R ==
            "anti-EGFR",
        PS_COMPLETE_PCT
    ]

    cov_bev <- ps_coverage[
        FIRST_BIOLOGIC_R ==
            "bevacizumab",
        PS_COMPLETE_PCT
    ]

    ps_coverage_pass <- (
        length(cov_anti) == 1L &&
        length(cov_bev) == 1L &&
        cov_anti >= 75 &&
        cov_bev >= 75
    )

    om_anti <- overlap_metrics[
        Treatment ==
            "anti-EGFR"
    ]

    om_bev <- overlap_metrics[
        Treatment ==
            "bevacizumab"
    ]

    support_pass <- (
        nrow(om_anti) == 1L &&
        nrow(om_bev) == 1L &&
        om_anti$COMMON_SUPPORT_PCT >=
            80 &&
        om_bev$COMMON_SUPPORT_PCT >=
            80
    )

    positive_support_pass <- (
        om_anti$BYPASS_POS_SUPPORT_PCT >=
            70 &&
        om_bev$BYPASS_POS_SUPPORT_PCT >=
            70
    )

    ess_pass <- (
        om_anti$OVERLAP_ESS >=
            50 &&
        om_bev$OVERLAP_ESS >=
            50
    )

    positive_ess_pass <- (
        om_anti$BYPASS_POS_OVERLAP_ESS >=
            8 &&
        om_bev$BYPASS_POS_OVERLAP_ESS >=
            8
    )
}

all_overlap_pass <- (
    ps_coverage_pass &&
    support_pass &&
    positive_support_pass &&
    ess_pass &&
    positive_ess_pass
)

if (
    !hard_counts_pass ||
    !ps_estimable ||
    !common_support_exists
) {

    gate_decision <-
        "NO-GO"

    gate_interpretation <- paste0(
        "Do NOT run treatment x bypass interaction. ",
        "Comparator feasibility is insufficient for a credible exploratory ",
        "active-comparator analysis. Report the feasibility limitation in ",
        "the main Discussion."
    )

} else if (
    all_overlap_pass
) {

    gate_decision <-
        "GO"

    gate_interpretation <- paste0(
        "An explicitly exploratory comparator analysis may be estimable. ",
        "This does not establish exchangeability or causal comparability. ",
        "Any next-stage analysis must retain an exploratory label and use ",
        "overlap-aware adjustment."
    )

} else {

    gate_decision <-
        "CONDITIONAL GO"

    gate_interpretation <- paste0(
        "Hard genomic cell counts are adequate and the PS model is estimable, ",
        "but at least one prespecified overlap threshold failed. ",
        "Do NOT automatically run an interaction. Review the overlap diagnostics ",
        "before deciding whether only a descriptive specificity analysis is defensible."
    )
}

gate_table <- data.table(
    Criterion = c(
        "High-confidence bypass+ >=10 in each arm",
        "Bypass- >=50 in each arm",
        "Outcome-blind PS model estimable",
        "PS complete-case coverage >=75% in each arm",
        "Common-support coverage >=80% in each arm",
        "Bypass+ common-support coverage >=70% in each arm",
        "Overlap-weight ESS >=50 in each arm",
        "Bypass+ overlap-weight ESS >=8 in each arm"
    ),
    Passed = c(
        anti_pos >=
            10L &&
            bev_pos >=
            10L,
        anti_neg >=
            50L &&
            bev_neg >=
            50L,
        ps_estimable,
        ps_coverage_pass,
        support_pass,
        positive_support_pass,
        ess_pass,
        positive_ess_pass
    )
)


# ============================================================
# 17. Save RDS
# ============================================================

saveRDS(
    list(
        analysis_status =
            "OUTCOME-BLIND COMPARATOR FEASIBILITY GATE ONLY",

        cohort_anchors =
            data.table(
                Metric = c(
                    "CRC",
                    "Eligible",
                    "First relevant biologic",
                    "All anti-EGFR first",
                    "All bevacizumab first",
                    "Sequencing-pre anti-EGFR",
                    "Sequencing-pre bevacizumab"
                ),
                N = c(
                    n_crc,
                    n_eligible,
                    n_biologic,
                    n_anti_all,
                    n_bev_all,
                    n_anti_seqpre,
                    n_bev_seqpre
                )
            ),

        patient_audit =
            audit,

        anti_validation =
            list(
                anchor_vector =
                    anti_anchor_vector,
                component_counts =
                    observed_component_counts,
                subsite_check =
                    anti_subsite_check
            ),

        genomic_summary =
            genomic_summary,

        component_summary =
            component_summary,

        context_summary =
            context_summary,

        prior_cells =
            prior_cells,

        backbone_cells =
            backbone_cells,

        timing_cells =
            timing_cells,

        ps_formula =
            ps_formula,

        ps_model =
            ps_fit,

        ps_coverage =
            ps_coverage,

        ps_common_support =
            c(
                lower =
                    common_lower,
                upper =
                    common_upper
            ),

        ps_quantiles =
            ps_quantiles,

        overlap_metrics =
            overlap_metrics,

        gate_table =
            gate_table,

        gate_decision =
            gate_decision,

        gate_interpretation =
            gate_interpretation
    ),
    result_file
)


# ============================================================
# 18. Audit TXT
# ============================================================

audit_lines <- c(
    "B1-14 BEVACIZUMAB COMPARATOR FEASIBILITY GATE",
    "=============================================",
    "",
    "STATUS",
    "------",
    "OUTCOME-BLIND comparator feasibility only.",
    "No rwPFS, OS, death, progression, HR, treatment effect, or interaction was analyzed.",
    "Frozen high-confidence genomic definition unchanged.",
    "",
    "COHORT ANCHORS",
    "--------------",
    paste0(
        "CRC = ",
        n_crc
    ),
    paste0(
        "Eligible = ",
        n_eligible
    ),
    paste0(
        "First relevant biologic cohort = ",
        n_biologic
    ),
    paste0(
        "All anti-EGFR / bevacizumab first = ",
        n_anti_all,
        " / ",
        n_bev_all
    ),
    paste0(
        "Sequencing-pre-biologic anti-EGFR / bevacizumab = ",
        n_anti_seqpre,
        " / ",
        n_bev_seqpre
    ),
    "",
    "EXACT B1-06 ANTI-EGFR GENOMIC VALIDATION",
    "----------------------------------------",
    "PASS",
    paste(
        names(
            anti_anchor_vector
        ),
        anti_anchor_vector,
        sep = "=",
        collapse = ", "
    ),
    paste(
        names(
            observed_component_counts
        ),
        observed_component_counts,
        sep = "=",
        collapse = ", "
    ),
    paste0(
        "Subsite left/rectum/missing/conflict = ",
        anti_subsite_check$LEFT_COLON, "/",
        anti_subsite_check$RECTUM, "/",
        anti_subsite_check$MISSING, "/",
        anti_subsite_check$CONFLICT
    ),
    "",
    "GENOMIC PREVALENCE BY ARM",
    "-------------------------",
    capture.output(
        print(
            genomic_summary
        )
    ),
    "",
    "COMPONENT COUNTS BY ARM",
    "-----------------------",
    capture.output(
        print(
            component_summary
        )
    ),
    "",
    "OUTCOME-BLIND TREATMENT CONTEXT",
    "-------------------------------",
    capture.output(
        print(
            context_summary
        )
    ),
    "",
    "TREATMENT x HIGH-CONFIDENCE BYPASS WITHIN PRIOR-MAJOR 0/1/2",
    "-----------------------------------------------------------",
    capture.output(
        print(
            prior_cells
        )
    ),
    "",
    "TREATMENT x HIGH-CONFIDENCE BYPASS WITHIN ACTIVE BACKBONE",
    "---------------------------------------------------------",
    capture.output(
        print(
            backbone_cells
        )
    ),
    "",
    "TREATMENT x HIGH-CONFIDENCE BYPASS WITHIN STAGE-IV TIMING BANDS",
    "--------------------------------------------------------------",
    capture.output(
        print(
            timing_cells
        )
    ),
    "",
    "PROPENSITY MODEL COVERAGE",
    "-------------------------",
    paste0(
        "PS model estimable: ",
        ps_estimable
    ),
    capture.output(
        print(
            ps_coverage
        )
    ),
    "",
    "PROPENSITY COMMON SUPPORT",
    "-------------------------",
    paste0(
        "Common support exists: ",
        common_support_exists
    ),
    paste0(
        "Common support interval: ",
        if (
            common_support_exists
        ) {
            sprintf(
                "[%.4f, %.4f]",
                common_lower,
                common_upper
            )
        } else {
            "NA"
        }
    ),
    capture.output(
        print(
            ps_quantiles
        )
    ),
    "",
    "OVERLAP METRICS",
    "---------------",
    capture.output(
        print(
            overlap_metrics
        )
    ),
    "",
    "PREDECLARED GATE CRITERIA",
    "-------------------------",
    capture.output(
        print(
            gate_table
        )
    ),
    "",
    paste0(
        "FINAL GATE DECISION: ",
        gate_decision
    ),
    gate_interpretation,
    "",
    "INTERPRETATION GUARDRAILS",
    "-------------------------",
    "- GO means only that an exploratory comparator model may be estimable.",
    "- GO does NOT imply causal exchangeability.",
    "- CONDITIONAL GO requires manual review before any interaction model.",
    "- NO-GO means no treatment x bypass interaction should be run.",
    "- If comparator modeling is not credible, this audit belongs in the main Discussion.",
    "- Do not relabel this as external validation or treatment-effect confirmation.",
    "",
    paste0(
        "Saved RDS: ",
        result_file
    ),
    paste0(
        "Saved Table S7: ",
        table_file
    )
)

writeLines(
    audit_lines,
    audit_file
)


# ============================================================
# 19. Word Table S7
# ============================================================

genomic_word <- copy(
    genomic_summary
)

genomic_word[
    ,
    `High-confidence+, %` :=
        sprintf(
            "%.1f",
            `High-confidence+, %`
        )
]

context_word <- copy(
    context_summary
)

setnames(
    context_word,
    "FIRST_BIOLOGIC_R",
    "Treatment"
)

overlap_word <- copy(
    overlap_metrics
)

if (
    nrow(
        overlap_word
    ) > 0L
) {

    numeric_fmt_cols <- intersect(
        c(
            "COMMON_SUPPORT_PCT",
            "BYPASS_POS_SUPPORT_PCT",
            "OVERLAP_ESS",
            "BYPASS_POS_OVERLAP_ESS"
        ),
        names(
            overlap_word
        )
    )

    for (v in numeric_fmt_cols) {

        overlap_word[
            ,
            (v) :=
                sprintf(
                    "%.1f",
                    get(v)
                )
        ]
    }
}

gate_word <- copy(
    gate_table
)

gate_word[
    ,
    Passed :=
        fifelse(
            Passed,
            "Yes",
            "No"
        )
]

ft_three_line <- function(
    df,
    font_size = 7.6
) {

    ft <- flextable(
        df
    )

    ft <- border_remove(
        ft
    )

    b_top <- fp_border(
        color = "black",
        width = 1.0
    )

    b_mid <- fp_border(
        color = "black",
        width = 0.6
    )

    ft <- hline_top(
        ft,
        border = b_top,
        part = "header"
    )

    ft <- hline_bottom(
        ft,
        border = b_mid,
        part = "header"
    )

    ft <- hline_bottom(
        ft,
        border = b_top,
        part = "body"
    )

    ft <- font(
        ft,
        fontname = "Arial",
        part = "all"
    )

    ft <- fontsize(
        ft,
        size = font_size,
        part = "all"
    )

    ft <- bold(
        ft,
        part = "header"
    )

    ft <- align(
        ft,
        align = "center",
        part = "all"
    )

    ft <- align(
        ft,
        j = 1,
        align = "left",
        part = "all"
    )

    ft <- valign(
        ft,
        valign = "center",
        part = "all"
    )

    ft <- autofit(
        ft
    )

    ft
}

doc <- read_docx()

sec <- prop_section(
    page_size =
        page_size(
            orient =
                "landscape"
        ),
    page_margins =
        page_mar(
            top = 0.45,
            bottom = 0.45,
            left = 0.45,
            right = 0.45
        )
)

doc <- body_set_default_section(
    doc,
    sec
)

doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Supplementary Table S7. Outcome-blind bevacizumab comparator feasibility audit",
            fp_text(
                font.family =
                    "Arial",
                font.size =
                    10.5,
                bold =
                    TRUE
            )
        )
    )
)

doc <- body_add_par(
    doc,
    "A. Genomic positivity by treatment arm"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        genomic_word,
        7.5
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "B. Outcome-blind treatment-context summary"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        context_word,
        7.2
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "C. Propensity-overlap diagnostics"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        overlap_word,
        7.4
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "D. Predeclared feasibility gate"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        gate_word,
        7.5
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            paste0(
                "Gate decision: ",
                gate_decision,
                ". ",
                gate_interpretation,
                " The propensity model and overlap diagnostics use no outcome information. ",
                "This feasibility audit does not estimate treatment effects or treatment-by-biomarker interaction."
            ),
            fp_text(
                font.family =
                    "Arial",
                font.size =
                    7.5
            )
        )
    )
)

print(
    doc,
    target =
        table_file
)


# ============================================================
# 20. Console summary
# ============================================================

cat("\nGENOMIC SUMMARY\n")
cat("---------------\n")
print(
    genomic_summary
)

cat("\nPRIOR-MAJOR 0/1/2 CELLS\n")
cat("-----------------------\n")
print(
    prior_cells
)

cat("\nPS COVERAGE / OVERLAP\n")
cat("---------------------\n")
print(
    ps_coverage
)
print(
    overlap_metrics
)

cat("\nGATE\n")
cat("----\n")
print(
    gate_table
)

cat(
    "\nFINAL GATE DECISION: ",
    gate_decision,
    "\n",
    sep = ""
)

cat(
    gate_interpretation,
    "\n"
)

cat("\n============================================================\n")
cat("B1-14 COMPLETE\n")
cat("============================================================\n")

cat(
    "\nAudit:\n",
    audit_file,
    "\n"
)

cat(
    "\nTable S7:\n",
    table_file,
    "\n"
)

cat(
    "\nRDS:\n",
    result_file,
    "\n"
)
