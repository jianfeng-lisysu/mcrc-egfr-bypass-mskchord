# ============================================================
# A7 / B-line
# 04_B1_04_first_biologic_salvage_audit.R
#
# PURPOSE
#   Audit whether the OLD 458-patient cohort can be reframed as
#   a "first relevant biologic exposure" active-comparator cohort
#   (first anti-EGFR vs first bevacizumab exposure), across lines.
#
# IMPORTANT
#   - No rwPFS / OS / treatment-effect model is run.
#   - Does NOT modify formal *_R.csv files.
#   - Reads the original MSK-CHORD TAR for raw treatment timelines.
#
# Main questions:
#   1) Are all four treatment x bypass cells adequately populated?
#   2) Was specimen acquisition / sequencing available before biologic?
#   3) What backbone was given around biologic initiation?
#   4) What prior cytotoxic exposure existed before biologic?
#   5) Does positivity remain within prior-treatment strata?
#
# Outputs:
#   06_logs_and_audit/B1_04_first_biologic_salvage_audit.txt
#   06_logs_and_audit/B1_04_first_biologic_patient_audit.csv
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
    stop("Missing main CSV:\n", main_file)
}

if (!file.exists(tar_file)) {
    stop("Missing MSK tar:\n", tar_file)
}

if (!dir.exists(audit_dir)) {
    dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}
library(data.table)

cat("\n============================================================\n")
cat("B1-04 FIRST BIOLOGIC SALVAGE AUDIT\n")
cat("============================================================\n\n")


# ============================================================
# Helpers
# ============================================================

min_or_na <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    min(x)
}

max_or_na <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    max(x)
}

pct <- function(n, d) {
    if (d == 0L) return("NA")
    sprintf("%d/%d (%.1f%%)", n, d, 100 * n / d)
}

median_iqr <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) return("NA")
    q <- quantile(
        x,
        probs = c(0.25, 0.50, 0.75),
        na.rm = TRUE,
        names = FALSE
    )
    sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

safe_max01 <- function(x) {
    if (length(x) == 0L) return(0L)
    as.integer(any(x %in% TRUE))
}

extract_site_code <- function(x) {

    x <- toupper(as.character(x))

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
# Formal main CSV: old 458 cohort
# ============================================================

main <- fread(
    main_file,
    check.names = FALSE,
    showProgress = TRUE
)

required_main <- c(
    "PATIENT_ID",
    "OLD_ELIGIBLE_803_R",
    "OLD_BIOLOGIC_458_R",
    "FIRST_BIOLOGIC_R",
    "BYPASS_ORIGINAL_R",
    "SPECIMEN_ACQUISITION_FIRST_DAY_R",
    "SEQUENCING_FIRST_DAY_R"
)

missing_main <- setdiff(
    required_main,
    names(main)
)

if (length(missing_main) > 0L) {
    stop(
        "Main CSV missing columns:\n",
        paste(
            missing_main,
            collapse = "\n"
        )
    )
}

cohort <- main[
    OLD_BIOLOGIC_458_R == 1L,
    .(
        PATIENT_ID,
        FIRST_BIOLOGIC_R_OLD = FIRST_BIOLOGIC_R,
        BYPASS_ORIGINAL_R,
        SPECIMEN_ACQUISITION_FIRST_DAY_R,
        SEQUENCING_FIRST_DAY_R
    )
]

if (nrow(cohort) != 458L) {
    stop(
        "Old biologic anchor failed: expected 458, got ",
        nrow(cohort)
    )
}

cat("Old first-relevant-biologic cohort:", nrow(cohort), "\n")


# ============================================================
# Temporary extraction
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_04_"
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

find_raw <- function(fn) {

    hit <- all_files[
        basename(all_files) == fn
    ]

    if (length(hit) == 0L) {
        stop("Missing raw MSK file: ", fn)
    }

    hit[1]
}


# ============================================================
# Raw diagnosis: explicit CRC Stage 4 anchor
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
    PATIENT_ID %chin% cohort$PATIENT_ID
]

diagnosis[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(START_DATE)
    )
]

diagnosis[
    ,
    SITE_CODE_R := extract_site_code(
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
    IS_EXPLICIT_STAGE4_CRC_R :=
        IS_CRC_DX_R &
        (
            STAGE_CDM_DERIVED == "Stage 4" |
            grepl(
                "Distant",
                SUMMARY,
                ignore.case = TRUE
            )
        )
]

stage4 <- diagnosis[
    IS_EXPLICIT_STAGE4_CRC_R == TRUE,
    .(
        EXPLICIT_STAGE4_CRC_DAY_R =
            min_or_na(
                START_DATE_R
            )
    ),
    by = PATIENT_ID
]


# ============================================================
# Raw treatment timeline
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
    PATIENT_ID %chin% cohort$PATIENT_ID
]

treatment[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(START_DATE)
    )
]

treatment[
    ,
    AGENT_UPPER_R := toupper(
        trimws(
            fifelse(
                is.na(AGENT),
                "",
                AGENT
            )
        )
    )
]

treatment[
    ,
    IS_5FU_R :=
        AGENT_UPPER_R %chin%
        c(
            "FLUOROURACIL",
            "5-FLUOROURACIL",
            "5-FU"
        )
]

treatment[
    ,
    IS_CAPE_R :=
        AGENT_UPPER_R ==
        "CAPECITABINE"
]

treatment[
    ,
    IS_FP_R :=
        IS_5FU_R |
        IS_CAPE_R
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

treatment[
    ,
    IS_RELEVANT_BIOLOGIC_R :=
        IS_EGFR_R |
        IS_BEV_R
]


# ============================================================
# Reconstruct first relevant biologic exactly from raw timeline
# ============================================================

first_bio <- treatment[
    IS_RELEVANT_BIOLOGIC_R == TRUE &
        is.finite(
            START_DATE_R
        ),
    {
        d <- min(
            START_DATE_R
        )

        same_day <- .SD[
            START_DATE_R == d
        ]

        egfr <- any(
            same_day$IS_EGFR_R
        )

        bev <- any(
            same_day$IS_BEV_R
        )

        group <- if (egfr && bev) {
            "dual"
        } else if (egfr) {
            "anti-EGFR"
        } else if (bev) {
            "bevacizumab"
        } else {
            NA_character_
        }

        .(
            FIRST_BIOLOGIC_DAY_R = d,
            FIRST_BIOLOGIC_R = group,
            FIRST_BIOLOGIC_AGENTS_R = paste(
                sort(
                    unique(
                        same_day$AGENT[
                            !is.na(
                                same_day$AGENT
                            )
                        ]
                    )
                ),
                collapse = "|"
            )
        )
    },
    by = PATIENT_ID
]

audit <- merge(
    cohort,
    first_bio,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

audit <- merge(
    audit,
    stage4,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

# Hard treatment-group audit
n_egfr <- audit[
    FIRST_BIOLOGIC_R ==
        "anti-EGFR",
    .N
]

n_bev <- audit[
    FIRST_BIOLOGIC_R ==
        "bevacizumab",
    .N
]

n_dual <- audit[
    FIRST_BIOLOGIC_R ==
        "dual",
    .N
]

if (
    n_egfr != 230L ||
    n_bev != 228L ||
    n_dual != 0L
) {

    stop(
        paste0(
            "First-biologic anchor failed. ",
            "Expected anti-EGFR 230 / BEV 228 / dual 0; got ",
            n_egfr,
            " / ",
            n_bev,
            " / ",
            n_dual
        )
    )
}


# ============================================================
# Around-biologic backbone: +/- 30 days
# ============================================================

build_bio_context <- function(pid, bio_day) {

    ptx <- treatment[
        PATIENT_ID == pid &
            is.finite(
                START_DATE_R
            )
    ]

    if (!is.finite(bio_day) ||
        nrow(ptx) == 0L) {

        return(
            data.table(
                PATIENT_ID = pid,
                BACKBONE_PM30_R = NA_character_,
                BACKBONE_DOUBLET_PM30_R = 0L,
                TRIPLET_PM30_R = 0L,
                AGENTS_PM30_R = NA_character_,
                PRIOR_FP_R = 0L,
                PRIOR_OX_R = 0L,
                PRIOR_IRI_R = 0L,
                PRIOR_OX_AND_IRI_R = 0L,
                PRIOR_CYTOTOXIC_PATTERN_R = NA_character_,
                PRIOR_RELEVANT_TX_EVENT_N_R = 0L
            )
        )
    }

    same_line <- ptx[
        START_DATE_R >= bio_day - 30 &
            START_DATE_R <= bio_day + 30
    ]

    prior <- ptx[
        START_DATE_R <
            bio_day - 30
    ]

    has_5fu <- any(
        same_line$IS_5FU_R
    )

    has_cape <- any(
        same_line$IS_CAPE_R
    )

    has_fp <- (
        has_5fu |
        has_cape
    )

    has_ox <- any(
        same_line$IS_OX_R
    )

    has_iri <- any(
        same_line$IS_IRI_R
    )

    triplet <- (
        has_fp &
        has_ox &
        has_iri
    )

    doublet <- (
        has_fp &
        xor(
            has_ox,
            has_iri
        )
    )

    backbone <- fcase(
        triplet,
        "triplet/FOLFOXIRI-like",

        has_5fu & has_ox & !has_iri,
        "FOLFOX-like",

        has_cape & has_ox & !has_iri,
        "CAPOX-like",

        has_5fu & has_iri & !has_ox,
        "FOLFIRI-like",

        has_cape & has_iri & !has_ox,
        "CAPIRI-like",

        has_fp & !has_ox & !has_iri,
        "fluoropyrimidine-only",

        !has_fp & has_ox & !has_iri,
        "oxaliplatin-only",

        !has_fp & !has_ox & has_iri,
        "irinotecan-only",

        default =
            "other/no cytotoxic backbone"
    )

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

    data.table(
        PATIENT_ID = pid,

        BACKBONE_PM30_R =
            backbone,

        BACKBONE_DOUBLET_PM30_R =
            as.integer(
                doublet
            ),

        TRIPLET_PM30_R =
            as.integer(
                triplet
            ),

        AGENTS_PM30_R =
            paste(
                sort(
                    unique(
                        same_line$AGENT[
                            !is.na(
                                same_line$AGENT
                            )
                        ]
                    )
                ),
                collapse = "|"
            ),

        PRIOR_FP_R =
            as.integer(
                prior_fp
            ),

        PRIOR_OX_R =
            as.integer(
                prior_ox
            ),

        PRIOR_IRI_R =
            as.integer(
                prior_iri
            ),

        PRIOR_OX_AND_IRI_R =
            as.integer(
                prior_ox &
                prior_iri
            ),

        PRIOR_CYTOTOXIC_PATTERN_R =
            prior_pattern,

        PRIOR_RELEVANT_TX_EVENT_N_R =
            prior[
                IS_FP_R |
                    IS_OX_R |
                    IS_IRI_R,
                .N
            ]
    )
}

bio_context <- rbindlist(
    Map(
        build_bio_context,
        audit$PATIENT_ID,
        audit$FIRST_BIOLOGIC_DAY_R
    ),
    fill = TRUE
)

audit <- merge(
    audit,
    bio_context,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)


# ============================================================
# Baseline timing
# ============================================================

audit[
    ,
    SPECIMEN_PRE_BIOLOGIC_R :=
        as.integer(
            is.finite(
                SPECIMEN_ACQUISITION_FIRST_DAY_R
            ) &
            is.finite(
                FIRST_BIOLOGIC_DAY_R
            ) &
            SPECIMEN_ACQUISITION_FIRST_DAY_R <=
                FIRST_BIOLOGIC_DAY_R
        )
]

audit[
    ,
    SEQUENCING_PRE_BIOLOGIC_R :=
        as.integer(
            is.finite(
                SEQUENCING_FIRST_DAY_R
            ) &
            is.finite(
                FIRST_BIOLOGIC_DAY_R
            ) &
            SEQUENCING_FIRST_DAY_R <=
                FIRST_BIOLOGIC_DAY_R
        )
]

audit[
    ,
    STAGE4_TO_BIOLOGIC_DAYS_R :=
        FIRST_BIOLOGIC_DAY_R -
        EXPLICIT_STAGE4_CRC_DAY_R
]

audit[
    ,
    STAGE4_TO_BIOLOGIC_BAND_R := cut(
        STAGE4_TO_BIOLOGIC_DAYS_R,
        breaks = c(
            -Inf,
            30,
            60,
            180,
            365,
            Inf
        ),
        labels = c(
            "<=30 d",
            "31-60 d",
            "61-180 d",
            "181-365 d",
            ">365 d"
        ),
        right = TRUE
    )
]


# ============================================================
# Core cell tables
# ============================================================

core_cells <- audit[
    ,
    .N,
    by = .(
        FIRST_BIOLOGIC_R,
        BYPASS_ORIGINAL_R
    )
][order(
    FIRST_BIOLOGIC_R,
    BYPASS_ORIGINAL_R
)]

specimen_cells <- audit[
    SPECIMEN_PRE_BIOLOGIC_R == 1L,
    .N,
    by = .(
        FIRST_BIOLOGIC_R,
        BYPASS_ORIGINAL_R
    )
][order(
    FIRST_BIOLOGIC_R,
    BYPASS_ORIGINAL_R
)]

sequencing_cells <- audit[
    SEQUENCING_PRE_BIOLOGIC_R == 1L,
    .N,
    by = .(
        FIRST_BIOLOGIC_R,
        BYPASS_ORIGINAL_R
    )
][order(
    FIRST_BIOLOGIC_R,
    BYPASS_ORIGINAL_R
)]

doublet_cells <- audit[
    BACKBONE_DOUBLET_PM30_R == 1L,
    .N,
    by = .(
        FIRST_BIOLOGIC_R,
        BYPASS_ORIGINAL_R
    )
][order(
    FIRST_BIOLOGIC_R,
    BYPASS_ORIGINAL_R
)]

backbone_tab <- audit[
    ,
    .N,
    by = .(
        FIRST_BIOLOGIC_R,
        BACKBONE_PM30_R
    )
][order(
    FIRST_BIOLOGIC_R,
    -N
)]

prior_tab <- audit[
    ,
    .N,
    by = .(
        FIRST_BIOLOGIC_R,
        PRIOR_CYTOTOXIC_PATTERN_R
    )
][order(
    FIRST_BIOLOGIC_R,
    -N
)]

prior_cells <- audit[
    ,
    .N,
    by = .(
        PRIOR_CYTOTOXIC_PATTERN_R,
        FIRST_BIOLOGIC_R,
        BYPASS_ORIGINAL_R
    )
][order(
    PRIOR_CYTOTOXIC_PATTERN_R,
    FIRST_BIOLOGIC_R,
    BYPASS_ORIGINAL_R
)]

timing_tab <- audit[
    !is.na(
        STAGE4_TO_BIOLOGIC_BAND_R
    ),
    .N,
    by = .(
        FIRST_BIOLOGIC_R,
        STAGE4_TO_BIOLOGIC_BAND_R
    )
][order(
    FIRST_BIOLOGIC_R,
    STAGE4_TO_BIOLOGIC_BAND_R
)]


# ============================================================
# Quantities for report
# ============================================================

n458 <- nrow(audit)

n_spec_pre <- audit[
    SPECIMEN_PRE_BIOLOGIC_R == 1L,
    .N
]

n_seq_pre <- audit[
    SEQUENCING_PRE_BIOLOGIC_R == 1L,
    .N
]

n_doublet <- audit[
    BACKBONE_DOUBLET_PM30_R == 1L,
    .N
]

n_triplet <- audit[
    TRIPLET_PM30_R == 1L,
    .N
]

n_stage4_anchor <- audit[
    is.finite(
        EXPLICIT_STAGE4_CRC_DAY_R
    ),
    .N
]

# positivity check helper
cell_min <- function(tab) {
    if (nrow(tab) == 0L) return(0L)
    min(tab$N)
}

core_min <- cell_min(
    core_cells
)

specimen_min <- cell_min(
    specimen_cells
)

sequencing_min <- cell_min(
    sequencing_cells
)

doublet_min <- cell_min(
    doublet_cells
)


# ============================================================
# Patient-level audit output
# ============================================================

patient_audit_file <- file.path(
    audit_dir,
    "B1_04_first_biologic_patient_audit.csv"
)

fwrite(
    audit,
    patient_audit_file,
    na = ""
)


# ============================================================
# Main text audit
# ============================================================

audit_file <- file.path(
    audit_dir,
    "B1_04_first_biologic_salvage_audit.txt"
)

lines <- c(
    "B1-04 FIRST RELEVANT BIOLOGIC SALVAGE AUDIT",
    "===========================================",
    "",
    "1. Cohort anchors",
    "-----------------",
    paste0(
        "Old first-relevant-biologic cohort: ",
        n458
    ),
    paste0(
        "anti-EGFR: ",
        n_egfr
    ),
    paste0(
        "bevacizumab: ",
        n_bev
    ),
    paste0(
        "dual on first day: ",
        n_dual
    ),
    "",
    "Treatment x bypass cells:",
    capture.output(
        print(
            core_cells
        )
    ),
    paste0(
        "Minimum overall cell size: ",
        core_min
    ),
    "",
    "2. Baseline genomic timing",
    "--------------------------",
    paste0(
        "Specimen acquisition <= first biologic: ",
        pct(
            n_spec_pre,
            n458
        )
    ),
    paste0(
        "Sequencing <= first biologic: ",
        pct(
            n_seq_pre,
            n458
        )
    ),
    "",
    "Treatment x bypass cells among specimen-pre-biologic patients:",
    capture.output(
        print(
            specimen_cells
        )
    ),
    paste0(
        "Minimum specimen-pre-biologic cell size: ",
        specimen_min
    ),
    "",
    "Treatment x bypass cells among sequencing-pre-biologic patients:",
    capture.output(
        print(
            sequencing_cells
        )
    ),
    paste0(
        "Minimum sequencing-pre-biologic cell size: ",
        sequencing_min
    ),
    "",
    "3. Backbone around biologic initiation (+/-30 days)",
    "---------------------------------------------------",
    paste0(
        "Doublet backbone: ",
        pct(
            n_doublet,
            n458
        )
    ),
    paste0(
        "Triplet/FOLFOXIRI-like: ",
        pct(
            n_triplet,
            n458
        )
    ),
    "",
    capture.output(
        print(
            backbone_tab
        )
    ),
    "",
    "Treatment x bypass cells restricted to doublet backbone:",
    capture.output(
        print(
            doublet_cells
        )
    ),
    paste0(
        "Minimum doublet-backbone cell size: ",
        doublet_min
    ),
    "",
    "4. Prior cytotoxic exposure before first biologic",
    "-------------------------------------------------",
    capture.output(
        print(
            prior_tab
        )
    ),
    "",
    "Treatment x bypass cells within prior-treatment strata:",
    capture.output(
        print(
            prior_cells
        )
    ),
    "",
    "5. Timing from explicit Stage 4 diagnosis",
    "----------------------------------------",
    paste0(
        "Explicit CRC Stage 4 anchor available: ",
        pct(
            n_stage4_anchor,
            n458
        )
    ),
    paste0(
        "Stage 4 to first biologic median [Q1,Q3]: ",
        median_iqr(
            audit$
                STAGE4_TO_BIOLOGIC_DAYS_R
        ),
        " days"
    ),
    "",
    capture.output(
        print(
            timing_tab
        )
    ),
    "",
    "6. Audit interpretation",
    "-----------------------",
    "This is an AUDIT only; no rwPFS / OS / treatment-effect estimate was calculated.",
    "The salvage design is viable only if treatment x bypass positivity remains acceptable after baseline-timing and prior-treatment restrictions.",
    "If cells collapse within prior-treatment strata, do not proceed with an active-comparator interaction model."
)

writeLines(
    lines,
    audit_file
)

cat("\n============================================================\n")
cat("B1-04 AUDIT COMPLETE\n")
cat("============================================================\n")
cat(
    "anti-EGFR / BEV:",
    n_egfr,
    "/",
    n_bev,
    "\n"
)
cat(
    "Specimen pre-biologic:",
    n_spec_pre,
    "/",
    n458,
    "\n"
)
cat(
    "Sequencing pre-biologic:",
    n_seq_pre,
    "/",
    n458,
    "\n"
)
cat(
    "Doublet backbone:",
    n_doublet,
    "/",
    n458,
    "\n"
)
cat("\nMain audit:\n", audit_file, "\n")
cat("\nPatient audit:\n", patient_audit_file, "\n")
