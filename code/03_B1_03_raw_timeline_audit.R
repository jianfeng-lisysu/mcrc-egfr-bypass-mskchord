# ============================================================
# A7 / B-line
# 03_B1_03_raw_timeline_audit.R
#
# Corrected B1-03 audit.
#
# IMPORTANT:
# - Does NOT parse collapsed DIAGNOSIS_SEQUENCE_R.
# - Reads the ORIGINAL MSK-CHORD timeline files from the TAR.
# - Uses explicit CRC-specific Stage 4 registry diagnosis as the
#   primary metastatic anchor.
# - Characterizes fallback evidence for patients without an
#   explicit CRC Stage 4 diagnosis, but does NOT silently include
#   them in the primary first-line cohort.
# - Does NOT run rwPFS / OS / treatment-effect models.
# - Does NOT modify either formal *_R.csv file.
#
# Outputs:
#   06_logs_and_audit/B1_03_raw_timeline_audit.txt
#   06_logs_and_audit/B1_03_raw_timeline_patient_audit.csv
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
cat("B1-03 CORRECTED RAW-TIMELINE AUDIT\n")
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

row_min_na <- function(a, b) {
    out <- pmin(
        a,
        b,
        na.rm = TRUE
    )
    out[is.infinite(out)] <- NA_real_
    out
}


# ============================================================
# Read formal main CSV and retain old eligible 803
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
    OLD_ELIGIBLE_803_R == 1L,
    .(
        PATIENT_ID,
        OLD_BIOLOGIC_458_R,
        BYPASS_ORIGINAL_R,
        SPECIMEN_ACQUISITION_FIRST_DAY_R,
        SEQUENCING_FIRST_DAY_R
    )
]

if (nrow(cohort) != 803L) {
    stop(
        "Old eligible anchor failed: expected 803, got ",
        nrow(cohort)
    )
}

cat("Old eligible cohort:", nrow(cohort), "\n")


# ============================================================
# Temporary extraction of only required raw files
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_03_raw_"
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
# Raw clinical sample table:
# detect patients with an additional non-CRC sequenced cancer.
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

sample_803 <- sample[
    PATIENT_ID %chin% cohort$PATIENT_ID
]

other_cancer <- sample_803[
    CANCER_TYPE != "Colorectal Cancer",
    .(
        OTHER_SEQUENCED_CANCER_N_R = .N,
        OTHER_SEQUENCED_CANCER_R = paste(
            unique(
                paste(
                    CANCER_TYPE,
                    CANCER_TYPE_DETAILED,
                    sep = ":"
                )
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]


# ============================================================
# RAW diagnosis timeline
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

# Include standard colon/rectum codes plus C218 ("Rectum, other parts")
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

crc_dx_sum <- diagnosis[
    IS_CRC_DX_R == TRUE,
    .(
        CRC_DX_FIRST_DAY_R =
            min_or_na(START_DATE_R),

        CRC_DX_EVENT_N_R = .N,

        CRC_DX_STAGE_SEQUENCE_R = paste(
            paste(
                START_DATE_R,
                SITE_CODE_R,
                STAGE_CDM_DERIVED,
                trimws(SUMMARY),
                sep = "~"
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]

explicit_stage4 <- diagnosis[
    IS_EXPLICIT_STAGE4_CRC_R == TRUE,
    .(
        EXPLICIT_STAGE4_CRC_DAY_R =
            min_or_na(START_DATE_R),

        EXPLICIT_STAGE4_CRC_EVENT_N_R = .N
    ),
    by = PATIENT_ID
]

non_crc_dx <- diagnosis[
    IS_CRC_DX_R == FALSE,
    .(
        NONCRC_DX_N_R = .N,
        NONCRC_DX_R = paste(
            unique(
                paste(
                    START_DATE_R,
                    DX_DESCRIPTION,
                    STAGE_CDM_DERIVED,
                    sep = "~"
                )
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]


# ============================================================
# RAW tumor-sites timeline
# Characterize fallback evidence only.
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
    PATIENT_ID %chin% cohort$PATIENT_ID
]

tumor_sites[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(START_DATE)
    )
]

tumor_site_sum <- tumor_sites[
    ,
    .(
        TUMOR_SITE_FIRST_DAY_ANY_R =
            min_or_na(START_DATE_R),

        TUMOR_SITE_EVENT_N_R = .N,

        TUMOR_SITE_TYPES_R = paste(
            sort(
                unique(
                    TUMOR_SITE[
                        !is.na(TUMOR_SITE)
                    ]
                )
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]


# ============================================================
# RAW progression timeline
# Characterize fallback evidence only.
# ============================================================

progression <- fread(
    find_raw(
        "data_timeline_progression.txt"
    ),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

progression <- progression[
    PATIENT_ID %chin% cohort$PATIENT_ID
]

progression[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(START_DATE)
    )
]

progression_y <- progression[
    PROGRESSION == "Y",
    .(
        PROGRESSION_Y_FIRST_DAY_R =
            min_or_na(START_DATE_R),

        PROGRESSION_Y_EVENT_N_R = .N
    ),
    by = PATIENT_ID
]


# ============================================================
# RAW treatment timeline
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
    IS_5FU_R := AGENT_UPPER_R %chin%
        c(
            "FLUOROURACIL",
            "5-FLUOROURACIL",
            "5-FU"
        )
]

treatment[
    ,
    IS_CAPE_R :=
        AGENT_UPPER_R == "CAPECITABINE"
]

treatment[
    ,
    IS_OX_R :=
        AGENT_UPPER_R == "OXALIPLATIN"
]

treatment[
    ,
    IS_IRI_R :=
        AGENT_UPPER_R == "IRINOTECAN"
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
        AGENT_UPPER_R == "BEVACIZUMAB"
]

treatment[
    ,
    IS_RELEVANT_CRC_AGENT_R :=
        IS_5FU_R |
        IS_CAPE_R |
        IS_OX_R |
        IS_IRI_R |
        IS_EGFR_R |
        IS_BEV_R
]


# ============================================================
# Assemble patient-level anchor table
# ============================================================

audit <- copy(
    cohort
)

for (obj in list(
    crc_dx_sum,
    explicit_stage4,
    other_cancer,
    non_crc_dx,
    tumor_site_sum,
    progression_y
)) {

    audit <- merge(
        audit,
        obj,
        by = "PATIENT_ID",
        all.x = TRUE,
        sort = FALSE
    )
}

audit[
    ,
    OTHER_CANCER_FLAG_R := as.integer(
        (
            !is.na(OTHER_SEQUENCED_CANCER_N_R) &
            OTHER_SEQUENCED_CANCER_N_R > 0
        ) |
        (
            !is.na(NONCRC_DX_N_R) &
            NONCRC_DX_N_R > 0
        )
    )
]

audit[
    ,
    EXPLICIT_STAGE4_CRC_R := as.integer(
        is.finite(
            EXPLICIT_STAGE4_CRC_DAY_R
        )
    )
]

# For descriptive fallback only:
# require tumor-site/progression evidence not to be >30 days before
# the CRC diagnosis, to reduce obvious pre-CRC contamination.
audit[
    ,
    FALLBACK_TUMOR_SITE_DAY_R := fifelse(
        is.finite(TUMOR_SITE_FIRST_DAY_ANY_R) &
            is.finite(CRC_DX_FIRST_DAY_R) &
            TUMOR_SITE_FIRST_DAY_ANY_R >=
                CRC_DX_FIRST_DAY_R - 30,
        TUMOR_SITE_FIRST_DAY_ANY_R,
        NA_real_
    )
]

audit[
    ,
    FALLBACK_PROGRESSION_Y_DAY_R := fifelse(
        is.finite(PROGRESSION_Y_FIRST_DAY_R) &
            is.finite(CRC_DX_FIRST_DAY_R) &
            PROGRESSION_Y_FIRST_DAY_R >=
                CRC_DX_FIRST_DAY_R - 30,
        PROGRESSION_Y_FIRST_DAY_R,
        NA_real_
    )
]

audit[
    ,
    FALLBACK_MET_EVIDENCE_DAY_R := row_min_na(
        FALLBACK_TUMOR_SITE_DAY_R,
        FALLBACK_PROGRESSION_Y_DAY_R
    )
]


# ============================================================
# First-line reconstruction using EXPLICIT CRC Stage 4 only
#
# Candidate line start:
# earliest relevant CRC systemic agent at/after
# explicit Stage 4 date - 30 days.
#
# Chemotherapy backbone is assembled within 30 days.
# Biologic is audited within both 30 and 60 days.
# ============================================================

build_line <- function(pid, anchor_day) {

    ptx <- treatment[
        PATIENT_ID == pid &
            IS_RELEVANT_CRC_AGENT_R == TRUE &
            is.finite(START_DATE_R)
    ]

    if (!is.finite(anchor_day) ||
        nrow(ptx) == 0L) {

        return(
            data.table(
                PATIENT_ID = pid,
                LINE_START_DAY_R = NA_real_,
                BACKBONE_30D_R = NA_character_,
                BIOLOGIC_30D_R = NA_character_,
                BIOLOGIC_60D_R = NA_character_,
                BIOLOGIC_FIRST_DAY_R = NA_real_,
                BIOLOGIC_DELAY_DAYS_R = NA_real_,
                TRIPLET_30D_R = NA_integer_,
                DOUBLET_30D_R = NA_integer_,
                QUALIFY_30D_R = 0L,
                QUALIFY_60D_R = 0L,
                LINE_AGENTS_30D_R = NA_character_,
                LINE_AGENTS_60D_R = NA_character_
            )
        )
    }

    setorder(
        ptx,
        START_DATE_R
    )

    candidate <- ptx[
        START_DATE_R >= anchor_day - 30
    ]

    if (nrow(candidate) == 0L) {

        return(
            data.table(
                PATIENT_ID = pid,
                LINE_START_DAY_R = NA_real_,
                BACKBONE_30D_R = NA_character_,
                BIOLOGIC_30D_R = NA_character_,
                BIOLOGIC_60D_R = NA_character_,
                BIOLOGIC_FIRST_DAY_R = NA_real_,
                BIOLOGIC_DELAY_DAYS_R = NA_real_,
                TRIPLET_30D_R = NA_integer_,
                DOUBLET_30D_R = NA_integer_,
                QUALIFY_30D_R = 0L,
                QUALIFY_60D_R = 0L,
                LINE_AGENTS_30D_R = NA_character_,
                LINE_AGENTS_60D_R = NA_character_
            )
        )
    }

    line_start <- min(
        candidate$START_DATE_R
    )

    w30 <- ptx[
        START_DATE_R >= line_start &
            START_DATE_R <= line_start + 30
    ]

    w60 <- ptx[
        START_DATE_R >= line_start &
            START_DATE_R <= line_start + 60
    ]

    has_5fu <- any(w30$IS_5FU_R)
    has_cape <- any(w30$IS_CAPE_R)
    has_fp <- has_5fu | has_cape
    has_ox <- any(w30$IS_OX_R)
    has_iri <- any(w30$IS_IRI_R)

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

        default =
            "other/non-doublet"
    )

    egfr30 <- any(w30$IS_EGFR_R)
    bev30 <- any(w30$IS_BEV_R)

    egfr60 <- any(w60$IS_EGFR_R)
    bev60 <- any(w60$IS_BEV_R)

    bio30 <- fcase(
        egfr30 & bev30,
        "dual",

        egfr30,
        "anti-EGFR",

        bev30,
        "bevacizumab",

        default =
            "none"
    )

    bio60 <- fcase(
        egfr60 & bev60,
        "dual",

        egfr60,
        "anti-EGFR",

        bev60,
        "bevacizumab",

        default =
            "none"
    )

    bio_rows <- ptx[
        (IS_EGFR_R | IS_BEV_R) &
            START_DATE_R >= line_start &
            START_DATE_R <= line_start + 60
    ]

    bio_day <- if (nrow(bio_rows) > 0L) {
        min(
            bio_rows$START_DATE_R
        )
    } else {
        NA_real_
    }

    bio_delay <- if (
        is.finite(bio_day)
    ) {
        bio_day - line_start
    } else {
        NA_real_
    }

    data.table(
        PATIENT_ID = pid,

        LINE_START_DAY_R =
            line_start,

        BACKBONE_30D_R =
            backbone,

        BIOLOGIC_30D_R =
            bio30,

        BIOLOGIC_60D_R =
            bio60,

        BIOLOGIC_FIRST_DAY_R =
            bio_day,

        BIOLOGIC_DELAY_DAYS_R =
            bio_delay,

        TRIPLET_30D_R =
            as.integer(
                triplet
            ),

        DOUBLET_30D_R =
            as.integer(
                doublet
            ),

        QUALIFY_30D_R =
            as.integer(
                doublet &
                bio30 %chin%
                    c(
                        "anti-EGFR",
                        "bevacizumab"
                    )
            ),

        QUALIFY_60D_R =
            as.integer(
                doublet &
                bio60 %chin%
                    c(
                        "anti-EGFR",
                        "bevacizumab"
                    )
            ),

        LINE_AGENTS_30D_R =
            paste(
                sort(
                    unique(
                        w30$AGENT[
                            !is.na(w30$AGENT)
                        ]
                    )
                ),
                collapse = "|"
            ),

        LINE_AGENTS_60D_R =
            paste(
                sort(
                    unique(
                        w60$AGENT[
                            !is.na(w60$AGENT)
                        ]
                    )
                ),
                collapse = "|"
            )
    )
}

line_audit <- rbindlist(
    Map(
        build_line,
        audit$PATIENT_ID,
        audit$EXPLICIT_STAGE4_CRC_DAY_R
    ),
    fill = TRUE
)

audit <- merge(
    audit,
    line_audit,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)


# ============================================================
# Baseline timing
# ============================================================

audit[
    ,
    SPECIMEN_PRE_LINE_R := fifelse(
        is.finite(
            SPECIMEN_ACQUISITION_FIRST_DAY_R
        ) &
        is.finite(
            LINE_START_DAY_R
        ),
        as.integer(
            SPECIMEN_ACQUISITION_FIRST_DAY_R <=
                LINE_START_DAY_R
        ),
        NA_integer_
    )
]

audit[
    ,
    SEQUENCING_PRE_LINE_R := fifelse(
        is.finite(
            SEQUENCING_FIRST_DAY_R
        ) &
        is.finite(
            LINE_START_DAY_R
        ),
        as.integer(
            SEQUENCING_FIRST_DAY_R <=
                LINE_START_DAY_R
        ),
        NA_integer_
    )
]

audit[
    ,
    SPECIMEN_PRE_BIOLOGIC_R := fifelse(
        is.finite(
            SPECIMEN_ACQUISITION_FIRST_DAY_R
        ) &
        is.finite(
            BIOLOGIC_FIRST_DAY_R
        ),
        as.integer(
            SPECIMEN_ACQUISITION_FIRST_DAY_R <=
                BIOLOGIC_FIRST_DAY_R
        ),
        NA_integer_
    )
]

audit[
    ,
    SEQUENCING_PRE_BIOLOGIC_R := fifelse(
        is.finite(
            SEQUENCING_FIRST_DAY_R
        ) &
        is.finite(
            BIOLOGIC_FIRST_DAY_R
        ),
        as.integer(
            SEQUENCING_FIRST_DAY_R <=
                BIOLOGIC_FIRST_DAY_R
        ),
        NA_integer_
    )
]


# ============================================================
# Old 458: first relevant biologic timing relative to explicit
# CRC Stage 4 anchor
# ============================================================

first_bio_all <- treatment[
    IS_EGFR_R | IS_BEV_R,
    .(
        OLD_FIRST_RELEVANT_BIOLOGIC_DAY_R =
            min_or_na(
                START_DATE_R
            )
    ),
    by = PATIENT_ID
]

audit <- merge(
    audit,
    first_bio_all,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

audit[
    ,
    OLD_BIOLOGIC_MINUS_STAGE4_DAYS_R :=
        OLD_FIRST_RELEVANT_BIOLOGIC_DAY_R -
        EXPLICIT_STAGE4_CRC_DAY_R
]


# ============================================================
# Summary counts
# ============================================================

n803 <- nrow(audit)

n_explicit <- audit[
    EXPLICIT_STAGE4_CRC_R == 1L,
    .N
]

n_no_explicit <- audit[
    EXPLICIT_STAGE4_CRC_R == 0L,
    .N
]

n_other_cancer <- audit[
    OTHER_CANCER_FLAG_R == 1L,
    .N
]

n_no_explicit_fallback <- audit[
    EXPLICIT_STAGE4_CRC_R == 0L &
        is.finite(
            FALLBACK_MET_EVIDENCE_DAY_R
        ),
    .N
]

n_line <- audit[
    EXPLICIT_STAGE4_CRC_R == 1L &
        is.finite(
            LINE_START_DAY_R
        ),
    .N
]

backbone_tab <- audit[
    EXPLICIT_STAGE4_CRC_R == 1L,
    .N,
    by = BACKBONE_30D_R
][order(-N)]

bio30_tab <- audit[
    EXPLICIT_STAGE4_CRC_R == 1L,
    .N,
    by = BIOLOGIC_30D_R
][order(-N)]

bio60_tab <- audit[
    EXPLICIT_STAGE4_CRC_R == 1L,
    .N,
    by = BIOLOGIC_60D_R
][order(-N)]

q30 <- audit[
    EXPLICIT_STAGE4_CRC_R == 1L &
        QUALIFY_30D_R == 1L
]

q60 <- audit[
    EXPLICIT_STAGE4_CRC_R == 1L &
        QUALIFY_60D_R == 1L
]

q30_tab <- q30[
    ,
    .N,
    by = .(
        BIOLOGIC_30D_R,
        BYPASS_ORIGINAL_R
    )
][order(
    BIOLOGIC_30D_R,
    BYPASS_ORIGINAL_R
)]

q60_tab <- q60[
    ,
    .N,
    by = .(
        BIOLOGIC_60D_R,
        BYPASS_ORIGINAL_R
    )
][order(
    BIOLOGIC_60D_R,
    BYPASS_ORIGINAL_R
)]

old458_explicit <- audit[
    OLD_BIOLOGIC_458_R == 1L &
        EXPLICIT_STAGE4_CRC_R == 1L &
        is.finite(
            OLD_BIOLOGIC_MINUS_STAGE4_DAYS_R
        )
]

old458_delay_bands <- old458_explicit[
    ,
    .(
        N = .N
    ),
    by = .(
        BAND = cut(
            OLD_BIOLOGIC_MINUS_STAGE4_DAYS_R,
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
    )
][order(BAND)]

old458_overlap30 <- audit[
    OLD_BIOLOGIC_458_R == 1L &
        EXPLICIT_STAGE4_CRC_R == 1L &
        QUALIFY_30D_R == 1L,
    .N
]

old458_overlap60 <- audit[
    OLD_BIOLOGIC_458_R == 1L &
        EXPLICIT_STAGE4_CRC_R == 1L &
        QUALIFY_60D_R == 1L,
    .N
]


# ============================================================
# Baseline timing helper among qualifying cohort
# ============================================================

baseline_line <- function(data, var, label) {

    den <- data[
        !is.na(
            get(var)
        ),
        .N
    ]

    num <- data[
        get(var) == 1L,
        .N
    ]

    paste0(
        label,
        ": ",
        pct(
            num,
            den
        )
    )
}


# ============================================================
# Write patient-level audit
# ============================================================

patient_audit_file <- file.path(
    audit_dir,
    "B1_03_raw_timeline_patient_audit.csv"
)

fwrite(
    audit,
    patient_audit_file,
    na = ""
)


# ============================================================
# Write main text audit
# ============================================================

audit_file <- file.path(
    audit_dir,
    "B1_03_raw_timeline_audit.txt"
)

lines <- c(
    "B1-03 CORRECTED RAW-TIMELINE AUDIT",
    "==================================",
    "",
    "1. Old feasibility cohort",
    "--------------------------",
    paste0(
        "Old eligible cohort: ",
        n803
    ),
    paste0(
        "Explicit CRC-specific Stage 4 registry diagnosis: ",
        pct(
            n_explicit,
            n803
        )
    ),
    paste0(
        "No explicit CRC-specific Stage 4 registry diagnosis: ",
        n_no_explicit
    ),
    paste0(
        "Patients with another sequenced / registry cancer: ",
        n_other_cancer
    ),
    paste0(
        "No-explicit-Stage4 patients with fallback tumor-site/progression evidence: ",
        n_no_explicit_fallback,
        "/",
        n_no_explicit
    ),
    "",
    "Primary B1-03 line reconstruction below uses EXPLICIT CRC-specific Stage 4 only.",
    "Fallback metastatic evidence is descriptive only and is not silently included.",
    "",
    "2. Candidate first metastatic line",
    "----------------------------------",
    paste0(
        "Explicit Stage 4 patients with a relevant CRC systemic line after anchor: ",
        pct(
            n_line,
            n_explicit
        )
    ),
    "",
    "Backbone within 30 days:",
    capture.output(
        print(
            backbone_tab
        )
    ),
    "",
    "Biologic within 30 days:",
    capture.output(
        print(
            bio30_tab
        )
    ),
    "",
    "Biologic within 60 days:",
    capture.output(
        print(
            bio60_tab
        )
    ),
    "",
    "3. Qualifying doublet + single biologic",
    "--------------------------------------",
    paste0(
        "30-day qualifying total: ",
        nrow(q30)
    ),
    paste0(
        "60-day qualifying total: ",
        nrow(q60)
    ),
    "",
    "30-day treatment x bypass cells:",
    capture.output(
        print(
            q30_tab
        )
    ),
    "",
    "60-day treatment x bypass cells:",
    capture.output(
        print(
            q60_tab
        )
    ),
    "",
    "4. Relation to old 458",
    "----------------------",
    paste0(
        "Old 458 overlapping corrected 30-day first-line cohort: ",
        old458_overlap30,
        "/458"
    ),
    paste0(
        "Old 458 overlapping corrected 60-day first-line cohort: ",
        old458_overlap60,
        "/458"
    ),
    paste0(
        "Among old 458 with explicit Stage 4 anchor, first relevant biologic delay median [Q1,Q3]: ",
        median_iqr(
            old458_explicit$
                OLD_BIOLOGIC_MINUS_STAGE4_DAYS_R
        ),
        " days"
    ),
    "",
    "Old 458 first relevant biologic timing after Stage 4 anchor:",
    capture.output(
        print(
            old458_delay_bands
        )
    ),
    "",
    "5. Baseline timing among corrected 30-day qualifying cohort",
    "-----------------------------------------------------------",
    baseline_line(
        q30,
        "SPECIMEN_PRE_LINE_R",
        "Specimen acquisition <= line start"
    ),
    baseline_line(
        q30,
        "SEQUENCING_PRE_LINE_R",
        "Sequencing <= line start"
    ),
    baseline_line(
        q30,
        "SPECIMEN_PRE_BIOLOGIC_R",
        "Specimen acquisition <= biologic start"
    ),
    baseline_line(
        q30,
        "SEQUENCING_PRE_BIOLOGIC_R",
        "Sequencing <= biologic start"
    ),
    "",
    "6. Interpretation",
    "-----------------",
    "This is still an AUDIT, not an outcome analysis.",
    "Do not run rwPFS/OS models until the treatment x bypass cell sizes and baseline timing are reviewed.",
    "No treatment-effect estimate was calculated."
)

writeLines(
    lines,
    audit_file
)

cat("\n============================================================\n")
cat("CORRECTED B1-03 AUDIT COMPLETE\n")
cat("============================================================\n")
cat(
    "Explicit CRC Stage 4:",
    n_explicit,
    "/",
    n803,
    "\n"
)
cat(
    "30-day qualifying:",
    nrow(q30),
    "\n"
)
cat(
    "60-day qualifying:",
    nrow(q60),
    "\n"
)
cat("\nMain audit:\n", audit_file, "\n")
cat("\nPatient audit:\n", patient_audit_file, "\n")
