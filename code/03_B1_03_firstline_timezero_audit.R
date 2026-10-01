# ============================================================
# A7 / B-line
# 03_B1_03_firstline_timezero_audit.R
#
# PURPOSE
#   Reconstruct and AUDIT first-line metastatic systemic therapy
#   among the OLD eligible 803 patients.
#
# IMPORTANT
#   - Reads only:
#       A7_MSK_CHORD_CRC_clean_R.csv
#   - Does NOT run rwPFS / OS / treatment-effect models.
#   - Does NOT overwrite the main R-ready CSV.
#   - Does NOT freeze the final first-line/time-zero definition.
#
# OUTPUTS
#   06_logs_and_audit/B1_03_firstline_timezero_audit.txt
#   06_logs_and_audit/B1_03_patient_level_audit.csv
#
# Candidate rules audited:
#   - metastatic diagnosis anchor from diagnosis timeline
#   - first systemic treatment at/after metastatic anchor
#     allowing a 30-day look-back
#   - backbone assembled within 30 days of line start
#   - biologic assignment within 30 days (primary candidate)
#   - biologic assignment within 60 days (sensitivity candidate)
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory
raw_dir <- file.path(project_root, "01_raw_data")
audit_dir <- file.path(project_root, "06_logs_and_audit")

input_file <- file.path(
    raw_dir,
    "A7_MSK_CHORD_CRC_clean_R.csv"
)

if (!file.exists(input_file)) {
    stop("Missing input file:\n", input_file)
}

if (!dir.exists(audit_dir)) {
    dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}
library(data.table)

cat("\n============================================================\n")
cat("B1-03 FIRST-LINE / TIME-ZERO AUDIT\n")
cat("============================================================\n\n")

dt <- fread(
    input_file,
    check.names = FALSE,
    showProgress = TRUE
)

required_cols <- c(
    "PATIENT_ID",
    "OLD_ELIGIBLE_803_R",
    "OLD_BIOLOGIC_458_R",
    "FIRST_BIOLOGIC_R",
    "TREATMENT_SEQUENCE_R",
    "DIAGNOSIS_SEQUENCE_R",
    "SPECIMEN_ACQUISITION_FIRST_DAY_R",
    "SEQUENCING_FIRST_DAY_R"
)

missing_cols <- setdiff(required_cols, names(dt))

if (length(missing_cols) > 0L) {
    stop(
        "Missing required columns:\n",
        paste(missing_cols, collapse = "\n")
    )
}

eligible <- dt[
    OLD_ELIGIBLE_803_R == 1L
]

if (nrow(eligible) != 803L) {
    stop(
        "Anchor failed: OLD_ELIGIBLE_803_R should contain 803 patients; got ",
        nrow(eligible)
    )
}

cat("Old eligible cohort:", nrow(eligible), "\n")


# ============================================================
# Helpers
# ============================================================

clean_token <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "<NA>", "NaN")] <- NA_character_
    x
}

num_or_na <- function(x) {
    suppressWarnings(as.numeric(clean_token(x)))
}

first_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(NA_real_)
    min(x)
}

median_iqr <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) return("NA")
    q <- quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE, names = FALSE)
    sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

pct_text <- function(n, d) {
    if (is.na(d) || d == 0L) return("NA")
    sprintf("%d/%d (%.1f%%)", n, d, 100 * n / d)
}


# ============================================================
# Parse diagnosis timeline
#
# DIAGNOSIS_SEQUENCE_R event format:
# START_DATE ~ SUBTYPE ~ SOURCE ~ DX_DESCRIPTION ~
# STAGE_CDM_DERIVED ~ SUMMARY
# ============================================================

parse_diagnosis <- function(patient_id, seq_string) {

    if (is.na(seq_string) || !nzchar(seq_string)) {
        return(
            data.table(
                PATIENT_ID = patient_id,
                META_ANCHOR_DAY_R = NA_real_,
                META_ANCHOR_SOURCE_R = NA_character_,
                META_DIAG_EVENT_N_R = 0L
            )
        )
    }

    events <- strsplit(
        seq_string,
        "\\|",
        perl = TRUE
    )[[1]]

    out <- rbindlist(
        lapply(
            events,
            function(ev) {

                p <- strsplit(
                    ev,
                    "~",
                    fixed = TRUE
                )[[1]]

                length(p) <- max(length(p), 6L)

                data.table(
                    START_DATE = num_or_na(p[1]),
                    SUBTYPE = clean_token(p[2]),
                    SOURCE = clean_token(p[3]),
                    DX_DESCRIPTION = clean_token(p[4]),
                    STAGE_CDM_DERIVED = clean_token(p[5]),
                    SUMMARY = clean_token(p[6])
                )
            }
        ),
        fill = TRUE
    )

    if (nrow(out) == 0L) {
        return(
            data.table(
                PATIENT_ID = patient_id,
                META_ANCHOR_DAY_R = NA_real_,
                META_ANCHOR_SOURCE_R = NA_character_,
                META_DIAG_EVENT_N_R = 0L
            )
        )
    }

    txt <- tolower(
        paste(
            out$DX_DESCRIPTION,
            out$STAGE_CDM_DERIVED,
            out$SUMMARY
        )
    )

    meta_flag <- grepl(
        "stage[[:space:]]*(4|iv)|distant|metasta",
        txt,
        perl = TRUE
    )

    meta_rows <- out[
        meta_flag %in% TRUE &
            is.finite(START_DATE)
    ]

    if (nrow(meta_rows) == 0L) {
        return(
            data.table(
                PATIENT_ID = patient_id,
                META_ANCHOR_DAY_R = NA_real_,
                META_ANCHOR_SOURCE_R = NA_character_,
                META_DIAG_EVENT_N_R = nrow(out)
            )
        )
    }

    setorder(
        meta_rows,
        START_DATE
    )

    first <- meta_rows[1]

    source_text <- paste(
        na.omit(
            c(
                first$STAGE_CDM_DERIVED,
                first$SUMMARY,
                first$DX_DESCRIPTION
            )
        ),
        collapse = " | "
    )

    data.table(
        PATIENT_ID = patient_id,
        META_ANCHOR_DAY_R = first$START_DATE,
        META_ANCHOR_SOURCE_R = source_text,
        META_DIAG_EVENT_N_R = nrow(out)
    )
}


diag_audit <- rbindlist(
    Map(
        parse_diagnosis,
        eligible$PATIENT_ID,
        eligible$DIAGNOSIS_SEQUENCE_R
    ),
    fill = TRUE
)


# ============================================================
# Parse treatment timeline
#
# TREATMENT_SEQUENCE_R event format:
# START_DATE ~ STOP_DATE ~ SUBTYPE ~ AGENT ~ RX_INVESTIGATIVE
# ============================================================

parse_treatment <- function(patient_id, seq_string) {

    if (is.na(seq_string) || !nzchar(seq_string)) {
        return(
            data.table(
                PATIENT_ID = character(),
                START_DATE = numeric(),
                STOP_DATE = numeric(),
                SUBTYPE = character(),
                AGENT = character(),
                RX_INVESTIGATIVE = character()
            )
        )
    }

    events <- strsplit(
        seq_string,
        "\\|",
        perl = TRUE
    )[[1]]

    out <- rbindlist(
        lapply(
            events,
            function(ev) {

                p <- strsplit(
                    ev,
                    "~",
                    fixed = TRUE
                )[[1]]

                length(p) <- max(length(p), 5L)

                data.table(
                    PATIENT_ID = patient_id,
                    START_DATE = num_or_na(p[1]),
                    STOP_DATE = num_or_na(p[2]),
                    SUBTYPE = clean_token(p[3]),
                    AGENT = clean_token(p[4]),
                    RX_INVESTIGATIVE = clean_token(p[5])
                )
            }
        ),
        fill = TRUE
    )

    out[
        is.finite(START_DATE)
    ]
}


tx_long <- rbindlist(
    Map(
        parse_treatment,
        eligible$PATIENT_ID,
        eligible$TREATMENT_SEQUENCE_R
    ),
    fill = TRUE
)

if (nrow(tx_long) == 0L) {
    stop("No treatment events parsed from TREATMENT_SEQUENCE_R.")
}

tx_long[
    ,
    AGENT_UPPER := toupper(
        trimws(
            fifelse(
                is.na(AGENT),
                "",
                AGENT
            )
        )
    )
]

# Agent families
tx_long[, IS_5FU := AGENT_UPPER %chin% c(
    "FLUOROURACIL",
    "5-FLUOROURACIL",
    "5-FU"
)]

tx_long[, IS_CAPE := AGENT_UPPER == "CAPECITABINE"]
tx_long[, IS_OX := AGENT_UPPER == "OXALIPLATIN"]
tx_long[, IS_IRI := AGENT_UPPER == "IRINOTECAN"]
tx_long[, IS_CETUX := AGENT_UPPER == "CETUXIMAB"]
tx_long[, IS_PANI := AGENT_UPPER == "PANITUMUMAB"]
tx_long[, IS_BEV := AGENT_UPPER == "BEVACIZUMAB"]

tx_long[, IS_FP := IS_5FU | IS_CAPE]
tx_long[, IS_EGFR := IS_CETUX | IS_PANI]
tx_long[, IS_BIOLOGIC := IS_EGFR | IS_BEV]


# ============================================================
# Reconstruct candidate first metastatic line
#
# Strict candidate:
#   earliest treatment START_DATE >= metastatic anchor - 30 days.
#
# Fallback candidate:
#   if metastatic anchor absent, earliest treatment START_DATE.
#
# We keep strict and fallback flags separately.
# ============================================================

base <- merge(
    eligible[
        ,
        .(
            PATIENT_ID,
            OLD_BIOLOGIC_458_R,
            OLD_FIRST_BIOLOGIC_R = FIRST_BIOLOGIC_R,
            SPECIMEN_ACQUISITION_FIRST_DAY_R,
            SEQUENCING_FIRST_DAY_R
        )
    ],
    diag_audit,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

firstline_list <- lapply(
    seq_len(nrow(base)),
    function(i) {

        pid <- base$PATIENT_ID[i]
        meta_day <- base$META_ANCHOR_DAY_R[i]

        ptx <- tx_long[
            PATIENT_ID == pid
        ]

        if (nrow(ptx) == 0L) {

            return(
                data.table(
                    PATIENT_ID = pid,
                    TX_EVENT_N_R = 0L,
                    META_ANCHOR_AVAILABLE_R = as.integer(
                        is.finite(meta_day)
                    ),
                    LINE_START_DAY_R = NA_real_,
                    LINE_START_RULE_R = NA_character_,
                    BACKBONE_30D_R = NA_character_,
                    BIOLOGIC_30D_R = NA_character_,
                    BIOLOGIC_60D_R = NA_character_,
                    BIOLOGIC_FIRST_DAY_R = NA_real_,
                    BIOLOGIC_DELAY_DAYS_R = NA_real_,
                    DUAL_BIOLOGIC_30D_R = 0L,
                    DUAL_BIOLOGIC_60D_R = 0L,
                    TRIPLET_30D_R = 0L,
                    DOUBLEt_30D_R = 0L,
                    QUALIFY_30D_R = 0L,
                    QUALIFY_60D_R = 0L,
                    LINE_AGENTS_30D_R = NA_character_,
                    LINE_AGENTS_60D_R = NA_character_
                )
            )
        }

        setorder(
            ptx,
            START_DATE
        )

        if (is.finite(meta_day)) {

            candidates <- ptx[
                START_DATE >= (meta_day - 30)
            ]

            if (nrow(candidates) == 0L) {
                line_start <- NA_real_
                line_rule <- "no treatment after metastatic anchor"
            } else {
                line_start <- min(candidates$START_DATE)
                line_rule <- "metastatic anchor -30d"
            }

        } else {

            line_start <- min(ptx$START_DATE)
            line_rule <- "fallback: earliest treatment"
        }

        if (!is.finite(line_start)) {

            return(
                data.table(
                    PATIENT_ID = pid,
                    TX_EVENT_N_R = nrow(ptx),
                    META_ANCHOR_AVAILABLE_R = as.integer(
                        is.finite(meta_day)
                    ),
                    LINE_START_DAY_R = NA_real_,
                    LINE_START_RULE_R = line_rule,
                    BACKBONE_30D_R = NA_character_,
                    BIOLOGIC_30D_R = NA_character_,
                    BIOLOGIC_60D_R = NA_character_,
                    BIOLOGIC_FIRST_DAY_R = NA_real_,
                    BIOLOGIC_DELAY_DAYS_R = NA_real_,
                    DUAL_BIOLOGIC_30D_R = 0L,
                    DUAL_BIOLOGIC_60D_R = 0L,
                    TRIPLET_30D_R = 0L,
                    DOUBLEt_30D_R = 0L,
                    QUALIFY_30D_R = 0L,
                    QUALIFY_60D_R = 0L,
                    LINE_AGENTS_30D_R = NA_character_,
                    LINE_AGENTS_60D_R = NA_character_
                )
            )
        }

        w30 <- ptx[
            START_DATE >= line_start &
                START_DATE <= line_start + 30
        ]

        w60 <- ptx[
            START_DATE >= line_start &
                START_DATE <= line_start + 60
        ]

        has_fp30 <- any(w30$IS_FP)
        has_ox30 <- any(w30$IS_OX)
        has_iri30 <- any(w30$IS_IRI)

        triplet30 <- (
            has_fp30 &&
                has_ox30 &&
                has_iri30
        )

        doublet30 <- (
            has_fp30 &&
                xor(has_ox30, has_iri30)
        )

        backbone30 <- if (triplet30) {

            "triplet/FOLFOXIRI-like"

        } else if (has_fp30 && has_ox30 && !has_iri30) {

            if (any(w30$IS_CAPE)) {
                "CAPOX/FOLFOX-like"
            } else {
                "FOLFOX-like"
            }

        } else if (has_fp30 && has_iri30 && !has_ox30) {

            "FOLFIRI-like"

        } else {

            "other/non-doublet"
        }

        egfr30 <- any(w30$IS_EGFR)
        bev30 <- any(w30$IS_BEV)

        egfr60 <- any(w60$IS_EGFR)
        bev60 <- any(w60$IS_BEV)

        biologic30 <- if (egfr30 && bev30) {
            "dual"
        } else if (egfr30) {
            "anti-EGFR"
        } else if (bev30) {
            "bevacizumab"
        } else {
            "none"
        }

        biologic60 <- if (egfr60 && bev60) {
            "dual"
        } else if (egfr60) {
            "anti-EGFR"
        } else if (bev60) {
            "bevacizumab"
        } else {
            "none"
        }

        bio_rows <- ptx[
            IS_BIOLOGIC &
                START_DATE >= line_start &
                START_DATE <= line_start + 60
        ]

        bio_first <- if (nrow(bio_rows) > 0L) {
            min(bio_rows$START_DATE)
        } else {
            NA_real_
        }

        bio_delay <- if (is.finite(bio_first)) {
            bio_first - line_start
        } else {
            NA_real_
        }

        qualify30 <- as.integer(
            doublet30 &&
                biologic30 %chin% c(
                    "anti-EGFR",
                    "bevacizumab"
                )
        )

        qualify60 <- as.integer(
            doublet30 &&
                biologic60 %chin% c(
                    "anti-EGFR",
                    "bevacizumab"
                )
        )

        data.table(
            PATIENT_ID = pid,
            TX_EVENT_N_R = nrow(ptx),

            META_ANCHOR_AVAILABLE_R = as.integer(
                is.finite(meta_day)
            ),

            LINE_START_DAY_R = line_start,
            LINE_START_RULE_R = line_rule,

            BACKBONE_30D_R = backbone30,
            BIOLOGIC_30D_R = biologic30,
            BIOLOGIC_60D_R = biologic60,

            BIOLOGIC_FIRST_DAY_R = bio_first,
            BIOLOGIC_DELAY_DAYS_R = bio_delay,

            DUAL_BIOLOGIC_30D_R = as.integer(
                egfr30 && bev30
            ),

            DUAL_BIOLOGIC_60D_R = as.integer(
                egfr60 && bev60
            ),

            TRIPLET_30D_R = as.integer(
                triplet30
            ),

            DOUBLEt_30D_R = as.integer(
                doublet30
            ),

            QUALIFY_30D_R = qualify30,
            QUALIFY_60D_R = qualify60,

            LINE_AGENTS_30D_R = paste(
                sort(
                    unique(
                        w30$AGENT[
                            !is.na(w30$AGENT)
                        ]
                    )
                ),
                collapse = "|"
            ),

            LINE_AGENTS_60D_R = paste(
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
)

firstline <- rbindlist(
    firstline_list,
    fill = TRUE
)

audit <- merge(
    base,
    firstline,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)


# ============================================================
# Baseline timing relative to candidate time-zero
# ============================================================

audit[
    ,
    SPECIMEN_PRE_LINE_R := fifelse(
        is.finite(SPECIMEN_ACQUISITION_FIRST_DAY_R) &
            is.finite(LINE_START_DAY_R),
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
        is.finite(SEQUENCING_FIRST_DAY_R) &
            is.finite(LINE_START_DAY_R),
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
        is.finite(SPECIMEN_ACQUISITION_FIRST_DAY_R) &
            is.finite(BIOLOGIC_FIRST_DAY_R),
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
        is.finite(SEQUENCING_FIRST_DAY_R) &
            is.finite(BIOLOGIC_FIRST_DAY_R),
        as.integer(
            SEQUENCING_FIRST_DAY_R <=
                BIOLOGIC_FIRST_DAY_R
        ),
        NA_integer_
    )
]


# ============================================================
# Counts
# ============================================================

n803 <- nrow(audit)
n_meta_anchor <- audit[META_ANCHOR_AVAILABLE_R == 1L, .N]
n_tx <- audit[TX_EVENT_N_R > 0L, .N]

backbone_tab <- audit[
    ,
    .N,
    by = BACKBONE_30D_R
][order(-N)]

bio30_tab <- audit[
    ,
    .N,
    by = BIOLOGIC_30D_R
][order(-N)]

bio60_tab <- audit[
    ,
    .N,
    by = BIOLOGIC_60D_R
][order(-N)]

n_q30 <- audit[QUALIFY_30D_R == 1L, .N]
n_q60 <- audit[QUALIFY_60D_R == 1L, .N]

q30_egfr <- audit[
    QUALIFY_30D_R == 1L &
        BIOLOGIC_30D_R == "anti-EGFR",
    .N
]

q30_bev <- audit[
    QUALIFY_30D_R == 1L &
        BIOLOGIC_30D_R == "bevacizumab",
    .N
]

q60_egfr <- audit[
    QUALIFY_60D_R == 1L &
        BIOLOGIC_60D_R == "anti-EGFR",
    .N
]

q60_bev <- audit[
    QUALIFY_60D_R == 1L &
        BIOLOGIC_60D_R == "bevacizumab",
    .N
]

n_triplet <- audit[
    TRIPLET_30D_R == 1L,
    .N
]

n_dual30 <- audit[
    DUAL_BIOLOGIC_30D_R == 1L,
    .N
]

n_dual60 <- audit[
    DUAL_BIOLOGIC_60D_R == 1L,
    .N
]

old458_overlap30 <- audit[
    OLD_BIOLOGIC_458_R == 1L &
        QUALIFY_30D_R == 1L,
    .N
]

old458_overlap60 <- audit[
    OLD_BIOLOGIC_458_R == 1L &
        QUALIFY_60D_R == 1L,
    .N
]

delayed_31_60 <- audit[
    is.finite(BIOLOGIC_DELAY_DAYS_R) &
        BIOLOGIC_DELAY_DAYS_R > 30 &
        BIOLOGIC_DELAY_DAYS_R <= 60,
    .N
]


# ============================================================
# Baseline timing counts among qualifying patients
# ============================================================

baseline_block <- function(x, label) {

    d <- audit[
        QUALIFY_30D_R == 1L &
            !is.na(get(x))
    ]

    n_yes <- d[
        get(x) == 1L,
        .N
    ]

    paste0(
        label,
        ": ",
        pct_text(
            n_yes,
            nrow(d)
        )
    )
}

baseline_line_1 <- baseline_block(
    "SPECIMEN_PRE_LINE_R",
    "Specimen acquisition <= line start"
)

baseline_line_2 <- baseline_block(
    "SEQUENCING_PRE_LINE_R",
    "Sequencing <= line start"
)

baseline_line_3 <- baseline_block(
    "SPECIMEN_PRE_BIOLOGIC_R",
    "Specimen acquisition <= biologic start"
)

baseline_line_4 <- baseline_block(
    "SEQUENCING_PRE_BIOLOGIC_R",
    "Sequencing <= biologic start"
)


# ============================================================
# Write patient-level audit CSV
# This is an audit artifact, NOT an R-ready analysis CSV.
# ============================================================

patient_audit_file <- file.path(
    audit_dir,
    "B1_03_patient_level_audit.csv"
)

fwrite(
    audit,
    patient_audit_file,
    na = ""
)


# ============================================================
# Write text audit
# ============================================================

audit_file <- file.path(
    audit_dir,
    "B1_03_firstline_timezero_audit.txt"
)

lines <- c(
    "B1-03 FIRST-LINE / TIME-ZERO AUDIT",
    "=================================",
    "",
    "Scope",
    "-----",
    paste0("Old eligible cohort: ", n803),
    "No rwPFS / OS / treatment-effect estimate was calculated.",
    "",
    "Metastatic anchor / treatment availability",
    "------------------------------------------",
    paste0(
        "Metastatic diagnosis anchor available: ",
        pct_text(n_meta_anchor, n803)
    ),
    paste0(
        "At least one parsed treatment event: ",
        pct_text(n_tx, n803)
    ),
    "",
    "Candidate first-line backbone within 30 days",
    "--------------------------------------------",
    capture.output(print(backbone_tab)),
    "",
    "Biologic within 30 days",
    "-----------------------",
    capture.output(print(bio30_tab)),
    "",
    "Biologic within 60 days",
    "-----------------------",
    capture.output(print(bio60_tab)),
    "",
    "Qualifying candidate first-line doublet + single biologic",
    "---------------------------------------------------------",
    paste0(
        "30-day definition: ",
        n_q30,
        " total; anti-EGFR ",
        q30_egfr,
        "; bevacizumab ",
        q30_bev
    ),
    paste0(
        "60-day definition: ",
        n_q60,
        " total; anti-EGFR ",
        q60_egfr,
        "; bevacizumab ",
        q60_bev
    ),
    "",
    "Important exclusions / timing",
    "-----------------------------",
    paste0(
        "Triplet/FOLFOXIRI-like within 30 days: ",
        n_triplet
    ),
    paste0(
        "Dual anti-EGFR + BEV within 30 days: ",
        n_dual30
    ),
    paste0(
        "Dual anti-EGFR + BEV within 60 days: ",
        n_dual60
    ),
    paste0(
        "Biologic first added on days 31-60: ",
        delayed_31_60
    ),
    paste0(
        "Biologic delay median [Q1, Q3], days: ",
        median_iqr(
            audit$BIOLOGIC_DELAY_DAYS_R
        )
    ),
    "",
    "Relation to OLD 458 feasibility cohort",
    "--------------------------------------",
    paste0(
        "OLD 458 overlapping 30-day qualifying first-line cohort: ",
        old458_overlap30,
        "/458"
    ),
    paste0(
        "OLD 458 overlapping 60-day qualifying first-line cohort: ",
        old458_overlap60,
        "/458"
    ),
    "",
    "Baseline timing among 30-day qualifying candidates",
    "--------------------------------------------------",
    baseline_line_1,
    baseline_line_2,
    baseline_line_3,
    baseline_line_4,
    "",
    "Interpretation rule",
    "-------------------",
    "This file is an AUDIT only.",
    "Do not freeze line-start versus biologic-start time zero until these counts are reviewed.",
    "Do not run treatment-effect models from B1-03 alone."
)

writeLines(
    lines,
    audit_file
)

cat("\n============================================================\n")
cat("B1-03 AUDIT COMPLETE\n")
cat("============================================================\n")
cat("Old eligible:", n803, "\n")
cat(
    "30-day qualifying:",
    n_q30,
    "(anti-EGFR",
    q30_egfr,
    "/ BEV",
    q30_bev,
    ")\n"
)
cat(
    "60-day qualifying:",
    n_q60,
    "(anti-EGFR",
    q60_egfr,
    "/ BEV",
    q60_bev,
    ")\n"
)
cat("\nMain audit file:\n", audit_file, "\n")
cat("\nPatient-level audit:\n", patient_audit_file, "\n")
