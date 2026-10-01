# ============================================================
# A7 / B-line
# 05_B1_05_antiEGFR_outcome_feasibility_audit.R
#
# PURPOSE
#   Audit outcome feasibility for the reframed B-line:
#
#   Baseline EGFR-bypass alterations and outcomes after first
#   anti-EGFR exposure in left-sided RAS/BRAF-WT MSS mCRC.
#
# PRIMARY CLEAN COHORT FOR THIS AUDIT:
#   - OLD_ELIGIBLE_803_R == 1
#   - FIRST_BIOLOGIC_R == "anti-EGFR"
#   - sequencing completed on/before first anti-EGFR start
#
# This automatically requires survival to sequencing before index
# and avoids post-treatment cohort ascertainment for the primary
# anti-EGFR biomarker analysis.
#
# IMPORTANT
#   - This is still an EVENT-FEASIBILITY AUDIT.
#   - No Cox model, HR, KM curve, interaction, rwPFS comparison,
#     TTNTD model, OS model, or treatment-effect estimate is run.
#   - Formal *_R.csv files are not modified.
#
# OUTPUTS
#   06_logs_and_audit/B1_05_antiEGFR_outcome_feasibility_audit.txt
#   06_logs_and_audit/B1_05_antiEGFR_patient_audit.csv
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
cat("B1-05 ANTI-EGFR OUTCOME FEASIBILITY AUDIT\n")
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
        c(0.25, 0.50, 0.75),
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


# ============================================================
# Read main patient-level R-ready CSV
# ============================================================

main <- fread(
    main_file,
    check.names = FALSE,
    showProgress = TRUE
)

required <- c(
    "PATIENT_ID",
    "OLD_ELIGIBLE_803_R",
    "FIRST_BIOLOGIC_R",
    "FIRST_BIOLOGIC_DAY_R",
    "BYPASS_ORIGINAL_R",
    "SPECIMEN_ACQUISITION_FIRST_DAY_R",
    "SEQUENCING_FIRST_DAY_R",
    "OS_MONTHS",
    "OS_STATUS"
)

missing_required <- setdiff(
    required,
    names(main)
)

if (length(missing_required) > 0L) {
    stop(
        "Main CSV missing required columns:\n",
        paste(
            missing_required,
            collapse = "\n"
        )
    )
}


# ============================================================
# Build anti-EGFR candidate cohort
# ============================================================

anti_all <- main[
    OLD_ELIGIBLE_803_R == 1L &
        FIRST_BIOLOGIC_R == "anti-EGFR"
]

if (nrow(anti_all) != 230L) {
    stop(
        "Anti-EGFR anchor failed: expected 230, got ",
        nrow(anti_all)
    )
}

anti_all[
    ,
    SPECIMEN_PRE_ANTI_EGFR_R := as.integer(
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

anti_all[
    ,
    SEQUENCING_PRE_ANTI_EGFR_R := as.integer(
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

primary <- anti_all[
    SEQUENCING_PRE_ANTI_EGFR_R == 1L
]

n_primary <- nrow(primary)

if (n_primary != 191L) {
    stop(
        "Primary clean anti-EGFR cohort anchor failed: ",
        "expected 191, got ",
        n_primary
    )
}

if (
    primary[
        SPECIMEN_PRE_ANTI_EGFR_R != 1L,
        .N
    ] != 0L
) {
    stop(
        "Unexpected: sequencing-pre cohort contains ",
        "post-treatment specimen acquisition."
    )
}

cell_counts <- primary[
    ,
    .N,
    by = BYPASS_ORIGINAL_R
][order(BYPASS_ORIGINAL_R)]

if (
    cell_counts[
        BYPASS_ORIGINAL_R == 0L,
        N
    ] != 156L ||
    cell_counts[
        BYPASS_ORIGINAL_R == 1L,
        N
    ] != 35L
) {
    stop(
        "Bypass cell anchor failed. ",
        "Expected 156 bypass- / 35 bypass+."
    )
}

cat(
    "Primary clean anti-EGFR cohort:",
    n_primary,
    "\n"
)
cat(
    "Bypass- / bypass+:",
    cell_counts[
        BYPASS_ORIGINAL_R == 0L,
        N
    ],
    "/",
    cell_counts[
        BYPASS_ORIGINAL_R == 1L,
        N
    ],
    "\n\n"
)


# ============================================================
# Extract raw progression timeline
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_05_"
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
        stop(
            "Missing raw MSK file: ",
            fn
        )
    }

    hit[1]
}

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
    PATIENT_ID %chin% primary$PATIENT_ID
]

progression[
    ,
    START_DATE_R := suppressWarnings(
        as.numeric(START_DATE)
    )
]

progression_y <- progression[
    PROGRESSION == "Y" &
        is.finite(
            START_DATE_R
        )
]


# ============================================================
# Patient-level outcome audit
#
# Timeline convention in MSK-CHORD:
#   sequencing date approximately = day 0.
#
# OS_MONTHS is defined from first sequenced sample to death or
# last contact. Therefore:
#   FOLLOWUP_DAY_R ~= OS_MONTHS * 30.4375
#
# In the PRIMARY cohort, sequencing is on/before anti-EGFR index,
# so cohort ascertainment is not post-index.
# ============================================================

audit_list <- lapply(
    seq_len(
        nrow(primary)
    ),
    function(i) {

        row <- primary[i]

        pid <- row$PATIENT_ID
        index_day <- suppressWarnings(
            as.numeric(
                row$FIRST_BIOLOGIC_DAY_R
            )
        )

        bypass <- as.integer(
            row$BYPASS_ORIGINAL_R
        )

        os_months <- suppressWarnings(
            as.numeric(
                row$OS_MONTHS
            )
        )

        followup_day <- if (
            is.finite(os_months)
        ) {
            os_months * 30.4375
        } else {
            NA_real_
        }

        deceased <- identical(
            as.character(
                row$OS_STATUS
            ),
            "1:DECEASED"
        )

        death_day <- if (
            deceased &&
            is.finite(
                followup_day
            )
        ) {
            followup_day
        } else {
            NA_real_
        }

        p <- progression_y[
            PATIENT_ID == pid &
                START_DATE_R >=
                    index_day
        ]

        progression_day <- if (
            nrow(p) > 0L
        ) {
            min(
                p$START_DATE_R
            )
        } else {
            NA_real_
        }

        progression_event <- as.integer(
            is.finite(
                progression_day
            )
        )

        death_after_index <- as.integer(
            is.finite(
                death_day
            ) &&
            death_day >=
                index_day
        )

        candidate_event_days <- c(
            progression_day,
            if (
                death_after_index == 1L
            ) {
                death_day
            } else {
                NA_real_
            }
        )

        candidate_event_days <- candidate_event_days[
            is.finite(
                candidate_event_days
            )
        ]

        rwpfs_event <- as.integer(
            length(
                candidate_event_days
            ) > 0L
        )

        rwpfs_event_day <- if (
            rwpfs_event == 1L
        ) {
            min(
                candidate_event_days
            )
        } else {
            NA_real_
        }

        valid_followup <- as.integer(
            is.finite(
                followup_day
            ) &&
            followup_day >=
                index_day
        )

        rwpfs_end_day <- if (
            rwpfs_event == 1L
        ) {
            rwpfs_event_day
        } else if (
            valid_followup == 1L
        ) {
            followup_day
        } else {
            NA_real_
        }

        rwpfs_days <- if (
            is.finite(
                rwpfs_end_day
            )
        ) {
            rwpfs_end_day -
                index_day
        } else {
            NA_real_
        }

        # OS from anti-EGFR index:
        # event if deceased after index;
        # otherwise censor at last contact.
        os_event_from_index <- death_after_index

        os_end_day <- if (
            death_after_index == 1L
        ) {
            death_day
        } else if (
            valid_followup == 1L
        ) {
            followup_day
        } else {
            NA_real_
        }

        os_days_from_index <- if (
            is.finite(
                os_end_day
            )
        ) {
            os_end_day -
                index_day
        } else {
            NA_real_
        }

        data.table(
            PATIENT_ID = pid,
            BYPASS_ORIGINAL_R = bypass,

            ANTI_EGFR_INDEX_DAY_R =
                index_day,

            SPECIMEN_ACQUISITION_FIRST_DAY_R =
                row$SPECIMEN_ACQUISITION_FIRST_DAY_R,

            SEQUENCING_FIRST_DAY_R =
                row$SEQUENCING_FIRST_DAY_R,

            FOLLOWUP_DAY_R =
                followup_day,

            VALID_FOLLOWUP_AFTER_INDEX_R =
                valid_followup,

            PROGRESSION_AFTER_INDEX_R =
                progression_event,

            FIRST_PROGRESSION_DAY_R =
                progression_day,

            DEATH_AFTER_INDEX_R =
                death_after_index,

            DEATH_DAY_R =
                death_day,

            RWPFS_EVENT_R =
                rwpfs_event,

            RWPFS_DAYS_R =
                rwpfs_days,

            OS_EVENT_FROM_INDEX_R =
                os_event_from_index,

            OS_DAYS_FROM_INDEX_R =
                os_days_from_index
        )
    }
)

audit <- rbindlist(
    audit_list,
    fill = TRUE
)


# ============================================================
# Data validity checks
# ============================================================

invalid_followup <- audit[
    VALID_FOLLOWUP_AFTER_INDEX_R != 1L,
    .N
]

negative_rwpfs <- audit[
    is.finite(
        RWPFS_DAYS_R
    ) &
        RWPFS_DAYS_R < 0,
    .N
]

negative_os <- audit[
    is.finite(
        OS_DAYS_FROM_INDEX_R
    ) &
        OS_DAYS_FROM_INDEX_R < 0,
    .N
]


# ============================================================
# Event tables
# ============================================================

event_table <- audit[
    ,
    .(
        N = .N,

        Progression_events =
            sum(
                PROGRESSION_AFTER_INDEX_R,
                na.rm = TRUE
            ),

        Deaths_after_index =
            sum(
                DEATH_AFTER_INDEX_R,
                na.rm = TRUE
            ),

        rwPFS_events =
            sum(
                RWPFS_EVENT_R,
                na.rm = TRUE
            ),

        OS_events =
            sum(
                OS_EVENT_FROM_INDEX_R,
                na.rm = TRUE
            ),

        rwPFS_event_pct =
            round(
                100 *
                    mean(
                        RWPFS_EVENT_R,
                        na.rm = TRUE
                    ),
                1
            ),

        OS_event_pct =
            round(
                100 *
                    mean(
                        OS_EVENT_FROM_INDEX_R,
                        na.rm = TRUE
                    ),
                1
            ),

        rwPFS_days_median_raw =
            median(
                RWPFS_DAYS_R,
                na.rm = TRUE
            ),

        followup_after_index_days_median =
            median(
                FOLLOWUP_DAY_R -
                    ANTI_EGFR_INDEX_DAY_R,
                na.rm = TRUE
            )
    ),
    by = BYPASS_ORIGINAL_R
][order(BYPASS_ORIGINAL_R)]


# ============================================================
# Early events: descriptive quality check only
# ============================================================

early_event_table <- audit[
    ,
    .(
        N = .N,

        rwPFS_event_le_14d =
            sum(
                RWPFS_EVENT_R == 1L &
                    RWPFS_DAYS_R <= 14,
                na.rm = TRUE
            ),

        rwPFS_event_le_30d =
            sum(
                RWPFS_EVENT_R == 1L &
                    RWPFS_DAYS_R <= 30,
                na.rm = TRUE
            ),

        rwPFS_event_le_60d =
            sum(
                RWPFS_EVENT_R == 1L &
                    RWPFS_DAYS_R <= 60,
                na.rm = TRUE
            )
    ),
    by = BYPASS_ORIGINAL_R
][order(BYPASS_ORIGINAL_R)]


# ============================================================
# Event sufficiency flags
# ============================================================

n_bypass_neg <- audit[
    BYPASS_ORIGINAL_R == 0L,
    .N
]

n_bypass_pos <- audit[
    BYPASS_ORIGINAL_R == 1L,
    .N
]

rwpfs_neg_events <- audit[
    BYPASS_ORIGINAL_R == 0L,
    sum(
        RWPFS_EVENT_R,
        na.rm = TRUE
    )
]

rwpfs_pos_events <- audit[
    BYPASS_ORIGINAL_R == 1L,
    sum(
        RWPFS_EVENT_R,
        na.rm = TRUE
    )
]

os_neg_events <- audit[
    BYPASS_ORIGINAL_R == 0L,
    sum(
        OS_EVENT_FROM_INDEX_R,
        na.rm = TRUE
    )
]

os_pos_events <- audit[
    BYPASS_ORIGINAL_R == 1L,
    sum(
        OS_EVENT_FROM_INDEX_R,
        na.rm = TRUE
    )
]

rwpfs_event_feasible <- (
    rwpfs_neg_events >= 30L &&
        rwpfs_pos_events >= 20L
)

os_event_feasible <- (
    os_neg_events >= 30L &&
        os_pos_events >= 20L
)


# ============================================================
# Write patient-level audit
# ============================================================

patient_file <- file.path(
    audit_dir,
    "B1_05_antiEGFR_patient_audit.csv"
)

fwrite(
    audit,
    patient_file,
    na = ""
)


# ============================================================
# Write text audit
# ============================================================

audit_file <- file.path(
    audit_dir,
    "B1_05_antiEGFR_outcome_feasibility_audit.txt"
)

lines <- c(
    "B1-05 ANTI-EGFR OUTCOME FEASIBILITY AUDIT",
    "========================================",
    "",
    "1. Reframed cohort",
    "------------------",
    "Clinical question:",
    "Among left-sided RAS/BRAF-WT MSS mCRC patients receiving first anti-EGFR exposure,",
    "are baseline EGFR-bypass alterations associated with subsequent outcomes?",
    "",
    paste0(
        "All old anti-EGFR first-biologic patients: ",
        nrow(
            anti_all
        )
    ),
    paste0(
        "Primary clean cohort with sequencing <= anti-EGFR index: ",
        n_primary
    ),
    paste0(
        "Bypass-: ",
        n_bypass_neg
    ),
    paste0(
        "Bypass+: ",
        n_bypass_pos
    ),
    paste0(
        "All primary-cohort specimens acquired <= anti-EGFR index: ",
        audit[
            SPECIMEN_ACQUISITION_FIRST_DAY_R <=
                ANTI_EGFR_INDEX_DAY_R,
            .N
        ],
        "/",
        n_primary
    ),
    "",
    "2. Follow-up validity",
    "---------------------",
    paste0(
        "Invalid/missing follow-up after anti-EGFR index: ",
        invalid_followup
    ),
    paste0(
        "Negative rwPFS durations: ",
        negative_rwpfs
    ),
    paste0(
        "Negative OS durations: ",
        negative_os
    ),
    "",
    "3. Outcome-event feasibility by bypass group",
    "--------------------------------------------",
    capture.output(
        print(
            event_table
        )
    ),
    "",
    paste0(
        "rwPFS event sufficiency flag: ",
        rwpfs_event_feasible
    ),
    paste0(
        "OS event sufficiency flag: ",
        os_event_feasible
    ),
    "",
    "4. Early rwPFS-event quality check",
    "----------------------------------",
    capture.output(
        print(
            early_event_table
        )
    ),
    "",
    "5. Interpretation",
    "-----------------",
    "This is an event-feasibility audit only.",
    "No HR, p-value, KM curve, Cox model, treatment-effect estimate, or interaction test was calculated.",
    "",
    "If rwPFS events remain adequate in bypass+ and bypass- groups,",
    "the next step is to freeze the bypass definition and covariate set before any outcome comparison.",
    "",
    "The bevacizumab cohort should NOT be used as a primary active comparator after B1-04.",
    "It may later be used only as an exploratory specificity / negative-control analysis."
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
cat("B1-05 AUDIT COMPLETE\n")
cat("============================================================\n")
cat(
    "Primary clean anti-EGFR cohort:",
    n_primary,
    "\n"
)
cat(
    "Bypass- / bypass+:",
    n_bypass_neg,
    "/",
    n_bypass_pos,
    "\n"
)
cat(
    "rwPFS events bypass- / bypass+:",
    rwpfs_neg_events,
    "/",
    rwpfs_pos_events,
    "\n"
)
cat(
    "OS events bypass- / bypass+:",
    os_neg_events,
    "/",
    os_pos_events,
    "\n"
)
cat("\nMain audit:\n", audit_file, "\n")
cat("\nPatient audit:\n", patient_file, "\n")
