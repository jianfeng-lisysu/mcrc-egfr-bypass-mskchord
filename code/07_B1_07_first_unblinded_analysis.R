# ============================================================
# A7 / B-line
# 07_B1_07_first_unblinded_analysis.R
#
# FIRST FORMAL UNBLINDED OUTCOME ANALYSIS
#
# Frozen cohort:
#   left-sided, RAS-WT, BRAF V600E-negative, MSS mCRC
#   first anti-EGFR exposure
#   sequencing completed on/before anti-EGFR index
#
# Frozen PRIMARY exposure:
#   BYPASS_HIGH_CONFIDENCE_R from B1-06
#
# Primary outcome:
#   rwPFS from first anti-EGFR exposure
#   event = first documented progression STRICTLY AFTER index
#           OR death after index
#   censor = last follow-up
#
# Secondary outcome:
#   OS from first anti-EGFR exposure
#
# Primary adjusted Cox covariates (frozen before outcome):
#   sex
#   frozen left subsite
#   anti-EGFR agent
#   log(Stage IV -> anti-EGFR interval)
#   prior exposure to BOTH oxaliplatin and irinotecan
#
# Pre-specified exposure sensitivities:
#   - PRESSING-like strict
#   - high-confidence excluding MAP2K1 D67E
#   - original broad composite
#
# Pre-specified cohort sensitivities:
#   - exclude BRAF class III
#   - specimen acquired <=730 days before anti-EGFR
#   - exclude structured-vs-diagnosis sidedness conflict
#   - BRAF class III exclusion + specimen <=730 days
#
# IMPORTANT:
#   This script does NOT modify the two formal *_R.csv files.
#
# Outputs:
#   03_intermediate/A7_B1_locked_analysis.rds
#   04_results/B1_07_primary_models.rds
#   06_logs_and_audit/B1_07_first_unblinded_results.txt
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

raw_dir <- file.path(
    project_root,
    "01_raw_data"
)

intermediate_dir <- file.path(
    project_root,
    "03_intermediate"
)

results_dir <- file.path(
    project_root,
    "04_results"
)

audit_dir <- file.path(
    project_root,
    "06_logs_and_audit"
)

for (d in c(
    intermediate_dir,
    results_dir,
    audit_dir
)) {
    if (!dir.exists(d)) {
        dir.create(
            d,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }
}

b106_file <- file.path(
    audit_dir,
    "B1_06_patient_genomic_covariate_audit.csv"
)

b105_file <- file.path(
    audit_dir,
    "B1_05_antiEGFR_patient_audit.csv"
)

tar_file <- file.path(
    raw_dir,
    "msk_chord_2024.tar"
)

for (f in c(
    b106_file,
    b105_file,
    tar_file
)) {
    if (!file.exists(f)) {
        stop(
            "Missing required input:\n",
            f
        )
    }
}

if (!requireNamespace(
    "data.table",
    quietly = TRUE
)) {
    install.packages(
        "data.table"
    )
}

if (!requireNamespace(
    "survival",
    quietly = TRUE
)) {
    install.packages(
        "survival"
    )
}

library(data.table)
library(survival)

cat("\n============================================================\n")
cat("B1-07 FIRST FORMAL UNBLINDED ANALYSIS\n")
cat("============================================================\n\n")


# ============================================================
# Helpers
# ============================================================

safe_chr <- function(x) {

    x <- as.character(x)
    x[is.na(x)] <- ""

    x
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

fmt_num <- function(x, digits = 2) {

    if (
        length(x) == 0L ||
        is.na(x) ||
        !is.finite(x)
    ) {
        return("NA")
    }

    formatC(
        x,
        format = "f",
        digits = digits
    )
}

fmt_p <- function(x) {

    if (
        length(x) == 0L ||
        is.na(x) ||
        !is.finite(x)
    ) {
        return("NA")
    }

    if (x < 0.001) {
        return("<0.001")
    }

    formatC(
        x,
        format = "f",
        digits = 3
    )
}

extract_site_codes <- function(x) {

    txt <- toupper(
        safe_chr(x)
    )

    m <- gregexpr(
        "C[0-9]{3}",
        txt,
        perl = TRUE
    )

    regmatches(
        txt,
        m
    )
}

cox_result <- function(
    fit,
    term,
    label,
    outcome,
    analysis
) {

    s <- summary(fit)

    if (!term %in%
        rownames(
            s$coefficients
        )) {

        stop(
            "Term not found in Cox model: ",
            term
        )
    }

    b <- s$coefficients[
        term,
        "coef"
    ]

    se <- s$coefficients[
        term,
        "se(coef)"
    ]

    p <- s$coefficients[
        term,
        "Pr(>|z|)"
    ]

    data.table(
        OUTCOME = outcome,
        ANALYSIS = analysis,
        EXPOSURE = label,
        N = fit$n,
        EVENTS = fit$nevent,
        HR = exp(b),
        LCL95 = exp(
            b - 1.96 * se
        ),
        UCL95 = exp(
            b + 1.96 * se
        ),
        P = p
    )
}

km_median <- function(
    fit,
    group_value
) {

    tab <- summary(fit)$table

    if (
        is.null(
            dim(tab)
        )
    ) {
        return(
            as.numeric(
                tab["median"]
            )
        )
    }

    rn <- rownames(tab)

    hit <- grep(
        paste0(
            "=",
            group_value,
            "$"
        ),
        rn
    )

    if (length(hit) != 1L) {
        return(NA_real_)
    }

    as.numeric(
        tab[
            hit,
            "median"
        ]
    )
}

surv_at_time <- function(
    fit,
    time
) {

    s <- summary(
        fit,
        times = time,
        extend = TRUE
    )

    data.table(
        STRATA = as.character(
            s$strata
        ),
        TIME = s$time,
        SURV = s$surv,
        LOWER = s$lower,
        UPPER = s$upper
    )
}


# ============================================================
# 1. Read frozen genomic/covariate audit and outcome scaffold
# ============================================================

b106 <- fread(
    b106_file,
    check.names = FALSE,
    showProgress = FALSE
)

b105 <- fread(
    b105_file,
    check.names = FALSE,
    showProgress = FALSE
)

required_106 <- c(
    "PATIENT_ID",
    "FIRST_BIOLOGIC_DAY_R",
    "BYPASS_PRESSING_STRICT_R",
    "BYPASS_HIGH_CONFIDENCE_R",
    "BYPASS_HIGH_CONFIDENCE_NO_D67E_R",
    "BYPASS_ORIGINAL_BROAD_R",
    "BRAF_CLASS3_R",
    "SPECIMEN_WITHIN_730D_R",
    "GENDER",
    "PRIMARY_SITE",
    "ANTI_EGFR_AGENT_R",
    "LOG_STAGE4_TO_ANTI_EGFR_R",
    "PRIOR_HEAVY_R",
    "SAMPLE_TYPE",
    "PANEL_OLD_R",
    "MET_SITE_COUNT_90D_R"
)

required_105 <- c(
    "PATIENT_ID",
    "FOLLOWUP_DAY_R",
    "DEATH_AFTER_INDEX_R",
    "DEATH_DAY_R"
)

missing_106 <- setdiff(
    required_106,
    names(b106)
)

missing_105 <- setdiff(
    required_105,
    names(b105)
)

if (length(missing_106) > 0L) {

    stop(
        "B1-06 audit missing:\n",
        paste(
            missing_106,
            collapse = "\n"
        )
    )
}

if (length(missing_105) > 0L) {

    stop(
        "B1-05 audit missing:\n",
        paste(
            missing_105,
            collapse = "\n"
        )
    )
}

dt <- merge(
    b106,
    b105[
        ,
        ..required_105
    ],
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

if (nrow(dt) != 191L) {

    stop(
        "Locked cohort should contain 191 patients; got ",
        nrow(dt)
    )
}


# ============================================================
# 2. Extract original diagnosis + progression timelines
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_07_"
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
        basename(all_files) ==
            filename
    ]

    if (length(hit) == 0L) {

        stop(
            "Missing raw MSK file: ",
            filename
        )
    }

    hit[1]
}

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
        dt$PATIENT_ID
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
    SITE_CODES_LIST_R :=
        extract_site_codes(
            DX_DESCRIPTION
        )
]

# Flatten site code to last ICD-O site code in each diagnosis row
diagnosis[
    ,
    SITE_CODE_R :=
        vapply(
            SITE_CODES_LIST_R,
            function(z) {

                if (
                    length(z) == 0L
                ) {
                    NA_character_
                } else {
                    tail(
                        z,
                        1L
                    )
                }
            },
            character(1)
        )
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


# ============================================================
# 3. Correct/freeze left subsite using:
#    structured specific site first, diagnosis salvage second.
# ============================================================

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
        DX_SIDE_R = dx_side,
        DX_SUBSITE_R =
            dx_subsite
    )
}

dx_patient <- rbindlist(
    lapply(
        dt$PATIENT_ID,
        derive_dx_patient
    )
)

dt <- merge(
    dt,
    dx_patient,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

primary_site_lower <- tolower(
    trimws(
        safe_chr(
            dt$PRIMARY_SITE
        )
    )
)

dt[
    ,
    STRUCTURED_SIDE_R := fcase(
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

dt[
    ,
    STRUCTURED_SUBSITE_R := fcase(
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

dt[
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

dt[
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

if (
    dt[
        is.na(
            SUBSITE_FROZEN_R
        ),
        .N
    ] != 0L
) {

    stop(
        "Frozen subsite contains missing values."
    )
}

subsite_counts <- dt[
    ,
    .N,
    by = SUBSITE_FROZEN_R
][order(
    SUBSITE_FROZEN_R
)]

n_side_conflict <- dt[
    SIDE_CONFLICT_R == 1L,
    .N
]

# Expected outcome-blind anchors from independent review
if (
    dt[
        SUBSITE_FROZEN_R ==
            "left colon",
        .N
    ] != 96L ||
    dt[
        SUBSITE_FROZEN_R ==
            "rectum/rectosigmoid",
        .N
    ] != 95L ||
    n_side_conflict != 1L
) {

    stop(
        paste0(
            "Subsite hard audit failed. ",
            "Expected left colon / rectum / conflict = 96 / 95 / 1."
        )
    )
}


# ============================================================
# 4. Reconstruct rwPFS from raw progression timeline
#    STRICTLY AFTER anti-EGFR index.
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
    PATIENT_ID %chin%
        dt$PATIENT_ID
]

progression[
    ,
    START_DATE_NUM_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

progression_y <- progression[
    PROGRESSION == "Y" &
        is.finite(
            START_DATE_NUM_R
        )
]

build_outcome <- function(
    pid,
    index_day,
    followup_day,
    death_after_index,
    death_day
) {

    p <- progression_y[
        PATIENT_ID == pid &
            START_DATE_NUM_R >
                index_day &
            START_DATE_NUM_R <=
                followup_day
    ]

    prog_day <- if (
        nrow(p) > 0L
    ) {
        min(
            p$START_DATE_NUM_R
        )
    } else {
        NA_real_
    }

    death_valid <- (
        !is.na(
            death_after_index
        ) &&
        death_after_index == 1L &&
        is.finite(
            death_day
        ) &&
        death_day >
            index_day &&
        death_day <=
            followup_day
    )

    death_event_day <- if (
        death_valid
    ) {
        death_day
    } else {
        NA_real_
    }

    event_days <- c(
        prog_day,
        death_event_day
    )

    event_days <- event_days[
        is.finite(
            event_days
        )
    ]

    rwpfs_event <- as.integer(
        length(
            event_days
        ) > 0L
    )

    rwpfs_end <- if (
        rwpfs_event == 1L
    ) {
        min(
            event_days
        )
    } else {
        followup_day
    }

    os_event <- as.integer(
        death_valid
    )

    os_end <- if (
        os_event == 1L
    ) {
        death_event_day
    } else {
        followup_day
    }

    data.table(
        PATIENT_ID = pid,

        FIRST_PROGRESSION_AFTER_INDEX_R =
            prog_day,

        RWPFS_EVENT_LOCKED_R =
            rwpfs_event,

        RWPFS_DAYS_LOCKED_R =
            rwpfs_end -
            index_day,

        OS_EVENT_LOCKED_R =
            os_event,

        OS_DAYS_LOCKED_R =
            os_end -
            index_day
    )
}

outcome_locked <- rbindlist(
    Map(
        build_outcome,
        dt$PATIENT_ID,
        dt$FIRST_BIOLOGIC_DAY_R,
        dt$FOLLOWUP_DAY_R,
        dt$DEATH_AFTER_INDEX_R,
        dt$DEATH_DAY_R
    ),
    fill = TRUE
)

dt <- merge(
    dt,
    outcome_locked,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

if (
    dt[
        !is.finite(
            RWPFS_DAYS_LOCKED_R
        ) |
        RWPFS_DAYS_LOCKED_R <= 0,
        .N
    ] != 0L
) {

    stop(
        "Locked rwPFS contains non-positive or missing durations."
    )
}

if (
    dt[
        !is.finite(
            OS_DAYS_LOCKED_R
        ) |
        OS_DAYS_LOCKED_R <= 0,
        .N
    ] != 0L
) {

    stop(
        "Locked OS contains non-positive or missing durations."
    )
}


# ============================================================
# 5. Hard analysis anchors
# ============================================================

anchors_ok <- (
    nrow(dt) == 191L &&
        dt[
            BYPASS_HIGH_CONFIDENCE_R == 1L,
            .N
        ] == 23L &&
        dt[
            BYPASS_PRESSING_STRICT_R == 1L,
            .N
        ] == 15L &&
        dt[
            BYPASS_HIGH_CONFIDENCE_NO_D67E_R == 1L,
            .N
        ] == 22L &&
        dt[
            BYPASS_ORIGINAL_BROAD_R == 1L,
            .N
        ] == 35L &&
        dt[
            BRAF_CLASS3_R == 1L,
            .N
        ] == 7L &&
        dt[
            SPECIMEN_WITHIN_730D_R == 1L,
            .N
        ] == 149L &&
        dt[
            RWPFS_EVENT_LOCKED_R == 1L,
            .N
        ] == 183L &&
        dt[
            OS_EVENT_LOCKED_R == 1L,
            .N
        ] == 129L
)

if (!anchors_ok) {

    stop(
        "B1-07 hard analysis anchors failed."
    )
}

cat(
    "Hard anchors PASS: ",
    "N=191; high-confidence=23; strict=15; noD67E=22; ",
    "broad=35; BRAF3=7; specimen<=730d=149; ",
    "rwPFS events=183; OS events=129\n",
    sep = ""
)


# ============================================================
# 6. Freeze factor coding
# ============================================================

dt[
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

dt[
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

dt[
    ,
    ANTI_EGFR_AGENT_F_R :=
        factor(
            ANTI_EGFR_AGENT_R,
            levels = c(
                "CETUXIMAB",
                "PANITUMUMAB"
            )
        )
]

dt[
    ,
    PRIOR_HEAVY_F_R :=
        factor(
            PRIOR_HEAVY_R,
            levels = c(
                0,
                1
            )
        )
]


# ============================================================
# 7. PRIMARY rwPFS analysis
# ============================================================

surv_rwpfs <- Surv(
    time =
        dt$RWPFS_DAYS_LOCKED_R,
    event =
        dt$RWPFS_EVENT_LOCKED_R
)

km_rwpfs <- survfit(
    surv_rwpfs ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = dt,
    conf.type = "log-log"
)

logrank_rwpfs <- survdiff(
    surv_rwpfs ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = dt
)

logrank_rwpfs_p <- pchisq(
    logrank_rwpfs$chisq,
    df = 1,
    lower.tail = FALSE
)

cox_rwpfs_crude <- coxph(
    surv_rwpfs ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = dt,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

cox_rwpfs_adj <- coxph(
    surv_rwpfs ~
        BYPASS_HIGH_CONFIDENCE_R +
        GENDER_F_R +
        SUBSITE_F_R +
        ANTI_EGFR_AGENT_F_R +
        LOG_STAGE4_TO_ANTI_EGFR_R +
        PRIOR_HEAVY_F_R,
    data = dt,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

ph_rwpfs <- cox.zph(
    cox_rwpfs_adj,
    transform = "km"
)


# ============================================================
# 8. SECONDARY OS analysis
# ============================================================

surv_os <- Surv(
    time =
        dt$OS_DAYS_LOCKED_R,
    event =
        dt$OS_EVENT_LOCKED_R
)

km_os <- survfit(
    surv_os ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = dt,
    conf.type = "log-log"
)

logrank_os <- survdiff(
    surv_os ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = dt
)

logrank_os_p <- pchisq(
    logrank_os$chisq,
    df = 1,
    lower.tail = FALSE
)

cox_os_crude <- coxph(
    surv_os ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = dt,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

cox_os_adj <- coxph(
    surv_os ~
        BYPASS_HIGH_CONFIDENCE_R +
        GENDER_F_R +
        SUBSITE_F_R +
        ANTI_EGFR_AGENT_F_R +
        LOG_STAGE4_TO_ANTI_EGFR_R +
        PRIOR_HEAVY_F_R,
    data = dt,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

ph_os <- cox.zph(
    cox_os_adj,
    transform = "km"
)


# ============================================================
# 9. Extract primary model results
# ============================================================

primary_results <- rbindlist(
    list(
        cox_result(
            cox_rwpfs_crude,
            "BYPASS_HIGH_CONFIDENCE_R",
            "High-confidence bypass",
            "rwPFS",
            "Crude"
        ),

        cox_result(
            cox_rwpfs_adj,
            "BYPASS_HIGH_CONFIDENCE_R",
            "High-confidence bypass",
            "rwPFS",
            "Primary adjusted"
        ),

        cox_result(
            cox_os_crude,
            "BYPASS_HIGH_CONFIDENCE_R",
            "High-confidence bypass",
            "OS",
            "Crude"
        ),

        cox_result(
            cox_os_adj,
            "BYPASS_HIGH_CONFIDENCE_R",
            "High-confidence bypass",
            "OS",
            "Primary adjusted"
        )
    )
)


# ============================================================
# 10. Exposure-definition sensitivity models
# ============================================================

fit_exposure_sensitivity <- function(
    exposure_col,
    exposure_label
) {

    formula_rwpfs <- as.formula(
        paste0(
            "Surv(RWPFS_DAYS_LOCKED_R, RWPFS_EVENT_LOCKED_R) ~ ",
            exposure_col,
            " + GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
            "LOG_STAGE4_TO_ANTI_EGFR_R + PRIOR_HEAVY_F_R"
        )
    )

    formula_os <- as.formula(
        paste0(
            "Surv(OS_DAYS_LOCKED_R, OS_EVENT_LOCKED_R) ~ ",
            exposure_col,
            " + GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
            "LOG_STAGE4_TO_ANTI_EGFR_R + PRIOR_HEAVY_F_R"
        )
    )

    fit_rwpfs <- coxph(
        formula_rwpfs,
        data = dt,
        ties = "efron"
    )

    fit_os <- coxph(
        formula_os,
        data = dt,
        ties = "efron"
    )

    rbindlist(
        list(
            cox_result(
                fit_rwpfs,
                exposure_col,
                exposure_label,
                "rwPFS",
                "Adjusted exposure sensitivity"
            ),

            cox_result(
                fit_os,
                exposure_col,
                exposure_label,
                "OS",
                "Adjusted exposure sensitivity"
            )
        )
    )
}

exposure_sensitivity <- rbindlist(
    list(
        fit_exposure_sensitivity(
            "BYPASS_PRESSING_STRICT_R",
            "PRESSING-like strict"
        ),

        fit_exposure_sensitivity(
            "BYPASS_HIGH_CONFIDENCE_NO_D67E_R",
            "High-confidence excluding MAP2K1 D67E"
        ),

        fit_exposure_sensitivity(
            "BYPASS_ORIGINAL_BROAD_R",
            "Original broad composite"
        )
    )
)


# ============================================================
# 11. Cohort-restriction sensitivity
# ============================================================

fit_restricted <- function(
    data_sub,
    label
) {

    if (
        data_sub[
            BYPASS_HIGH_CONFIDENCE_R == 1L,
            .N
        ] < 10L
    ) {

        return(
            data.table(
                OUTCOME = c(
                    "rwPFS",
                    "OS"
                ),
                ANALYSIS =
                    paste0(
                        "Restricted: ",
                        label
                    ),
                EXPOSURE =
                    "High-confidence bypass",
                N =
                    nrow(
                        data_sub
                    ),
                EVENTS = NA_integer_,
                HR = NA_real_,
                LCL95 = NA_real_,
                UCL95 = NA_real_,
                P = NA_real_
            )
        )
    }

    fit_rwpfs <- coxph(
        Surv(
            RWPFS_DAYS_LOCKED_R,
            RWPFS_EVENT_LOCKED_R
        ) ~
            BYPASS_HIGH_CONFIDENCE_R +
            GENDER_F_R +
            SUBSITE_F_R +
            ANTI_EGFR_AGENT_F_R +
            LOG_STAGE4_TO_ANTI_EGFR_R +
            PRIOR_HEAVY_F_R,
        data = data_sub,
        ties = "efron"
    )

    fit_os <- coxph(
        Surv(
            OS_DAYS_LOCKED_R,
            OS_EVENT_LOCKED_R
        ) ~
            BYPASS_HIGH_CONFIDENCE_R +
            GENDER_F_R +
            SUBSITE_F_R +
            ANTI_EGFR_AGENT_F_R +
            LOG_STAGE4_TO_ANTI_EGFR_R +
            PRIOR_HEAVY_F_R,
        data = data_sub,
        ties = "efron"
    )

    rbindlist(
        list(
            cox_result(
                fit_rwpfs,
                "BYPASS_HIGH_CONFIDENCE_R",
                "High-confidence bypass",
                "rwPFS",
                paste0(
                    "Restricted: ",
                    label
                )
            ),

            cox_result(
                fit_os,
                "BYPASS_HIGH_CONFIDENCE_R",
                "High-confidence bypass",
                "OS",
                paste0(
                    "Restricted: ",
                    label
                )
            )
        )
    )
}

restriction_sensitivity <- rbindlist(
    list(
        fit_restricted(
            dt[
                BRAF_CLASS3_R == 0L
            ],
            "exclude BRAF class III"
        ),

        fit_restricted(
            dt[
                SPECIMEN_WITHIN_730D_R == 1L
            ],
            "specimen <=730 days"
        ),

        fit_restricted(
            dt[
                SIDE_CONFLICT_R == 0L
            ],
            "exclude sidedness conflict"
        ),

        fit_restricted(
            dt[
                BRAF_CLASS3_R == 0L &
                SPECIMEN_WITHIN_730D_R == 1L
            ],
            "exclude BRAF class III + specimen <=730 days"
        )
    )
)


# ============================================================
# 12. Additional covariate sensitivity
#     Avoid overloading the primary model.
# ============================================================

# A. Add assay/specimen context
cox_rwpfs_assay <- coxph(
    Surv(
        RWPFS_DAYS_LOCKED_R,
        RWPFS_EVENT_LOCKED_R
    ) ~
        BYPASS_HIGH_CONFIDENCE_R +
        GENDER_F_R +
        SUBSITE_F_R +
        ANTI_EGFR_AGENT_F_R +
        LOG_STAGE4_TO_ANTI_EGFR_R +
        PRIOR_HEAVY_F_R +
        SAMPLE_TYPE +
        PANEL_OLD_R +
        scale(
            SPECIMEN_TO_ANTI_EGFR_DAYS_R
        ),
    data = dt,
    ties = "efron"
)

assay_sensitivity <- cox_result(
    cox_rwpfs_assay,
    "BYPASS_HIGH_CONFIDENCE_R",
    "High-confidence bypass",
    "rwPFS",
    "Adjusted + sample type/panel/specimen age"
)

# B. Add metastatic-site burden (near-complete)
dt_met <- dt[
    !is.na(
        MET_SITE_COUNT_90D_R
    )
]

cox_rwpfs_met <- coxph(
    Surv(
        RWPFS_DAYS_LOCKED_R,
        RWPFS_EVENT_LOCKED_R
    ) ~
        BYPASS_HIGH_CONFIDENCE_R +
        GENDER_F_R +
        SUBSITE_F_R +
        ANTI_EGFR_AGENT_F_R +
        LOG_STAGE4_TO_ANTI_EGFR_R +
        PRIOR_HEAVY_F_R +
        MET_SITE_COUNT_90D_R,
    data = dt_met,
    ties = "efron"
)

met_sensitivity <- cox_result(
    cox_rwpfs_met,
    "BYPASS_HIGH_CONFIDENCE_R",
    "High-confidence bypass",
    "rwPFS",
    "Adjusted + recent metastatic-site burden"
)


# ============================================================
# 13. KM summaries
# ============================================================

km_rwpfs_median_0 <- km_median(
    km_rwpfs,
    0
)

km_rwpfs_median_1 <- km_median(
    km_rwpfs,
    1
)

km_os_median_0 <- km_median(
    km_os,
    0
)

km_os_median_1 <- km_median(
    km_os,
    1
)

rwpfs_180 <- surv_at_time(
    km_rwpfs,
    180
)

rwpfs_365 <- surv_at_time(
    km_rwpfs,
    365
)

os_365 <- surv_at_time(
    km_os,
    365
)

os_730 <- surv_at_time(
    km_os,
    730
)


# ============================================================
# 14. PH diagnostics
# ============================================================

ph_rwpfs_table <- as.data.table(
    as.data.frame(
        ph_rwpfs$table
    ),
    keep.rownames = "TERM"
)

ph_os_table <- as.data.table(
    as.data.frame(
        ph_os$table
    ),
    keep.rownames = "TERM"
)


# ============================================================
# 15. Save locked dataset and model objects
# ============================================================

locked_data_file <- file.path(
    intermediate_dir,
    "A7_B1_locked_analysis.rds"
)

models_file <- file.path(
    results_dir,
    "B1_07_primary_models.rds"
)

saveRDS(
    dt,
    locked_data_file
)

saveRDS(
    list(
        primary_results =
            primary_results,

        exposure_sensitivity =
            exposure_sensitivity,

        restriction_sensitivity =
            restriction_sensitivity,

        assay_sensitivity =
            assay_sensitivity,

        met_sensitivity =
            met_sensitivity,

        km_rwpfs =
            km_rwpfs,

        km_os =
            km_os,

        cox_rwpfs_crude =
            cox_rwpfs_crude,

        cox_rwpfs_adj =
            cox_rwpfs_adj,

        cox_os_crude =
            cox_os_crude,

        cox_os_adj =
            cox_os_adj,

        ph_rwpfs =
            ph_rwpfs,

        ph_os =
            ph_os
    ),
    models_file
)


# ============================================================
# 16. Main text result
# ============================================================

audit_file <- file.path(
    audit_dir,
    "B1_07_first_unblinded_results.txt"
)

lines <- c(
    "B1-07 FIRST FORMAL UNBLINDED RESULTS",
    "===================================",
    "",
    "1. Locked analysis design",
    "-------------------------",
    "Population: left-sided, RAS-WT, BRAF V600E-negative, MSS mCRC",
    "Index: first anti-EGFR exposure",
    "Genomic ascertainment: sequencing completed on/before index",
    "Primary exposure: frozen high-confidence EGFR-bypass panel",
    "Primary outcome: rwPFS; progression must occur STRICTLY AFTER index",
    "",
    paste0(
        "N = ",
        nrow(dt)
    ),
    paste0(
        "High-confidence bypass+ = ",
        dt[
            BYPASS_HIGH_CONFIDENCE_R == 1L,
            .N
        ]
    ),
    paste0(
        "High-confidence bypass- = ",
        dt[
            BYPASS_HIGH_CONFIDENCE_R == 0L,
            .N
        ]
    ),
    paste0(
        "rwPFS events = ",
        dt[
            RWPFS_EVENT_LOCKED_R == 1L,
            .N
        ]
    ),
    paste0(
        "OS events = ",
        dt[
            OS_EVENT_LOCKED_R == 1L,
            .N
        ]
    ),
    "",
    "Frozen subsite:",
    capture.output(
        print(
            subsite_counts
        )
    ),
    paste0(
        "Structured-vs-diagnosis sidedness conflicts: ",
        n_side_conflict
    ),
    "",
    "2. PRIMARY RESULTS",
    "------------------",
    capture.output(
        print(
            primary_results
        )
    ),
    "",
    paste0(
        "rwPFS log-rank p = ",
        fmt_p(
            logrank_rwpfs_p
        )
    ),
    paste0(
        "OS log-rank p = ",
        fmt_p(
            logrank_os_p
        )
    ),
    "",
    "KM median rwPFS, days:",
    paste0(
        "  bypass-: ",
        fmt_num(
            km_rwpfs_median_0,
            1
        )
    ),
    paste0(
        "  bypass+: ",
        fmt_num(
            km_rwpfs_median_1,
            1
        )
    ),
    "",
    "KM median OS, days:",
    paste0(
        "  bypass-: ",
        fmt_num(
            km_os_median_0,
            1
        )
    ),
    paste0(
        "  bypass+: ",
        fmt_num(
            km_os_median_1,
            1
        )
    ),
    "",
    "rwPFS survival at 180 days:",
    capture.output(
        print(
            rwpfs_180
        )
    ),
    "",
    "rwPFS survival at 365 days:",
    capture.output(
        print(
            rwpfs_365
        )
    ),
    "",
    "OS survival at 365 days:",
    capture.output(
        print(
            os_365
        )
    ),
    "",
    "OS survival at 730 days:",
    capture.output(
        print(
            os_730
        )
    ),
    "",
    "3. Exposure-definition sensitivity",
    "----------------------------------",
    capture.output(
        print(
            exposure_sensitivity
        )
    ),
    "",
    "4. Cohort-restriction sensitivity",
    "---------------------------------",
    capture.output(
        print(
            restriction_sensitivity
        )
    ),
    "",
    "5. Additional rwPFS adjustment sensitivity",
    "------------------------------------------",
    capture.output(
        print(
            rbindlist(
                list(
                    assay_sensitivity,
                    met_sensitivity
                )
            )
        )
    ),
    "",
    "6. Proportional-hazards diagnostics",
    "-----------------------------------",
    "rwPFS adjusted Cox:",
    capture.output(
        print(
            ph_rwpfs_table
        )
    ),
    "",
    "OS adjusted Cox:",
    capture.output(
        print(
            ph_os_table
        )
    ),
    "",
    "7. Interpretation guardrails",
    "----------------------------",
    "This is an association study within anti-EGFR-treated patients.",
    "Do NOT call the bypass effect predictive without a valid untreated/alternative-treatment interaction.",
    "Do NOT re-define the genomic exposure after seeing these results.",
    "Component-level genes are too sparse for inferential gene-by-gene Cox models.",
    "",
    "Next step:",
    "Review effect size, CI, PH diagnostics and sensitivity consistency before generating final figures/tables."
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
cat("B1-07 COMPLETE\n")
cat("============================================================\n")
cat(
    "Locked N / bypass+ / rwPFS events / OS events = ",
    nrow(dt),
    " / ",
    dt[
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        .N
    ],
    " / ",
    dt[
        RWPFS_EVENT_LOCKED_R == 1L,
        .N
    ],
    " / ",
    dt[
        OS_EVENT_LOCKED_R == 1L,
        .N
    ],
    "\n",
    sep = ""
)

cat(
    "\nPrimary adjusted rwPFS result:\n"
)

print(
    primary_results[
        OUTCOME == "rwPFS" &
            ANALYSIS == "Primary adjusted"
    ]
)

cat(
    "\nPrimary adjusted OS result:\n"
)

print(
    primary_results[
        OUTCOME == "OS" &
            ANALYSIS == "Primary adjusted"
    ]
)

cat(
    "\nMain result file:\n",
    audit_file,
    "\n"
)

cat(
    "\nLocked analysis RDS:\n",
    locked_data_file,
    "\n"
)
