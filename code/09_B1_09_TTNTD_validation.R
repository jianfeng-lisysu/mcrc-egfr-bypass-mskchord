# ============================================================
# A7 / B-line
# 09_B1_09_TTNTD_validation.R
#
# INDEPENDENT TREATMENT-BASED ENDPOINT VALIDATION
#
# PURPOSE
#   Validate the frozen high-confidence EGFR-bypass association
#   using a treatment-based endpoint that is independent of the
#   radiology/NLP progression timeline.
#
# Endpoint:
#   pragmatic TTNTD-like endpoint from first anti-EGFR exposure
#   to:
#     1) initiation of a NEW antineoplastic agent not belonging
#        to the assembled index regimen, OR
#     2) death,
#   whichever occurs first.
#
# Primary regimen-assembly grace window:
#   +/- context around index, with new treatment required >30 d
#   after index.
#
# Sensitivity:
#   new treatment required >60 d after index.
#
# IMPORTANT
#   - Genomic exposure remains frozen from B1-06.
#   - No redefinition based on B1-07/B1-08 results.
#   - This endpoint is treatment-based and does not use the
#     progression timeline.
#   - Does not modify formal *_R.csv files.
#
# INPUTS
#   03_intermediate/A7_B1_locked_analysis.rds
#   01_raw_data/msk_chord_2024.tar
#
# OUTPUTS
#   04_results/B1_09_TTNTD_validation.rds
#   06_logs_and_audit/B1_09_TTNTD_validation.txt
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

locked_file <- file.path(
    project_root,
    "03_intermediate",
    "A7_B1_locked_analysis.rds"
)

tar_file <- file.path(
    project_root,
    "01_raw_data",
    "msk_chord_2024.tar"
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

for (f in c(
    locked_file,
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
cat("B1-09 TREATMENT-BASED ENDPOINT VALIDATION\n")
cat("============================================================\n\n")


# ============================================================
# Helpers
# ============================================================

safe_chr <- function(x) {

    x <- as.character(x)
    x[is.na(x)] <- ""

    x
}

normalize_agent <- function(x) {

    z <- toupper(
        trimws(
            safe_chr(x)
        )
    )

    z <- gsub(
        "[[:space:]]+",
        " ",
        z
    )

    # Common synonym harmonization
    z[z %chin% c(
        "5-FU",
        "5-FLUOROURACIL"
    )] <- "FLUOROURACIL"

    z[z %chin% c(
        "TRIFLURIDINE/TIPIRACIL",
        "TRIFLURIDINE-TIPIRACIL",
        "TAS-102"
    )] <- "TRIFLURIDINE-TIPIRACIL"

    z
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

median_iqr <- function(x) {

    x <- suppressWarnings(
        as.numeric(x)
    )

    x <- x[
        is.finite(x)
    ]

    if (
        length(x) == 0L
    ) {
        return("NA")
    }

    q <- quantile(
        x,
        c(
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

cox_result <- function(
    fit,
    term,
    endpoint,
    analysis
) {

    s <- summary(fit)

    if (
        !term %in%
            rownames(
                s$coefficients
            )
    ) {
        stop(
            "Exposure term missing from Cox model: ",
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
        ENDPOINT = endpoint,
        ANALYSIS = analysis,
        N = fit$n,
        EVENTS = fit$nevent,
        HR = exp(b),
        LCL95 = exp(
            b -
                1.96 *
                se
        ),
        UCL95 = exp(
            b +
                1.96 *
                se
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

    if (
        length(hit) != 1L
    ) {
        return(NA_real_)
    }

    as.numeric(
        tab[
            hit,
            "median"
        ]
    )
}


# ============================================================
# 1. Locked cohort
# ============================================================

dt <- as.data.table(
    readRDS(
        locked_file
    )
)

required_locked <- c(
    "PATIENT_ID",
    "FIRST_BIOLOGIC_DAY_R",
    "FOLLOWUP_DAY_R",
    "DEATH_AFTER_INDEX_R",
    "DEATH_DAY_R",
    "BYPASS_HIGH_CONFIDENCE_R",
    "BYPASS_PRESSING_STRICT_R",
    "BYPASS_ORIGINAL_BROAD_R",
    "GENDER_F_R",
    "SUBSITE_F_R",
    "ANTI_EGFR_AGENT_F_R",
    "LOG_STAGE4_TO_ANTI_EGFR_R",
    "PRIOR_HEAVY_F_R"
)

missing_locked <- setdiff(
    required_locked,
    names(dt)
)

if (
    length(
        missing_locked
    ) > 0L
) {
    stop(
        "Locked analysis dataset missing:\n",
        paste(
            missing_locked,
            collapse = "\n"
        )
    )
}

if (
    nrow(dt) != 191L ||
    dt[
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        .N
    ] != 23L
) {
    stop(
        "Locked cohort anchors failed."
    )
}


# ============================================================
# 2. Extract raw treatment timeline
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_09_"
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
        ) ==
            filename
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

required_tx <- c(
    "PATIENT_ID",
    "START_DATE",
    "STOP_DATE",
    "AGENT"
)

if (
    !all(
        required_tx %in%
            names(treatment)
    )
) {
    stop(
        "Treatment timeline header mismatch."
    )
}

treatment <- treatment[
    PATIENT_ID %chin%
        dt$PATIENT_ID
]

treatment[
    ,
    START_DAY_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

treatment[
    ,
    STOP_DAY_R :=
        suppressWarnings(
            as.numeric(
                STOP_DATE
            )
        )
]

treatment[
    ,
    AGENT_KEY_R :=
        normalize_agent(
            AGENT
        )
]

treatment <- treatment[
    is.finite(
        START_DAY_R
    ) &
        nzchar(
            AGENT_KEY_R
        )
]


# ============================================================
# 3. Patient-level pragmatic TTNTD construction
#
# Index regimen:
#   any agent that
#     A) starts from index-30 through index+30, OR
#     B) started earlier and has a recorded stop >= index.
#
# Next-treatment event:
#   first start of a NEW agent after grace window that was not
#   part of the assembled index regimen.
#
# Death competes only by occurring first; both next treatment
# and death count as endpoint events.
# ============================================================

build_ttntd <- function(
    pid,
    index_day,
    followup_day,
    death_after_index,
    death_day,
    grace_days
) {

    p <- treatment[
        PATIENT_ID == pid
    ]

    if (
        !is.finite(
            index_day
        ) ||
        !is.finite(
            followup_day
        ) ||
        followup_day <=
            index_day
    ) {
        stop(
            "Invalid index/follow-up for patient ",
            pid
        )
    }

    index_rows <- p[
        (
            START_DAY_R >=
                index_day - 30 &
            START_DAY_R <=
                index_day +
                    grace_days
        ) |
        (
            START_DAY_R <
                index_day - 30 &
            is.finite(
                STOP_DAY_R
            ) &
            STOP_DAY_R >=
                index_day
        )
    ]

    index_agents <- unique(
        index_rows$AGENT_KEY_R
    )

    index_agents <- index_agents[
        nzchar(
            index_agents
        )
    ]

    next_rows <- p[
        START_DAY_R >
            index_day +
                grace_days &
        START_DAY_R <=
            followup_day &
        !AGENT_KEY_R %chin%
            index_agents
    ]

    if (
        nrow(
            next_rows
        ) > 0L
    ) {
        setorder(
            next_rows,
            START_DAY_R
        )

        next_day <-
            next_rows$START_DAY_R[1]

        next_agent <-
            next_rows$AGENT_KEY_R[1]

    } else {

        next_day <-
            NA_real_

        next_agent <-
            NA_character_
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

    candidate_days <- c(
        next_day,
        death_event_day
    )

    candidate_days <- candidate_days[
        is.finite(
            candidate_days
        )
    ]

    event <- as.integer(
        length(
            candidate_days
        ) > 0L
    )

    end_day <- if (
        event == 1L
    ) {
        min(
            candidate_days
        )
    } else {
        followup_day
    }

    event_type <- if (
        event == 0L
    ) {

        "censored"

    } else if (
        is.finite(
            next_day
        ) &&
        is.finite(
            death_event_day
        ) &&
        next_day ==
            death_event_day
    ) {

        "next treatment + death same day"

    } else if (
        is.finite(
            next_day
        ) &&
        (
            !is.finite(
                death_event_day
            ) ||
            next_day <
                death_event_day
        )
    ) {

        "next treatment"

    } else {

        "death"
    }

    data.table(
        PATIENT_ID = pid,

        GRACE_DAYS_R =
            grace_days,

        INDEX_REGIMEN_AGENTS_R =
            if (
                length(
                    index_agents
                ) == 0L
            ) {
                NA_character_
            } else {
                paste(
                    sort(
                        index_agents
                    ),
                    collapse = "|"
                )
            },

        INDEX_AGENT_COUNT_R =
            length(
                index_agents
            ),

        NEXT_TREATMENT_DAY_R =
            next_day,

        NEXT_TREATMENT_AGENT_R =
            next_agent,

        TTNTD_EVENT_R =
            event,

        TTNTD_EVENT_TYPE_R =
            event_type,

        TTNTD_DAYS_R =
            end_day -
                index_day
    )
}


ttntd30 <- rbindlist(
    Map(
        function(
            pid,
            index_day,
            followup_day,
            death_after_index,
            death_day
        ) {
            build_ttntd(
                pid,
                index_day,
                followup_day,
                death_after_index,
                death_day,
                30
            )
        },
        dt$PATIENT_ID,
        dt$FIRST_BIOLOGIC_DAY_R,
        dt$FOLLOWUP_DAY_R,
        dt$DEATH_AFTER_INDEX_R,
        dt$DEATH_DAY_R
    ),
    fill = TRUE
)

ttntd60 <- rbindlist(
    Map(
        function(
            pid,
            index_day,
            followup_day,
            death_after_index,
            death_day
        ) {
            build_ttntd(
                pid,
                index_day,
                followup_day,
                death_after_index,
                death_day,
                60
            )
        },
        dt$PATIENT_ID,
        dt$FIRST_BIOLOGIC_DAY_R,
        dt$FOLLOWUP_DAY_R,
        dt$DEATH_AFTER_INDEX_R,
        dt$DEATH_DAY_R
    ),
    fill = TRUE
)

setnames(
    ttntd30,
    c(
        "TTNTD_EVENT_R",
        "TTNTD_EVENT_TYPE_R",
        "TTNTD_DAYS_R"
    ),
    c(
        "TTNTD30_EVENT_R",
        "TTNTD30_EVENT_TYPE_R",
        "TTNTD30_DAYS_R"
    )
)

setnames(
    ttntd60,
    c(
        "TTNTD_EVENT_R",
        "TTNTD_EVENT_TYPE_R",
        "TTNTD_DAYS_R"
    ),
    c(
        "TTNTD60_EVENT_R",
        "TTNTD60_EVENT_TYPE_R",
        "TTNTD60_DAYS_R"
    )
)

ttntd30_keep <- ttntd30[
    ,
    .(
        PATIENT_ID,
        INDEX_REGIMEN_AGENTS_30D_R =
            INDEX_REGIMEN_AGENTS_R,
        INDEX_AGENT_COUNT_30D_R =
            INDEX_AGENT_COUNT_R,
        NEXT_TREATMENT_DAY_30D_R =
            NEXT_TREATMENT_DAY_R,
        NEXT_TREATMENT_AGENT_30D_R =
            NEXT_TREATMENT_AGENT_R,
        TTNTD30_EVENT_R,
        TTNTD30_EVENT_TYPE_R,
        TTNTD30_DAYS_R
    )
]

ttntd60_keep <- ttntd60[
    ,
    .(
        PATIENT_ID,
        INDEX_REGIMEN_AGENTS_60D_R =
            INDEX_REGIMEN_AGENTS_R,
        INDEX_AGENT_COUNT_60D_R =
            INDEX_AGENT_COUNT_R,
        NEXT_TREATMENT_DAY_60D_R =
            NEXT_TREATMENT_DAY_R,
        NEXT_TREATMENT_AGENT_60D_R =
            NEXT_TREATMENT_AGENT_R,
        TTNTD60_EVENT_R,
        TTNTD60_EVENT_TYPE_R,
        TTNTD60_DAYS_R
    )
]

d <- merge(
    dt,
    ttntd30_keep,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

d <- merge(
    d,
    ttntd60_keep,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)


# ============================================================
# 4. Validity audit
# ============================================================

invalid_30 <- d[
    !is.finite(
        TTNTD30_DAYS_R
    ) |
        TTNTD30_DAYS_R <= 0,
    .N
]

invalid_60 <- d[
    !is.finite(
        TTNTD60_DAYS_R
    ) |
        TTNTD60_DAYS_R <= 0,
    .N
]

if (
    invalid_30 != 0L ||
    invalid_60 != 0L
) {
    stop(
        "Non-positive/missing TTNTD durations detected."
    )
}

event_type_30 <- d[
    ,
    .N,
    by = TTNTD30_EVENT_TYPE_R
][order(-N)]

event_type_60 <- d[
    ,
    .N,
    by = TTNTD60_EVENT_TYPE_R
][order(-N)]

next_agent_30 <- d[
    !is.na(
        NEXT_TREATMENT_AGENT_30D_R
    ),
    .N,
    by = NEXT_TREATMENT_AGENT_30D_R
][order(-N)]

next_agent_60 <- d[
    !is.na(
        NEXT_TREATMENT_AGENT_60D_R
    ),
    .N,
    by = NEXT_TREATMENT_AGENT_60D_R
][order(-N)]


# ============================================================
# 5. PRIMARY 30-day TTNTD-like analysis
# ============================================================

surv30 <- Surv(
    d$TTNTD30_DAYS_R,
    d$TTNTD30_EVENT_R
)

km30 <- survfit(
    surv30 ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = d,
    conf.type = "log-log"
)

logrank30 <- survdiff(
    surv30 ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = d
)

logrank30_p <- pchisq(
    logrank30$chisq,
    df = 1,
    lower.tail = FALSE
)

cox30_crude <- coxph(
    surv30 ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = d,
    ties = "efron"
)

cox30_adj <- coxph(
    surv30 ~
        BYPASS_HIGH_CONFIDENCE_R +
        GENDER_F_R +
        SUBSITE_F_R +
        ANTI_EGFR_AGENT_F_R +
        LOG_STAGE4_TO_ANTI_EGFR_R +
        PRIOR_HEAVY_F_R,
    data = d,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

ph30 <- cox.zph(
    cox30_adj,
    transform = "km"
)


# ============================================================
# 6. 60-day grace sensitivity
# ============================================================

surv60 <- Surv(
    d$TTNTD60_DAYS_R,
    d$TTNTD60_EVENT_R
)

km60 <- survfit(
    surv60 ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = d,
    conf.type = "log-log"
)

cox60_adj <- coxph(
    surv60 ~
        BYPASS_HIGH_CONFIDENCE_R +
        GENDER_F_R +
        SUBSITE_F_R +
        ANTI_EGFR_AGENT_F_R +
        LOG_STAGE4_TO_ANTI_EGFR_R +
        PRIOR_HEAVY_F_R,
    data = d,
    ties = "efron"
)


# ============================================================
# 7. Exposure-definition sensitivity on 30-day endpoint
# ============================================================

fit_exposure <- function(
    exposure_col,
    label
) {

    form <- as.formula(
        paste0(
            "Surv(TTNTD30_DAYS_R, TTNTD30_EVENT_R) ~ ",
            exposure_col,
            " + GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
            "LOG_STAGE4_TO_ANTI_EGFR_R + PRIOR_HEAVY_F_R"
        )
    )

    fit <- coxph(
        form,
        data = d,
        ties = "efron"
    )

    out <- cox_result(
        fit,
        exposure_col,
        "TTNTD30",
        label
    )

    out
}

exposure_sensitivity <- rbindlist(
    list(
        fit_exposure(
            "BYPASS_PRESSING_STRICT_R",
            "PRESSING-like strict"
        ),

        fit_exposure(
            "BYPASS_ORIGINAL_BROAD_R",
            "Original broad composite"
        )
    )
)


# ============================================================
# 8. Results tables
# ============================================================

primary_results <- rbindlist(
    list(
        cox_result(
            cox30_crude,
            "BYPASS_HIGH_CONFIDENCE_R",
            "TTNTD30",
            "Crude high-confidence"
        ),

        cox_result(
            cox30_adj,
            "BYPASS_HIGH_CONFIDENCE_R",
            "TTNTD30",
            "Adjusted high-confidence"
        ),

        cox_result(
            cox60_adj,
            "BYPASS_HIGH_CONFIDENCE_R",
            "TTNTD60",
            "Adjusted high-confidence"
        )
    )
)

km_median_30_neg <- km_median(
    km30,
    0
)

km_median_30_pos <- km_median(
    km30,
    1
)

km_median_60_neg <- km_median(
    km60,
    0
)

km_median_60_pos <- km_median(
    km60,
    1
)

event_by_bypass_30 <- d[
    ,
    .(
        N = .N,
        EVENTS =
            sum(
                TTNTD30_EVENT_R,
                na.rm = TRUE
            ),
        NEXT_TREATMENT_EVENTS =
            sum(
                TTNTD30_EVENT_TYPE_R ==
                    "next treatment",
                na.rm = TRUE
            ) +
            sum(
                TTNTD30_EVENT_TYPE_R ==
                    "next treatment + death same day",
                na.rm = TRUE
            ),
        DEATH_FIRST_EVENTS =
            sum(
                TTNTD30_EVENT_TYPE_R ==
                    "death",
                na.rm = TRUE
            ),
        MEDIAN_TTNTD_DAYS_RAW =
            median(
                TTNTD30_DAYS_R,
                na.rm = TRUE
            )
    ),
    by = BYPASS_HIGH_CONFIDENCE_R
][order(
    BYPASS_HIGH_CONFIDENCE_R
)]

ph30_table <- as.data.table(
    as.data.frame(
        ph30$table
    ),
    keep.rownames = "TERM"
)


# ============================================================
# 9. Save
# ============================================================

result_file <- file.path(
    results_dir,
    "B1_09_TTNTD_validation.rds"
)

audit_file <- file.path(
    audit_dir,
    "B1_09_TTNTD_validation.txt"
)

saveRDS(
    list(
        patient_data = d,
        primary_results =
            primary_results,
        exposure_sensitivity =
            exposure_sensitivity,
        km30 = km30,
        km60 = km60,
        ph30 = ph30,
        event_type_30 =
            event_type_30,
        event_type_60 =
            event_type_60,
        next_agent_30 =
            next_agent_30,
        next_agent_60 =
            next_agent_60
    ),
    result_file
)

lines <- c(
    "B1-09 TREATMENT-BASED ENDPOINT VALIDATION",
    "=========================================",
    "",
    "IMPORTANT",
    "---------",
    "This endpoint does NOT use the progression timeline.",
    "The genomic exposure remains frozen from B1-06.",
    "",
    "1. Pragmatic endpoint definition",
    "-------------------------------",
    "Index: first anti-EGFR exposure.",
    "Event: first new antineoplastic agent after regimen-assembly grace period OR death, whichever occurs first.",
    "Primary grace period: 30 days.",
    "Sensitivity grace period: 60 days.",
    "",
    "2. Endpoint validity",
    "--------------------",
    paste0(
        "N = ",
        nrow(d)
    ),
    paste0(
        "Invalid/non-positive TTNTD30 durations = ",
        invalid_30
    ),
    paste0(
        "Invalid/non-positive TTNTD60 durations = ",
        invalid_60
    ),
    "",
    "TTNTD30 event types:",
    capture.output(
        print(
            event_type_30
        )
    ),
    "",
    "TTNTD60 event types:",
    capture.output(
        print(
            event_type_60
        )
    ),
    "",
    "Most common next agents, 30-day definition:",
    capture.output(
        print(
            head(
                next_agent_30,
                20
            )
        )
    ),
    "",
    "3. TTNTD30 events by bypass group",
    "--------------------------------",
    capture.output(
        print(
            event_by_bypass_30
        )
    ),
    "",
    "4. Primary treatment-based validation",
    "-------------------------------------",
    capture.output(
        print(
            primary_results
        )
    ),
    "",
    paste0(
        "TTNTD30 log-rank p = ",
        fmt_p(
            logrank30_p
        )
    ),
    "",
    paste0(
        "KM median TTNTD30 bypass- = ",
        median_iqr(
            km_median_30_neg
        ),
        " days"
    ),
    paste0(
        "KM median TTNTD30 bypass+ = ",
        median_iqr(
            km_median_30_pos
        ),
        " days"
    ),
    "",
    paste0(
        "KM median TTNTD60 bypass- = ",
        median_iqr(
            km_median_60_neg
        ),
        " days"
    ),
    paste0(
        "KM median TTNTD60 bypass+ = ",
        median_iqr(
            km_median_60_pos
        ),
        " days"
    ),
    "",
    "5. Exposure-definition sensitivity on TTNTD30",
    "---------------------------------------------",
    capture.output(
        print(
            exposure_sensitivity
        )
    ),
    "",
    "6. PH diagnostics for adjusted TTNTD30",
    "--------------------------------------",
    capture.output(
        print(
            ph30_table
        )
    ),
    "",
    "Interpretation guardrails",
    "-------------------------",
    "- This is an independent treatment-based endpoint, not external validation.",
    "- A new-agent start is a pragmatic real-world proxy for subsequent treatment change and may not perfectly reconstruct formal line of therapy.",
    "- Concordant rwPFS and TTNTD effects strengthen robustness but do not establish predictive treatment-effect modification.",
    "",
    "Next gate:",
    "If high-confidence bypass remains associated with shorter TTNTD under both 30- and 60-day definitions, the core results are ready for manuscript figure/table production."
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
cat("B1-09 COMPLETE\n")
cat("============================================================\n")

cat("\nPrimary treatment-based validation:\n")
print(
    primary_results
)

cat(
    "\nMain audit:\n",
    audit_file,
    "\n"
)
