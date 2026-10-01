# ============================================================
# A7 / B-line
# 10_B1_10_refined_TTNTD_PH_validation.R
#
# POST-HOC TECHNICAL REFINEMENT OF B1-09
#
# PURPOSE
#   1) Refine the treatment-based endpoint so clearly non-CRC /
#      supportive agents do not trigger "next treatment".
#   2) Re-estimate high-confidence bypass association.
#   3) Address nuisance-covariate PH violations using a
#      stratified Cox model.
#
# IMPORTANT
#   - Frozen genomic exposure is NOT changed.
#   - rwPFS primary result is NOT changed.
#   - This is a technical endpoint-refinement sensitivity.
#
# Refined endpoint:
#   first new CRC-compatible cancer-directed agent after the
#   regimen-assembly grace period OR death, whichever occurs first.
#
# Two treatment filters:
#   A. Conservative CRC-compatible systemic/targeted whitelist
#   B. Broad cancer-directed filter excluding a priori clearly
#      non-CRC/supportive/endocrine agents
#
# Grace windows:
#   30 and 60 days
#
# PH-addressed sensitivity:
#   strata(SUBSITE_F_R, PRIOR_HEAVY_F_R)
#
# INPUTS
#   03_intermediate/A7_B1_locked_analysis.rds
#   01_raw_data/msk_chord_2024.tar
#
# OUTPUTS
#   04_results/B1_10_refined_TTNTD_PH_validation.rds
#   06_logs_and_audit/B1_10_refined_TTNTD_PH_validation.txt
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
cat("B1-10 REFINED TTNTD / PH VALIDATION\n")
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

    z[z %chin% c(
        "5-FU",
        "5-FLUOROURACIL"
    )] <- "FLUOROURACIL"

    z[z %chin% c(
        "TRIFLURIDINE/TIPIRACIL",
        "TRIFLURIDINE-TIPIRACIL",
        "TIPIRACIL-TRIFLURIDINE",
        "TAS-102"
    )] <- "TRIFLURIDINE-TIPIRACIL"

    z
}

cox_extract <- function(
    fit,
    term,
    endpoint,
    model
) {

    s <- summary(fit)

    if (
        !term %in%
            rownames(
                s$coefficients
            )
    ) {
        stop(
            "Exposure term not found: ",
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
        MODEL = model,
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


# ============================================================
# 1. Locked cohort
# ============================================================

dt <- as.data.table(
    readRDS(
        locked_file
    )
)

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
# 2. Extract treatment timeline
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_10_"
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

tx <- fread(
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

tx <- tx[
    PATIENT_ID %chin%
        dt$PATIENT_ID
]

tx[
    ,
    START_DAY_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

tx[
    ,
    STOP_DAY_R :=
        suppressWarnings(
            as.numeric(
                STOP_DATE
            )
        )
]

tx[
    ,
    AGENT_KEY_R :=
        normalize_agent(
            AGENT
        )
]

tx <- tx[
    is.finite(
        START_DAY_R
    ) &
        nzchar(
            AGENT_KEY_R
        )
]


# ============================================================
# 3. A priori treatment classification
# ============================================================

# Conservative whitelist of CRC-compatible systemic / targeted
# therapies plus trial therapy.
crc_compatible_agents <- c(
    "FLUOROURACIL",
    "CAPECITABINE",
    "OXALIPLATIN",
    "IRINOTECAN",
    "FLOXURIDINE",
    "BEVACIZUMAB",
    "CETUXIMAB",
    "PANITUMUMAB",
    "AFLIBERCEPT",
    "ZIV-AFLIBERCEPT",
    "RAMUCIRUMAB",
    "REGORAFENIB",
    "TRIFLURIDINE-TIPIRACIL",
    "FRUQUINTINIB",
    "ENCORAFENIB",
    "TRASTUZUMAB",
    "PERTUZUMAB",
    "TUCATINIB",
    "LAPATINIB",
    "TRASTUZUMAB DERUXTECAN",
    "FAM-TRASTUZUMAB DERUXTECAN",
    "PEMBROLIZUMAB",
    "NIVOLUMAB",
    "IPILIMUMAB",
    "SOTORASIB",
    "ADAGRASIB",
    "INVESTIGATIONAL"
)

# A priori clearly non-CRC/supportive/endocrine agents that should
# not trigger a CRC next-treatment event in the broad filter.
excluded_non_crc_supportive <- c(
    "ZOLEDRONIC ACID",
    "DENOSUMAB",
    "PAMIDRONATE",
    "MEGESTROL",
    "LETROZOLE",
    "ANASTROZOLE",
    "EXEMESTANE",
    "TAMOXIFEN",
    "FULVESTRANT",
    "ISOTRETINOIN",
    "LEUPROLIDE",
    "GOSERELIN",
    "DEGARELIX",
    "ENZALUTAMIDE",
    "ABIRATERONE",
    "BICALUTAMIDE"
)

# Locoregional therapies are not part of the conservative systemic
# whitelist. In the broad cancer-directed filter they remain eligible
# unless explicitly excluded below.
locoregional_agents <- c(
    "YTTRIUM Y-90 THERASPHERES",
    "YTTRIUM Y-90 MICROSPHERES"
)

tx[
    ,
    IS_CRC_CONSERVATIVE_R :=
        AGENT_KEY_R %chin%
            crc_compatible_agents
]

tx[
    ,
    IS_BROAD_CANCER_DIRECTED_R :=
        !AGENT_KEY_R %chin%
            excluded_non_crc_supportive
]


# ============================================================
# 4. Endpoint builder
# ============================================================

build_endpoint <- function(
    pid,
    index_day,
    followup_day,
    death_after_index,
    death_day,
    grace_days,
    filter_col
) {

    p <- tx[
        PATIENT_ID == pid
    ]

    if (
        !is.finite(index_day) ||
        !is.finite(followup_day) ||
        followup_day <= index_day
    ) {
        stop(
            "Invalid follow-up for ",
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

    candidate <- p[
        START_DAY_R >
            index_day +
                grace_days &
        START_DAY_R <=
            followup_day &
        !AGENT_KEY_R %chin%
            index_agents &
        get(filter_col) == TRUE
    ]

    if (
        nrow(candidate) > 0L
    ) {

        setorder(
            candidate,
            START_DAY_R
        )

        next_day <-
            candidate$START_DAY_R[1]

        next_agent <-
            candidate$AGENT_KEY_R[1]

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
        is.finite(next_day) &&
        (
            !is.finite(
                death_event_day
            ) ||
            next_day <
                death_event_day
        )
    ) {
        "next treatment"
    } else if (
        is.finite(next_day) &&
        is.finite(
            death_event_day
        ) &&
        next_day ==
            death_event_day
    ) {
        "next treatment + death same day"
    } else {
        "death"
    }

    data.table(
        PATIENT_ID = pid,
        EVENT_R = event,
        DAYS_R =
            end_day - index_day,
        EVENT_TYPE_R =
            event_type,
        NEXT_AGENT_R =
            next_agent
    )
}

make_endpoint <- function(
    grace_days,
    filter_col,
    prefix
) {

    out <- rbindlist(
        Map(
            function(
                pid,
                index_day,
                followup_day,
                death_after_index,
                death_day
            ) {
                build_endpoint(
                    pid,
                    index_day,
                    followup_day,
                    death_after_index,
                    death_day,
                    grace_days,
                    filter_col
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
        out,
        c(
            "EVENT_R",
            "DAYS_R",
            "EVENT_TYPE_R",
            "NEXT_AGENT_R"
        ),
        paste0(
            prefix,
            c(
                "_EVENT_R",
                "_DAYS_R",
                "_EVENT_TYPE_R",
                "_NEXT_AGENT_R"
            )
        )
    )

    out
}


# ============================================================
# 5. Build four refined endpoints
# ============================================================

ep_cons30 <- make_endpoint(
    30,
    "IS_CRC_CONSERVATIVE_R",
    "CRC30"
)

ep_cons60 <- make_endpoint(
    60,
    "IS_CRC_CONSERVATIVE_R",
    "CRC60"
)

ep_broad30 <- make_endpoint(
    30,
    "IS_BROAD_CANCER_DIRECTED_R",
    "BROAD30"
)

ep_broad60 <- make_endpoint(
    60,
    "IS_BROAD_CANCER_DIRECTED_R",
    "BROAD60"
)

d <- copy(dt)

for (obj in list(
    ep_cons30,
    ep_cons60,
    ep_broad30,
    ep_broad60
)) {
    d <- merge(
        d,
        obj,
        by = "PATIENT_ID",
        all.x = TRUE,
        sort = FALSE
    )
}

endpoint_cols <- c(
    "CRC30_DAYS_R",
    "CRC60_DAYS_R",
    "BROAD30_DAYS_R",
    "BROAD60_DAYS_R"
)

for (v in endpoint_cols) {

    if (
        d[
            !is.finite(
                get(v)
            ) |
            get(v) <= 0,
            .N
        ] != 0L
    ) {
        stop(
            "Invalid refined endpoint durations: ",
            v
        )
    }
}


# ============================================================
# 6. Model helper
# ============================================================

fit_endpoint <- function(
    days_col,
    event_col,
    endpoint_label
) {

    # Standard adjusted Cox, for comparability with B1-09
    f_standard <- as.formula(
        paste0(
            "Surv(",
            days_col,
            ", ",
            event_col,
            ") ~ BYPASS_HIGH_CONFIDENCE_R + ",
            "GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
            "LOG_STAGE4_TO_ANTI_EGFR_R + PRIOR_HEAVY_F_R"
        )
    )

    fit_standard <- coxph(
        f_standard,
        data = d,
        ties = "efron",
        x = TRUE,
        y = TRUE
    )

    # PH-addressed model:
    # nuisance variables with detected B1-09 PH issues become strata.
    f_stratified <- as.formula(
        paste0(
            "Surv(",
            days_col,
            ", ",
            event_col,
            ") ~ BYPASS_HIGH_CONFIDENCE_R + ",
            "GENDER_F_R + ANTI_EGFR_AGENT_F_R + ",
            "LOG_STAGE4_TO_ANTI_EGFR_R + ",
            "strata(SUBSITE_F_R, PRIOR_HEAVY_F_R)"
        )
    )

    fit_stratified <- coxph(
        f_stratified,
        data = d,
        ties = "efron",
        x = TRUE,
        y = TRUE
    )

    ph_standard <- cox.zph(
        fit_standard,
        transform = "km"
    )

    ph_stratified <- cox.zph(
        fit_stratified,
        transform = "km"
    )

    list(
        results = rbindlist(
            list(
                cox_extract(
                    fit_standard,
                    "BYPASS_HIGH_CONFIDENCE_R",
                    endpoint_label,
                    "Standard adjusted"
                ),

                cox_extract(
                    fit_stratified,
                    "BYPASS_HIGH_CONFIDENCE_R",
                    endpoint_label,
                    "PH-addressed stratified"
                )
            )
        ),
        standard =
            fit_standard,
        stratified =
            fit_stratified,
        ph_standard =
            ph_standard,
        ph_stratified =
            ph_stratified
    )
}


# ============================================================
# 7. Fit all endpoints
# ============================================================

m_crc30 <- fit_endpoint(
    "CRC30_DAYS_R",
    "CRC30_EVENT_R",
    "Conservative CRC TTNTD30"
)

m_crc60 <- fit_endpoint(
    "CRC60_DAYS_R",
    "CRC60_EVENT_R",
    "Conservative CRC TTNTD60"
)

m_broad30 <- fit_endpoint(
    "BROAD30_DAYS_R",
    "BROAD30_EVENT_R",
    "Broad cancer-directed TTNTD30"
)

m_broad60 <- fit_endpoint(
    "BROAD60_DAYS_R",
    "BROAD60_EVENT_R",
    "Broad cancer-directed TTNTD60"
)

model_results <- rbindlist(
    list(
        m_crc30$results,
        m_crc60$results,
        m_broad30$results,
        m_broad60$results
    )
)


# ============================================================
# 8. Event audit
# ============================================================

event_summary <- rbindlist(
    list(
        d[
            ,
            .(
                ENDPOINT =
                    "Conservative CRC TTNTD30",
                N = .N,
                EVENTS =
                    sum(
                        CRC30_EVENT_R
                    ),
                BYPASS_POS_EVENTS =
                    sum(
                        CRC30_EVENT_R[
                            BYPASS_HIGH_CONFIDENCE_R == 1L
                        ]
                    )
            )
        ],

        d[
            ,
            .(
                ENDPOINT =
                    "Conservative CRC TTNTD60",
                N = .N,
                EVENTS =
                    sum(
                        CRC60_EVENT_R
                    ),
                BYPASS_POS_EVENTS =
                    sum(
                        CRC60_EVENT_R[
                            BYPASS_HIGH_CONFIDENCE_R == 1L
                        ]
                    )
            )
        ],

        d[
            ,
            .(
                ENDPOINT =
                    "Broad cancer-directed TTNTD30",
                N = .N,
                EVENTS =
                    sum(
                        BROAD30_EVENT_R
                    ),
                BYPASS_POS_EVENTS =
                    sum(
                        BROAD30_EVENT_R[
                            BYPASS_HIGH_CONFIDENCE_R == 1L
                        ]
                    )
            )
        ],

        d[
            ,
            .(
                ENDPOINT =
                    "Broad cancer-directed TTNTD60",
                N = .N,
                EVENTS =
                    sum(
                        BROAD60_EVENT_R
                    ),
                BYPASS_POS_EVENTS =
                    sum(
                        BROAD60_EVENT_R[
                            BYPASS_HIGH_CONFIDENCE_R == 1L
                        ]
                    )
            )
        ]
    ),
    fill = TRUE
)

next_agent_crc30 <- d[
    !is.na(
        CRC30_NEXT_AGENT_R
    ),
    .N,
    by = CRC30_NEXT_AGENT_R
][order(-N)]

next_agent_broad30 <- d[
    !is.na(
        BROAD30_NEXT_AGENT_R
    ),
    .N,
    by = BROAD30_NEXT_AGENT_R
][order(-N)]

ignored_agent_counts <- tx[
    AGENT_KEY_R %chin%
        excluded_non_crc_supportive,
    .N,
    by = AGENT_KEY_R
][order(-N)]

all_agent_classification <- unique(
    tx[
        ,
        .(
            AGENT_KEY_R,
            IS_CRC_CONSERVATIVE_R,
            IS_BROAD_CANCER_DIRECTED_R
        )
    ]
)

setorder(
    all_agent_classification,
    AGENT_KEY_R
)


# ============================================================
# 9. PH diagnostics tables
# ============================================================

ph_crc30_standard <- as.data.table(
    as.data.frame(
        m_crc30$ph_standard$table
    ),
    keep.rownames = "TERM"
)

ph_crc30_stratified <- as.data.table(
    as.data.frame(
        m_crc30$ph_stratified$table
    ),
    keep.rownames = "TERM"
)


# ============================================================
# 10. Save results
# ============================================================

result_file <- file.path(
    results_dir,
    "B1_10_refined_TTNTD_PH_validation.rds"
)

audit_file <- file.path(
    audit_dir,
    "B1_10_refined_TTNTD_PH_validation.txt"
)

saveRDS(
    list(
        patient_data =
            d,
        model_results =
            model_results,
        event_summary =
            event_summary,
        next_agent_crc30 =
            next_agent_crc30,
        next_agent_broad30 =
            next_agent_broad30,
        ignored_agent_counts =
            ignored_agent_counts,
        agent_classification =
            all_agent_classification,
        models = list(
            crc30 = m_crc30,
            crc60 = m_crc60,
            broad30 = m_broad30,
            broad60 = m_broad60
        )
    ),
    result_file
)

lines <- c(
    "B1-10 REFINED TTNTD / PH VALIDATION",
    "===================================",
    "",
    "IMPORTANT",
    "---------",
    "This is a post-hoc technical refinement of the independent treatment-based endpoint.",
    "The frozen high-confidence genomic exposure is unchanged.",
    "",
    "1. Endpoint refinements",
    "-----------------------",
    "Conservative filter: new CRC-compatible systemic/targeted agent OR death.",
    "Broad filter: new cancer-directed agent excluding a priori clearly non-CRC/supportive/endocrine agents OR death.",
    "Both 30-day and 60-day regimen-assembly grace periods are evaluated.",
    "",
    "2. Event counts",
    "---------------",
    capture.output(
        print(
            event_summary
        )
    ),
    "",
    "Top next agents under conservative CRC 30-day definition:",
    capture.output(
        print(
            head(
                next_agent_crc30,
                20
            )
        )
    ),
    "",
    "Top next agents under broad cancer-directed 30-day definition:",
    capture.output(
        print(
            head(
                next_agent_broad30,
                20
            )
        )
    ),
    "",
    "A priori excluded non-CRC/supportive agent occurrences:",
    capture.output(
        print(
            ignored_agent_counts
        )
    ),
    "",
    "3. Refined treatment-based validation",
    "-------------------------------------",
    capture.output(
        print(
            model_results
        )
    ),
    "",
    "4. PH diagnostics: conservative CRC TTNTD30 standard model",
    "---------------------------------------------------------",
    capture.output(
        print(
            ph_crc30_standard
        )
    ),
    "",
    "5. PH diagnostics: conservative CRC TTNTD30 PH-addressed stratified model",
    "-------------------------------------------------------------------",
    capture.output(
        print(
            ph_crc30_stratified
        )
    ),
    "",
    "Interpretation gate",
    "-------------------",
    "The treatment-based validation is considered robust only if:",
    "  - high-confidence bypass HR remains >1 under conservative and broad filters;",
    "  - direction is concordant under 30- and 60-day grace windows;",
    "  - PH-addressed stratified estimates remain materially similar to standard adjusted estimates.",
    "",
    "If these criteria are met, stop adding primary analyses and proceed to manuscript figures/tables."
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
cat("B1-10 COMPLETE\n")
cat("============================================================\n")

print(
    model_results
)

cat(
    "\nAudit:\n",
    audit_file,
    "\n"
)
