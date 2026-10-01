# ============================================================
# A7 / B-line
# 14B_B1_14_exploratory_comparator_rwPFS.R
#
# EXPLORATORY ACTIVE-COMPARATOR OUTCOME ANALYSIS
#
# PURPOSE
#   After B1-14 outcome-blind comparator feasibility returned GO,
#   estimate whether the association between the frozen
#   high-confidence EGFR-bypass composite and rwPFS differs between
#   first anti-EGFR exposure and first bevacizumab exposure.
#
# THIS ANALYSIS IS EXPLORATORY.
# It is NOT a primary causal comparative-effectiveness analysis.
# It is NOT external validation.
# It does NOT change the locked anti-EGFR primary analysis.
#
# -----------------------------------------------------------------
# SEQUENTIAL SAFETY / FREEZE LOGIC
# -----------------------------------------------------------------
# 1) Load B1-14 outcome-blind RDS and require gate_decision == "GO".
# 2) Refit the EXACT B1-14 outcome-blind propensity model.
# 3) Verify overlap-weight balance BEFORE reading progression data.
#    If max absolute weighted SMD >0.10, STOP before outcome access.
# 4) Only after balance passes, reconstruct rwPFS from raw progression
#    timeline + OS follow-up using the same time convention as B1-07.
# 5) HARD-VALIDATE the anti-EGFR subset:
#       N = 191
#       rwPFS events = 183
#       OS events = 129
#       locked B1-07 primary adjusted rwPFS HR ~= 2.334372
#    If this validation fails, STOP before interpreting comparator
#    results.
# 6) Run the frozen exploratory comparator models below.
#
# -----------------------------------------------------------------
# PROPENSITY MODEL
# -----------------------------------------------------------------
# EXACTLY inherited from B1-14:
#
# anti-EGFR vs bevacizumab ~
#   sex
#   + frozen left subsite
#   + log(Stage IV -> first biologic)
#   + prior major cytotoxic classes (0/1/2)
#   + active irinotecan
#   + active oxaliplatin
#   + old panel
#   + metastatic specimen
#   + log(specimen -> first biologic interval)
#   + recent metastatic-site count
#
# Overlap weights:
#   anti-EGFR: 1 - PS(anti-EGFR)
#   bevacizumab: PS(anti-EGFR)
#
# -----------------------------------------------------------------
# PRIMARY EXPLORATORY ESTIMAND
# -----------------------------------------------------------------
# rwPFS treatment-by-bypass interaction in the propensity-score
# overlap population.
#
# M0: unweighted crude interaction (descriptive only)
# M1: overlap-weighted interaction (PRIMARY exploratory comparator)
# M2: overlap-weighted doubly adjusted interaction:
#       treatment * bypass
#       + sex
#       + subsite
#       + log(Stage IV -> biologic)
#       + prior major classes 0/1/2
#       + active irinotecan
#       + active oxaliplatin
#
# Robust sandwich variance is used for weighted Cox models.
#
# -----------------------------------------------------------------
# INTERACTION INTERPRETATION — FROZEN BEFORE OUTCOME
# -----------------------------------------------------------------
# Treatment coding:
#   0 = bevacizumab
#   1 = anti-EGFR
#
# Interaction HR = exp(beta_treatment:bypass)
#
#   Interaction HR >1:
#     the anti-EGFR vs bevacizumab hazard ratio is less favorable
#     among bypass+ than among bypass- patients.
#
#   PRIMARY M1 interpretation:
#     - HR >1 and 95% CI excludes 1:
#         "exploratory evidence of treatment-specific heterogeneity"
#     - HR >1 but 95% CI includes 1:
#         "directional but inconclusive heterogeneity"
#     - HR <=1:
#         "no support for the hypothesized treatment-specific pattern"
#
# Regardless of result:
#   DO NOT use "predictive biomarker" as an established conclusion.
#   DO NOT use causal treatment-selection language.
#
# M2 is a doubly adjusted sensitivity analysis and does not replace M1.
#
# -----------------------------------------------------------------
# SUPPORTIVE OUTPUTS
# -----------------------------------------------------------------
# - Four treatment x bypass cell event counts and KM medians.
# - Treatment HR (anti-EGFR vs bevacizumab) within bypass- and bypass+.
# - Bypass HR within bevacizumab and anti-EGFR derived from M1/M2.
# - Separate shared-covariate-adjusted within-arm bypass Cox models.
# - PH diagnostics are reported only; no outcome-driven PH remodeling
#   is performed in this script.
#
# -----------------------------------------------------------------
# OUTPUTS
# -----------------------------------------------------------------
# 04_results/B1_14B_exploratory_comparator_rwPFS.rds
# 06_logs_and_audit/B1_14B_exploratory_comparator_rwPFS.txt
# 07_tables/TableS8_exploratory_comparator_rwPFS.docx
#
# Same filenames are overwritten on rerun.
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

b114_file <- file.path(
    results_dir,
    "B1_14_bevacizumab_comparator_feasibility.rds"
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
    "B1_14B_exploratory_comparator_rwPFS.rds"
)

audit_file <- file.path(
    audit_dir,
    "B1_14B_exploratory_comparator_rwPFS.txt"
)

table_file <- file.path(
    table_dir,
    "TableS8_exploratory_comparator_rwPFS.docx"
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
    b114_file,
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
    "survival",
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
library(survival)
library(officer)
library(flextable)

cat("\n============================================================\n")
cat("B1-14B EXPLORATORY COMPARATOR rwPFS ANALYSIS\n")
cat("============================================================\n\n")


# ============================================================
# 1. Helpers
# ============================================================

fmt_p <- function(p) {

    if (
        length(p) == 0L ||
        is.na(p) ||
        !is.finite(p)
    ) {
        return("NA")
    }

    if (p < 0.001) {
        return("<0.001")
    }

    sprintf("%.3f", p)
}

fmt_hr <- function(
    hr,
    lo,
    hi
) {

    if (
        any(
            !is.finite(
                c(
                    hr,
                    lo,
                    hi
                )
            )
        )
    ) {
        return("NA")
    }

    sprintf(
        "%.2f (%.2f–%.2f)",
        hr,
        lo,
        hi
    )
}

weighted_mean_safe <- function(
    x,
    w
) {

    ok <- is.finite(x) &
        is.finite(w) &
        w >= 0

    if (
        sum(ok) == 0L ||
        sum(
            w[ok]
        ) <= 0
    ) {
        return(NA_real_)
    }

    sum(
        x[ok] *
            w[ok]
    ) /
        sum(
            w[ok]
        )
}

smd_one <- function(
    x,
    trt,
    weights = NULL
) {

    x <- suppressWarnings(
        as.numeric(x)
    )

    trt <- suppressWarnings(
        as.integer(trt)
    )

    if (is.null(weights)) {
        weights <- rep(
            1,
            length(x)
        )
    }

    weights <- suppressWarnings(
        as.numeric(weights)
    )

    ok <- is.finite(x) &
        is.finite(trt) &
        trt %in%
            c(
                0L,
                1L
            ) &
        is.finite(weights) &
        weights >= 0

    x <- x[ok]
    trt <- trt[ok]
    weights <- weights[ok]

    if (
        length(
            unique(trt)
        ) < 2L
    ) {
        return(NA_real_)
    }

    m0 <- weighted_mean_safe(
        x[
            trt == 0L
        ],
        weights[
            trt == 0L
        ]
    )

    m1 <- weighted_mean_safe(
        x[
            trt == 1L
        ],
        weights[
            trt == 1L
        ]
    )

    v0 <- stats::var(
        x[
            trt == 0L
        ],
        na.rm = TRUE
    )

    v1 <- stats::var(
        x[
            trt == 1L
        ],
        na.rm = TRUE
    )

    denom <- sqrt(
        (
            v0 +
                v1
        ) /
            2
    )

    if (
        !is.finite(
            denom
        ) ||
        denom == 0
    ) {

        if (
            isTRUE(
                all.equal(
                    m1,
                    m0
                )
            )
        ) {
            return(0)
        }

        return(
            sign(
                m1 -
                    m0
            ) *
                Inf
        )
    }

    (
        m1 -
            m0
    ) /
        denom
}

extract_coef <- function(
    fit,
    term
) {

    if (is.null(fit)) {

        return(
            data.table(
                HR = NA_real_,
                LCL95 = NA_real_,
                UCL95 = NA_real_,
                P = NA_real_,
                BETA = NA_real_,
                SE = NA_real_
            )
        )
    }

    s <- summary(fit)

    if (
        !term %in%
            rownames(
                s$coefficients
            )
    ) {

        return(
            data.table(
                HR = NA_real_,
                LCL95 = NA_real_,
                UCL95 = NA_real_,
                P = NA_real_,
                BETA = NA_real_,
                SE = NA_real_
            )
        )
    }

    b <- s$coefficients[
        term,
        "coef"
    ]

    se_col <- if (
        "robust se" %in%
            colnames(
                s$coefficients
            )
    ) {
        "robust se"
    } else {
        "se(coef)"
    }

    se <- s$coefficients[
        term,
        se_col
    ]

    p_col <- if (
        "Pr(>|z|)" %in%
            colnames(
                s$coefficients
            )
    ) {
        "Pr(>|z|)"
    } else {
        tail(
            colnames(
                s$coefficients
            ),
            1
        )
    }

    p <- s$coefficients[
        term,
        p_col
    ]

    data.table(
        HR = exp(b),
        LCL95 =
            exp(
                b -
                    1.96 *
                        se
            ),
        UCL95 =
            exp(
                b +
                    1.96 *
                        se
            ),
        P = p,
        BETA = b,
        SE = se
    )
}

linear_combo <- function(
    fit,
    terms,
    multipliers = NULL
) {

    if (
        is.null(fit)
    ) {
        return(
            data.table(
                HR = NA_real_,
                LCL95 = NA_real_,
                UCL95 = NA_real_,
                P = NA_real_,
                BETA = NA_real_,
                SE = NA_real_
            )
        )
    }

    cf <- coef(
        fit
    )

    V <- vcov(
        fit
    )

    if (
        is.null(
            multipliers
        )
    ) {
        multipliers <- rep(
            1,
            length(
                terms
            )
        )
    }

    if (
        length(
            terms
        ) !=
        length(
            multipliers
        )
    ) {
        stop(
            "linear_combo terms/multipliers length mismatch."
        )
    }

    if (
        !all(
            terms %in%
                names(cf)
        )
    ) {
        return(
            data.table(
                HR = NA_real_,
                LCL95 = NA_real_,
                UCL95 = NA_real_,
                P = NA_real_,
                BETA = NA_real_,
                SE = NA_real_
            )
        )
    }

    L <- rep(
        0,
        length(cf)
    )

    names(L) <- names(cf)

    L[
        terms
    ] <- multipliers

    b <- sum(
        L *
            cf
    )

    se <- sqrt(
        as.numeric(
            t(L) %*%
                V %*%
                L
        )
    )

    z <- b /
        se

    p <- 2 *
        pnorm(
            abs(z),
            lower.tail = FALSE
        )

    data.table(
        HR = exp(b),
        LCL95 =
            exp(
                b -
                    1.96 *
                        se
            ),
        UCL95 =
            exp(
                b +
                    1.96 *
                        se
            ),
        P = p,
        BETA = b,
        SE = se
    )
}

km_median_one <- function(
    days,
    event
) {

    ok <- is.finite(days) &
        days > 0 &
        event %in%
            c(
                0L,
                1L
            )

    if (
        sum(ok) ==
            0L
    ) {
        return(NA_real_)
    }

    sf <- survfit(
        Surv(
            days[ok],
            event[ok]
        ) ~ 1
    )

    tab <- summary(
        sf
    )$table

    med <- suppressWarnings(
        as.numeric(
            tab[
                "median"
            ]
        )
    )

    med
}


# ============================================================
# 2. Load B1-14 outcome-blind object; require GO
# ============================================================

b114 <- readRDS(
    b114_file
)

required_slots <- c(
    "patient_audit",
    "ps_formula",
    "gate_decision"
)

missing_slots <- setdiff(
    required_slots,
    names(
        b114
    )
)

if (
    length(
        missing_slots
    ) >
        0L
) {
    stop(
        "B1-14 RDS missing required slots:\n",
        paste(
            missing_slots,
            collapse = "\n"
        )
    )
}

if (
    !identical(
        as.character(
            b114$gate_decision
        ),
        "GO"
    )
) {
    stop(
        paste0(
            "B1-14 gate is not GO. Observed: ",
            as.character(
                b114$gate_decision
            ),
            "\nDo NOT run comparator outcome analysis."
        )
    )
}

d <- as.data.table(
    copy(
        b114$patient_audit
    )
)

if (
    nrow(d) !=
        328L ||
    d[
        FIRST_BIOLOGIC_R ==
            "anti-EGFR",
        .N
    ] !=
        191L ||
    d[
        FIRST_BIOLOGIC_R ==
            "bevacizumab",
        .N
    ] !=
        137L ||
    d[
        FIRST_BIOLOGIC_R ==
            "anti-EGFR" &
        BYPASS_HIGH_CONFIDENCE_R ==
            1L,
        .N
    ] !=
        23L ||
    d[
        FIRST_BIOLOGIC_R ==
            "bevacizumab" &
        BYPASS_HIGH_CONFIDENCE_R ==
            1L,
        .N
    ] !=
        16L
) {
    stop(
        "B1-14 patient-level anchors failed."
    )
}

cat(
    "B1-14 GO object loaded: 191 anti-EGFR / 137 bevacizumab; ",
    "HC+ 23 / 16.\n"
)


# ============================================================
# 3. Refit EXACT B1-14 PS and verify balance BEFORE outcome read
# ============================================================

ps_formula <- b114$ps_formula

ps_vars <- all.vars(
    ps_formula
)

ps_covars <- setdiff(
    ps_vars,
    "FIRST_BIOLOGIC_R"
)

required_ps_vars <- unique(
    c(
        "FIRST_BIOLOGIC_R",
        ps_covars
    )
)

missing_ps <- setdiff(
    required_ps_vars,
    names(d)
)

if (
    length(
        missing_ps
    ) >
        0L
) {
    stop(
        "B1-14 patient audit missing PS variables:\n",
        paste(
            missing_ps,
            collapse = "\n"
        )
    )
}

d[
    ,
    ANTI_EGFR_TREAT_R :=
        as.integer(
            FIRST_BIOLOGIC_R ==
                "anti-EGFR"
        )
]

ps_complete <- complete.cases(
    d[
        ,
        ..ps_covars
    ]
)

d[
    ,
    PS_COMPLETE_B14B_R :=
        ps_complete
]

dps <- d[
    PS_COMPLETE_B14B_R ==
        TRUE
]

if (
    dps[
        FIRST_BIOLOGIC_R ==
            "anti-EGFR",
        .N
    ] !=
        188L ||
    dps[
        FIRST_BIOLOGIC_R ==
            "bevacizumab",
        .N
    ] !=
        131L
) {
    stop(
        paste0(
            "B1-14 PS complete-case anchor changed.\n",
            "Observed anti/BEV = ",
            dps[
                FIRST_BIOLOGIC_R ==
                    "anti-EGFR",
                .N
            ],
            "/",
            dps[
                FIRST_BIOLOGIC_R ==
                    "bevacizumab",
                .N
            ],
            "; expected 188/131."
        )
    )
}

ps_fit <- glm(
    ps_formula,
    data = dps,
    family = binomial()
)

if (
    !isTRUE(
        ps_fit$converged
    ) ||
    any(
        !is.finite(
            coef(
                ps_fit
            )
        )
    )
) {
    stop(
        "Fixed B1-14 propensity model failed."
    )
}

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
    ) ||
    any(
        dps$PS_ANTI_R <=
            0 |
        dps$PS_ANTI_R >=
            1
    )
) {
    stop(
        "Invalid propensity scores."
    )
}

dps[
    ,
    OVERLAP_WEIGHT_R :=
        fifelse(
            ANTI_EGFR_TREAT_R ==
                1L,
            1 -
                PS_ANTI_R,
            PS_ANTI_R
        )
]

# Use the exact PS model matrix, excluding the intercept.
X <- model.matrix(
    ps_fit
)

if (
    "(Intercept)" %in%
        colnames(X)
) {
    X <- X[
        ,
        setdiff(
            colnames(X),
            "(Intercept)"
        ),
        drop = FALSE
    ]
}

balance <- rbindlist(
    lapply(
        colnames(X),
        function(v) {

            x <- X[
                ,
                v
            ]

            data.table(
                Covariate =
                    v,
                SMD_unweighted =
                    smd_one(
                        x,
                        dps$ANTI_EGFR_TREAT_R
                    ),
                SMD_overlap_weighted =
                    smd_one(
                        x,
                        dps$ANTI_EGFR_TREAT_R,
                        dps$OVERLAP_WEIGHT_R
                    )
            )
        }
    )
)

balance[
    ,
    ABS_SMD_unweighted :=
        abs(
            SMD_unweighted
        )
]

balance[
    ,
    ABS_SMD_overlap_weighted :=
        abs(
            SMD_overlap_weighted
        )
]

setorder(
    balance,
    -ABS_SMD_overlap_weighted
)

max_weighted_smd <- max(
    balance$ABS_SMD_overlap_weighted,
    na.rm = TRUE
)

balance_gate_pass <- (
    is.finite(
        max_weighted_smd
    ) &&
    max_weighted_smd <=
        0.10
)

cat("\nPRE-OUTCOME BALANCE CHECK\n")
cat("-------------------------\n")
print(
    balance
)

cat(
    "\nMaximum absolute overlap-weighted SMD = ",
    sprintf(
        "%.4f",
        max_weighted_smd
    ),
    "\n",
    sep = ""
)

if (
    !balance_gate_pass
) {

    preoutcome_lines <- c(
        "B1-14B PRE-OUTCOME BALANCE GATE FAILED",
        "======================================",
        "",
        paste0(
            "Maximum absolute overlap-weighted SMD = ",
            sprintf(
                "%.4f",
                max_weighted_smd
            )
        ),
        "Required <=0.10 before outcome access.",
        "",
        capture.output(
            print(
                balance
            )
        ),
        "",
        "Progression data were NOT read and no outcome model was run."
    )

    writeLines(
        preoutcome_lines,
        audit_file
    )

    stop(
        paste0(
            "Pre-outcome overlap-balance gate failed: max |SMD|=",
            sprintf(
                "%.4f",
                max_weighted_smd
            ),
            " >0.10.\nProgression data were not read."
        )
    )
}

cat(
    "Pre-outcome balance gate: PASS.\n"
)


# ============================================================
# 4. NOW attach follow-up/death information
#    This occurs only AFTER overlap balance passed.
# ============================================================

main <- fread(
    main_file,
    check.names = FALSE,
    showProgress = FALSE
)

required_main_outcome <- c(
    "PATIENT_ID",
    "OS_MONTHS",
    "OS_STATUS"
)

missing_main_outcome <- setdiff(
    required_main_outcome,
    names(
        main
    )
)

if (
    length(
        missing_main_outcome
    ) >
        0L
) {
    stop(
        "Main CSV missing outcome-support fields:\n",
        paste(
            missing_main_outcome,
            collapse = "\n"
        )
    )
}

os_info <- main[
    PATIENT_ID %chin%
        d$PATIENT_ID,
    .(
        PATIENT_ID,
        OS_MONTHS,
        OS_STATUS
    )
]

if (
    uniqueN(
        os_info$PATIENT_ID
    ) !=
        328L
) {
    stop(
        "OS/follow-up map did not return all 328 comparator patients."
    )
}

d <- merge(
    d,
    os_info,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

d[
    ,
    FOLLOWUP_DAY_B14B_R :=
        suppressWarnings(
            as.numeric(
                OS_MONTHS
            )
        ) *
        30.4375
]

d[
    ,
    DECEASED_B14B_R :=
        as.integer(
            OS_STATUS ==
                "1:DECEASED"
        )
]

d[
    ,
    DEATH_DAY_B14B_R :=
        fifelse(
            DECEASED_B14B_R ==
                1L &
            is.finite(
                FOLLOWUP_DAY_B14B_R
            ),
            FOLLOWUP_DAY_B14B_R,
            NA_real_
        )
]


# ============================================================
# 5. Extract raw progression timeline and reconstruct rwPFS
#    EXACT B1-07 rule:
#       first progression STRICTLY AFTER index OR death after index;
#       censor at last follow-up.
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_14B_"
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
        length(hit) ==
            0L
    ) {
        stop(
            "Missing raw MSK file: ",
            filename
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
    PATIENT_ID %chin%
        d$PATIENT_ID
]

progression[
    ,
    START_DATE_NUM_B14B_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

progression_y <- progression[
    PROGRESSION ==
        "Y" &
    is.finite(
        START_DATE_NUM_B14B_R
    )
]

build_outcome <- function(
    pid,
    index_day,
    followup_day,
    death_day
) {

    if (
        !is.finite(
            index_day
        ) ||
        !is.finite(
            followup_day
        ) ||
        followup_day <
            index_day
    ) {
        return(
            data.table(
                PATIENT_ID = pid,
                VALID_FOLLOWUP_B14B_R =
                    0L,
                FIRST_PROGRESSION_DAY_B14B_R =
                    NA_real_,
                DEATH_AFTER_INDEX_B14B_R =
                    0L,
                RWPFS_EVENT_B14B_R =
                    NA_integer_,
                RWPFS_DAYS_B14B_R =
                    NA_real_,
                OS_EVENT_B14B_R =
                    NA_integer_,
                OS_DAYS_B14B_R =
                    NA_real_
            )
        )
    }

    p <- progression_y[
        PATIENT_ID ==
            pid &
        START_DATE_NUM_B14B_R >
            index_day &
        START_DATE_NUM_B14B_R <=
            followup_day
    ]

    prog_day <- if (
        nrow(p) >
            0L
    ) {
        min(
            p$START_DATE_NUM_B14B_R
        )
    } else {
        NA_real_
    }

    death_after_index <- as.integer(
        is.finite(
            death_day
        ) &&
        death_day >
            index_day &&
        death_day <=
            followup_day
    )

    death_event_day <- if (
        death_after_index ==
            1L
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
        ) >
            0L
    )

    rwpfs_end_day <- if (
        rwpfs_event ==
            1L
    ) {
        min(
            event_days
        )
    } else {
        followup_day
    }

    os_event <- death_after_index

    os_end_day <- if (
        os_event ==
            1L
    ) {
        death_day
    } else {
        followup_day
    }

    data.table(
        PATIENT_ID = pid,
        VALID_FOLLOWUP_B14B_R =
            1L,
        FIRST_PROGRESSION_DAY_B14B_R =
            prog_day,
        DEATH_AFTER_INDEX_B14B_R =
            death_after_index,
        RWPFS_EVENT_B14B_R =
            rwpfs_event,
        RWPFS_DAYS_B14B_R =
            rwpfs_end_day -
            index_day,
        OS_EVENT_B14B_R =
            os_event,
        OS_DAYS_B14B_R =
            os_end_day -
            index_day
    )
}

outcome <- rbindlist(
    Map(
        build_outcome,
        d$PATIENT_ID,
        d$FIRST_BIOLOGIC_DAY_R,
        d$FOLLOWUP_DAY_B14B_R,
        d$DEATH_DAY_B14B_R
    ),
    fill = TRUE
)

d <- merge(
    d,
    outcome,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

if (
    d[
        VALID_FOLLOWUP_B14B_R !=
            1L,
        .N
    ] >
        0L ||
    d[
        !is.finite(
            RWPFS_DAYS_B14B_R
        ) |
        RWPFS_DAYS_B14B_R <=
            0,
        .N
    ] >
        0L
) {
    stop(
        "Invalid follow-up or non-positive rwPFS durations in comparator cohort."
    )
}


# ============================================================
# 6. HARD anti-EGFR outcome/model validation
#    Must reproduce locked B1-07 before comparator interpretation.
# ============================================================

anti <- d[
    FIRST_BIOLOGIC_R ==
        "anti-EGFR"
]

if (
    nrow(anti) !=
        191L ||
    sum(
        anti$RWPFS_EVENT_B14B_R
    ) !=
        183L ||
    sum(
        anti$OS_EVENT_B14B_R
    ) !=
        129L
) {
    stop(
        paste0(
            "Anti-EGFR outcome anchors failed.\n",
            "Observed N/rwPFS events/OS events = ",
            nrow(anti), "/",
            sum(
                anti$RWPFS_EVENT_B14B_R
            ), "/",
            sum(
                anti$OS_EVENT_B14B_R
            ),
            "\nExpected 191/183/129."
        )
    )
}

anti[
    ,
    ANTI_EGFR_AGENT_F_B14B_R :=
        factor(
            INDEX_BIOLOGIC_AGENT_R,
            levels = c(
                "CETUXIMAB",
                "PANITUMUMAB"
            )
        )
]

anti[
    ,
    PRIOR_HEAVY_F_B14B_R :=
        factor(
            PRIOR_HEAVY_R,
            levels = c(
                0,
                1
            )
        )
]

anti[
    ,
    GENDER_F_B14B_R :=
        factor(
            GENDER,
            levels = c(
                "Female",
                "Male"
            )
        )
]

anti[
    ,
    SUBSITE_F_B14B_R :=
        factor(
            SUBSITE_FROZEN_R,
            levels = c(
                "left colon",
                "rectum/rectosigmoid"
            )
        )
]

fit_anti_validation <- coxph(
    Surv(
        RWPFS_DAYS_B14B_R,
        RWPFS_EVENT_B14B_R
    ) ~
        BYPASS_HIGH_CONFIDENCE_R +
        GENDER_F_B14B_R +
        SUBSITE_F_B14B_R +
        ANTI_EGFR_AGENT_F_B14B_R +
        LOG_STAGE4_TO_BIOLOGIC_R +
        PRIOR_HEAVY_F_B14B_R,
    data = anti,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

anti_validation_result <- extract_coef(
    fit_anti_validation,
    "BYPASS_HIGH_CONFIDENCE_R"
)

locked_hr <- 2.334372

if (
    !is.finite(
        anti_validation_result$HR
    ) ||
    abs(
        anti_validation_result$HR -
            locked_hr
    ) >
        0.02
) {
    stop(
        paste0(
            "Anti-EGFR locked B1-07 HR validation failed.\n",
            "Observed HR=",
            sprintf(
                "%.6f",
                anti_validation_result$HR
            ),
            "; expected approximately ",
            locked_hr,
            ".\nDo not interpret comparator results."
        )
    )
}

cat(
    "\nHARD anti-EGFR B1-07 validation: PASS\n"
)

cat(
    "N=191; rwPFS events=183; OS events=129; validation HR=",
    sprintf(
        "%.4f",
        anti_validation_result$HR
    ),
    "\n",
    sep = ""
)


# ============================================================
# 7. Merge frozen PS / overlap weights back after outcome validation
# ============================================================

weight_map <- dps[
    ,
    .(
        PATIENT_ID,
        PS_ANTI_R,
        OVERLAP_WEIGHT_R
    )
]

d <- merge(
    d,
    weight_map,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

d[
    ,
    PS_COMPLETE_B14B_R :=
        as.integer(
            is.finite(
                PS_ANTI_R
            ) &
            is.finite(
                OVERLAP_WEIGHT_R
            )
        )
]

dw <- d[
    PS_COMPLETE_B14B_R ==
        1L
]

if (
    nrow(dw) !=
        319L ||
    dw[
        FIRST_BIOLOGIC_R ==
            "anti-EGFR",
        .N
    ] !=
        188L ||
    dw[
        FIRST_BIOLOGIC_R ==
            "bevacizumab",
        .N
    ] !=
        131L
) {
    stop(
        "Weighted analysis population anchor failed."
    )
}

dw[
    ,
    ANTI_EGFR_TREAT_R :=
        as.integer(
            FIRST_BIOLOGIC_R ==
                "anti-EGFR"
        )
]

# Recreate shared factors deterministically.
dw[
    ,
    GENDER_F_SHARED_R :=
        factor(
            GENDER,
            levels = c(
                "Female",
                "Male"
            )
        )
]

dw[
    ,
    SUBSITE_F_SHARED_R :=
        factor(
            SUBSITE_FROZEN_R,
            levels = c(
                "left colon",
                "rectum/rectosigmoid"
            )
        )
]

dw[
    ,
    PRIOR_MAJOR_N_F_SHARED_R :=
        factor(
            PRIOR_MAJOR_N_R,
            levels = c(
                0,
                1,
                2
            )
        )
]

dw[
    ,
    ACTIVE_IRI_F_SHARED_R :=
        factor(
            ACTIVE_IRI_R,
            levels = c(
                0,
                1
            )
        )
]

dw[
    ,
    ACTIVE_OX_F_SHARED_R :=
        factor(
            ACTIVE_OX_R,
            levels = c(
                0,
                1
            )
        )
]


# ============================================================
# 8. Four-cell descriptive endpoint audit
# ============================================================

four_cell <- dw[
    ,
    .(
        N = .N,
        EVENTS =
            sum(
                RWPFS_EVENT_B14B_R
            ),
        EVENT_PCT =
            100 *
            mean(
                RWPFS_EVENT_B14B_R ==
                    1L
            ),
        KM_MEDIAN_DAYS =
            km_median_one(
                RWPFS_DAYS_B14B_R,
                RWPFS_EVENT_B14B_R
            ),
        MEDIAN_RAW_DAYS =
            median(
                RWPFS_DAYS_B14B_R,
                na.rm = TRUE
            )
    ),
    by = .(
        Treatment =
            FIRST_BIOLOGIC_R,
        Bypass =
            BYPASS_HIGH_CONFIDENCE_R
    )
][
    order(
        Treatment,
        Bypass
    )
]


# ============================================================
# 9. Frozen comparator models
# ============================================================

# M0: unweighted crude interaction, descriptive only.
fit_M0 <- coxph(
    Surv(
        RWPFS_DAYS_B14B_R,
        RWPFS_EVENT_B14B_R
    ) ~
        ANTI_EGFR_TREAT_R *
        BYPASS_HIGH_CONFIDENCE_R,
    data = dw,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

# M1: PRIMARY exploratory comparator model.
fit_M1 <- coxph(
    Surv(
        RWPFS_DAYS_B14B_R,
        RWPFS_EVENT_B14B_R
    ) ~
        ANTI_EGFR_TREAT_R *
        BYPASS_HIGH_CONFIDENCE_R,
    data = dw,
    weights =
        OVERLAP_WEIGHT_R,
    cluster =
        PATIENT_ID,
    robust = TRUE,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

# M2: doubly adjusted overlap-weighted sensitivity.
fit_M2 <- coxph(
    Surv(
        RWPFS_DAYS_B14B_R,
        RWPFS_EVENT_B14B_R
    ) ~
        ANTI_EGFR_TREAT_R *
        BYPASS_HIGH_CONFIDENCE_R +
        GENDER_F_SHARED_R +
        SUBSITE_F_SHARED_R +
        LOG_STAGE4_TO_BIOLOGIC_R +
        PRIOR_MAJOR_N_F_SHARED_R +
        ACTIVE_IRI_F_SHARED_R +
        ACTIVE_OX_F_SHARED_R,
    data = dw,
    weights =
        OVERLAP_WEIGHT_R,
    cluster =
        PATIENT_ID,
    robust = TRUE,
    ties = "efron",
    x = TRUE,
    y = TRUE
)

term_treatment <-
    "ANTI_EGFR_TREAT_R"

term_bypass <-
    "BYPASS_HIGH_CONFIDENCE_R"

term_interaction <-
    "ANTI_EGFR_TREAT_R:BYPASS_HIGH_CONFIDENCE_R"


# ============================================================
# 10. Extract interaction and clinically interpretable contrasts
# ============================================================

extract_model_contrasts <- function(
    fit,
    model_id,
    model_label
) {

    interaction <- extract_coef(
        fit,
        term_interaction
    )

    treatment_bypass_neg <- linear_combo(
        fit,
        term_treatment
    )

    treatment_bypass_pos <- linear_combo(
        fit,
        c(
            term_treatment,
            term_interaction
        )
    )

    bypass_bev <- linear_combo(
        fit,
        term_bypass
    )

    bypass_anti <- linear_combo(
        fit,
        c(
            term_bypass,
            term_interaction
        )
    )

    data.table(
        MODEL_ID =
            model_id,
        MODEL =
            model_label,
        N =
            fit$n,
        EVENTS =
            fit$nevent,

        INTERACTION_HR =
            interaction$HR,
        INTERACTION_LCL95 =
            interaction$LCL95,
        INTERACTION_UCL95 =
            interaction$UCL95,
        INTERACTION_P =
            interaction$P,

        ANTI_VS_BEV_BYPASS_NEG_HR =
            treatment_bypass_neg$HR,
        ANTI_VS_BEV_BYPASS_NEG_LCL95 =
            treatment_bypass_neg$LCL95,
        ANTI_VS_BEV_BYPASS_NEG_UCL95 =
            treatment_bypass_neg$UCL95,
        ANTI_VS_BEV_BYPASS_NEG_P =
            treatment_bypass_neg$P,

        ANTI_VS_BEV_BYPASS_POS_HR =
            treatment_bypass_pos$HR,
        ANTI_VS_BEV_BYPASS_POS_LCL95 =
            treatment_bypass_pos$LCL95,
        ANTI_VS_BEV_BYPASS_POS_UCL95 =
            treatment_bypass_pos$UCL95,
        ANTI_VS_BEV_BYPASS_POS_P =
            treatment_bypass_pos$P,

        BYPASS_HR_IN_BEV =
            bypass_bev$HR,
        BYPASS_HR_IN_BEV_LCL95 =
            bypass_bev$LCL95,
        BYPASS_HR_IN_BEV_UCL95 =
            bypass_bev$UCL95,
        BYPASS_HR_IN_BEV_P =
            bypass_bev$P,

        BYPASS_HR_IN_ANTI =
            bypass_anti$HR,
        BYPASS_HR_IN_ANTI_LCL95 =
            bypass_anti$LCL95,
        BYPASS_HR_IN_ANTI_UCL95 =
            bypass_anti$UCL95,
        BYPASS_HR_IN_ANTI_P =
            bypass_anti$P
    )
}

model_contrasts <- rbindlist(
    list(
        extract_model_contrasts(
            fit_M0,
            "M0",
            "Unweighted crude interaction"
        ),

        extract_model_contrasts(
            fit_M1,
            "M1",
            "Overlap-weighted interaction (primary exploratory)"
        ),

        extract_model_contrasts(
            fit_M2,
            "M2",
            "Overlap-weighted doubly adjusted interaction"
        )
    ),
    fill = TRUE
)


# ============================================================
# 11. Shared-covariate adjusted within-arm bypass models
#     Supportive only; not the primary interaction test.
# ============================================================

fit_within_arm <- function(
    arm
) {

    ds <- dw[
        FIRST_BIOLOGIC_R ==
            arm
    ]

    fit <- coxph(
        Surv(
            RWPFS_DAYS_B14B_R,
            RWPFS_EVENT_B14B_R
        ) ~
            BYPASS_HIGH_CONFIDENCE_R +
            GENDER_F_SHARED_R +
            SUBSITE_F_SHARED_R +
            LOG_STAGE4_TO_BIOLOGIC_R +
            PRIOR_MAJOR_N_F_SHARED_R +
            ACTIVE_IRI_F_SHARED_R +
            ACTIVE_OX_F_SHARED_R,
        data = ds,
        ties = "efron",
        x = TRUE,
        y = TRUE
    )

    r <- extract_coef(
        fit,
        "BYPASS_HIGH_CONFIDENCE_R"
    )

    data.table(
        Treatment =
            arm,
        N =
            fit$n,
        EVENTS =
            fit$nevent,
        BYPASS_POS =
            ds[
                BYPASS_HIGH_CONFIDENCE_R ==
                    1L,
                .N
            ],
        HR =
            r$HR,
        LCL95 =
            r$LCL95,
        UCL95 =
            r$UCL95,
        P =
            r$P
    )
}

within_arm <- rbindlist(
    list(
        fit_within_arm(
            "anti-EGFR"
        ),
        fit_within_arm(
            "bevacizumab"
        )
    )
)


# ============================================================
# 12. PH diagnostics — REPORT ONLY, no remodeling
# ============================================================

ph_M1 <- tryCatch(
    cox.zph(
        fit_M1,
        transform = "km",
        terms = TRUE,
        singledf = TRUE
    ),
    error = function(e) NULL
)

ph_M2 <- tryCatch(
    cox.zph(
        fit_M2,
        transform = "km",
        terms = TRUE,
        singledf = TRUE
    ),
    error = function(e) NULL
)

ph_table <- function(
    z,
    model_id
) {

    if (is.null(z)) {
        return(
            data.table(
                MODEL_ID =
                    model_id,
                TERM =
                    NA_character_,
                CHISQ =
                    NA_real_,
                DF =
                    NA_real_,
                P =
                    NA_real_
            )
        )
    }

    zz <- as.data.table(
        as.data.frame(
            z$table
        ),
        keep.rownames =
            "TERM"
    )

    setnames(
        zz,
        names(zz),
        toupper(
            names(zz)
        )
    )

    if (
        !"DF" %in%
            names(zz)
    ) {
        zz[
            ,
            DF :=
                NA_real_
        ]
    }

    zz[
        ,
        MODEL_ID :=
            model_id
    ]

    setcolorder(
        zz,
        c(
            "MODEL_ID",
            "TERM",
            "CHISQ",
            "DF",
            "P"
        )
    )

    zz
}

ph_results <- rbindlist(
    list(
        ph_table(
            ph_M1,
            "M1"
        ),
        ph_table(
            ph_M2,
            "M2"
        )
    ),
    fill = TRUE
)


# ============================================================
# 13. Frozen interpretation gate based on PRIMARY M1 only
# ============================================================

m1 <- model_contrasts[
    MODEL_ID ==
        "M1"
]

if (
    nrow(m1) !=
        1L ||
    !is.finite(
        m1$INTERACTION_HR
    )
) {
    interaction_gate <-
        "NOT INTERPRETABLE"

    interaction_wording <- paste0(
        "The primary exploratory overlap-weighted interaction ",
        "could not be reliably estimated."
    )

} else if (
    m1$INTERACTION_HR >
        1 &&
    m1$INTERACTION_LCL95 >
        1
) {

    interaction_gate <-
        "EXPLORATORY HETEROGENEITY SUPPORT"

    interaction_wording <- paste0(
        "The overlap-weighted rwPFS analysis provides exploratory ",
        "evidence of treatment-specific heterogeneity. This remains ",
        "noncausal and requires external replication."
    )

} else if (
    m1$INTERACTION_HR >
        1
) {

    interaction_gate <-
        "DIRECTIONAL BUT INCONCLUSIVE"

    interaction_wording <- paste0(
        "The interaction was directionally consistent with a less ",
        "favorable anti-EGFR association among bypass-positive patients, ",
        "but the confidence interval included the null."
    )

} else {

    interaction_gate <-
        "NO SUPPORT FOR HYPOTHESIZED SPECIFICITY"

    interaction_wording <- paste0(
        "The primary exploratory interaction did not support the ",
        "hypothesized treatment-specific pattern."
    )
}

m2 <- model_contrasts[
    MODEL_ID ==
        "M2"
]

m2_direction_consistent <- (
    nrow(m2) ==
        1L &&
    is.finite(
        m2$INTERACTION_HR
    ) &&
    is.finite(
        m1$INTERACTION_HR
    ) &&
    sign(
        log(
            m2$INTERACTION_HR
        )
    ) ==
    sign(
        log(
            m1$INTERACTION_HR
        )
    )
)


# ============================================================
# 14. Save RDS
# ============================================================

saveRDS(
    list(
        analysis_status =
            "EXPLORATORY ACTIVE-COMPARATOR rwPFS ANALYSIS",

        frozen_interpretation_rules = c(
            "M1 interaction HR >1 and 95% CI excludes 1: exploratory evidence of treatment-specific heterogeneity.",
            "M1 interaction HR >1 with CI including 1: directional but inconclusive.",
            "M1 interaction HR <=1: no support for hypothesized treatment-specific pattern.",
            "No causal or established predictive-biomarker language."
        ),

        preoutcome_balance =
            balance,

        max_overlap_weighted_abs_smd =
            max_weighted_smd,

        anti_locked_validation =
            list(
                N =
                    nrow(anti),
                rwPFS_events =
                    sum(
                        anti$RWPFS_EVENT_B14B_R
                    ),
                OS_events =
                    sum(
                        anti$OS_EVENT_B14B_R
                    ),
                validation_model =
                    fit_anti_validation,
                validation_result =
                    anti_validation_result
            ),

        analysis_data =
            dw,

        four_cell_endpoint_summary =
            four_cell,

        model_contrasts =
            model_contrasts,

        within_arm_shared_adjusted =
            within_arm,

        models =
            list(
                M0_unweighted =
                    fit_M0,
                M1_overlap_weighted =
                    fit_M1,
                M2_overlap_weighted_doubly_adjusted =
                    fit_M2
            ),

        ph_diagnostics =
            ph_results,

        interaction_gate =
            interaction_gate,

        interaction_wording =
            interaction_wording,

        M2_direction_consistent =
            m2_direction_consistent
    ),
    result_file
)


# ============================================================
# 15. Audit TXT
# ============================================================

audit_lines <- c(
    "B1-14B EXPLORATORY ACTIVE-COMPARATOR rwPFS ANALYSIS",
    "===================================================",
    "",
    "STATUS",
    "------",
    "Exploratory comparator specificity analysis only.",
    "Not a primary causal comparative-effectiveness analysis.",
    "Not external validation.",
    "Frozen high-confidence bypass definition unchanged.",
    "",
    "SEQUENTIAL VALIDATION",
    "---------------------",
    "B1-14 gate decision: GO",
    paste0(
        "PS complete-case cohort: ",
        nrow(dw),
        " (anti-EGFR ",
        dw[
            FIRST_BIOLOGIC_R ==
                "anti-EGFR",
            .N
        ],
        "; bevacizumab ",
        dw[
            FIRST_BIOLOGIC_R ==
                "bevacizumab",
            .N
        ],
        ")"
    ),
    paste0(
        "Maximum absolute overlap-weighted SMD before outcome access: ",
        sprintf(
            "%.4f",
            max_weighted_smd
        ),
        " [required <=0.10; PASS]"
    ),
    "",
    "PRE-OUTCOME BALANCE TABLE",
    "-------------------------",
    capture.output(
        print(
            balance
        )
    ),
    "",
    "HARD ANTI-EGFR B1-07 VALIDATION",
    "-------------------------------",
    paste0(
        "N = ",
        nrow(anti)
    ),
    paste0(
        "rwPFS events = ",
        sum(
            anti$RWPFS_EVENT_B14B_R
        ),
        " [expected 183]"
    ),
    paste0(
        "OS events = ",
        sum(
            anti$OS_EVENT_B14B_R
        ),
        " [expected 129]"
    ),
    paste0(
        "Locked-model reproduced bypass HR = ",
        sprintf(
            "%.6f",
            anti_validation_result$HR
        ),
        " [B1-07 expected 2.334372]"
    ),
    "",
    "FOUR TREATMENT x BYPASS CELLS",
    "-----------------------------",
    capture.output(
        print(
            four_cell
        )
    ),
    "",
    "FROZEN COMPARATOR MODEL RESULTS",
    "-------------------------------",
    capture.output(
        print(
            model_contrasts
        )
    ),
    "",
    "SUPPORTIVE SHARED-COVARIATE WITHIN-ARM BYPASS MODELS",
    "----------------------------------------------------",
    capture.output(
        print(
            within_arm
        )
    ),
    "",
    "PH DIAGNOSTICS — REPORT ONLY",
    "----------------------------",
    capture.output(
        print(
            ph_results
        )
    ),
    "",
    "FROZEN INTERPRETATION GATE",
    "--------------------------",
    paste0(
        "PRIMARY M1 DECISION: ",
        interaction_gate
    ),
    interaction_wording,
    paste0(
        "M2 interaction direction consistent with M1: ",
        m2_direction_consistent
    ),
    "",
    "INTERPRETATION GUARDRAILS",
    "-------------------------",
    "- M1 is the primary exploratory comparator interaction.",
    "- M2 is a doubly adjusted sensitivity and does not replace M1.",
    "- Interaction HR >1 means anti-EGFR vs bevacizumab is relatively less favorable in bypass+ than bypass-.",
    "- Even a statistically significant interaction is exploratory evidence, not proof of treatment predictiveness.",
    "- Do not use causal treatment-selection language.",
    "- Do not call bevacizumab analysis independent validation.",
    "",
    paste0(
        "Saved RDS: ",
        result_file
    ),
    paste0(
        "Saved Table S8: ",
        table_file
    )
)

writeLines(
    audit_lines,
    audit_file
)


# ============================================================
# 16. Word Table S8
# ============================================================

four_cell_word <- four_cell[
    ,
    .(
        Treatment,
        `Bypass status` =
            fifelse(
                Bypass ==
                    1L,
                "High-confidence bypass+",
                "High-confidence bypass-"
            ),
        N,
        Events =
            EVENTS,
        `Event, %` =
            sprintf(
                "%.1f",
                EVENT_PCT
            ),
        `KM median rwPFS, days` =
            ifelse(
                is.finite(
                    KM_MEDIAN_DAYS
                ),
                sprintf(
                    "%.1f",
                    KM_MEDIAN_DAYS
                ),
                "Not reached"
            )
    )
]

model_word <- model_contrasts[
    ,
    .(
        Model =
            MODEL,
        N,
        Events =
            EVENTS,
        `Interaction HR (95% CI)` =
            mapply(
                fmt_hr,
                INTERACTION_HR,
                INTERACTION_LCL95,
                INTERACTION_UCL95
            ),
        `Interaction P` =
            vapply(
                INTERACTION_P,
                fmt_p,
                character(1)
            ),
        `Anti-EGFR vs BEV HR, bypass-` =
            mapply(
                fmt_hr,
                ANTI_VS_BEV_BYPASS_NEG_HR,
                ANTI_VS_BEV_BYPASS_NEG_LCL95,
                ANTI_VS_BEV_BYPASS_NEG_UCL95
            ),
        `Anti-EGFR vs BEV HR, bypass+` =
            mapply(
                fmt_hr,
                ANTI_VS_BEV_BYPASS_POS_HR,
                ANTI_VS_BEV_BYPASS_POS_LCL95,
                ANTI_VS_BEV_BYPASS_POS_UCL95
            ),
        `Bypass HR within BEV` =
            mapply(
                fmt_hr,
                BYPASS_HR_IN_BEV,
                BYPASS_HR_IN_BEV_LCL95,
                BYPASS_HR_IN_BEV_UCL95
            ),
        `Bypass HR within anti-EGFR` =
            mapply(
                fmt_hr,
                BYPASS_HR_IN_ANTI,
                BYPASS_HR_IN_ANTI_LCL95,
                BYPASS_HR_IN_ANTI_UCL95
            )
    )
]

within_word <- within_arm[
    ,
    .(
        Treatment,
        N,
        Events =
            EVENTS,
        `Bypass+, n` =
            BYPASS_POS,
        `Adjusted bypass HR (95% CI)` =
            mapply(
                fmt_hr,
                HR,
                LCL95,
                UCL95
            ),
        `P value` =
            vapply(
                P,
                fmt_p,
                character(1)
            )
    )
]

balance_word <- balance[
    order(
        -ABS_SMD_unweighted
    ),
    .(
        Covariate,
        `Unweighted SMD` =
            sprintf(
                "%.3f",
                SMD_unweighted
            ),
        `Overlap-weighted SMD` =
            sprintf(
                "%.3f",
                SMD_overlap_weighted
            )
    )
]

ft_three_line <- function(
    df,
    font_size = 7.4
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
            orient = "landscape"
        ),
    page_margins =
        page_mar(
            top = 0.40,
            bottom = 0.40,
            left = 0.40,
            right = 0.40
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
            "Supplementary Table S8. Exploratory overlap-weighted bevacizumab comparator analysis of rwPFS",
            fp_text(
                font.family = "Arial",
                font.size = 10.5,
                bold = TRUE
            )
        )
    )
)

doc <- body_add_par(
    doc,
    "A. Treatment-by-bypass endpoint counts"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        four_cell_word,
        7.5
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "B. Treatment-by-bypass interaction models"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        model_word,
        6.8
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "C. Supportive within-arm bypass associations using shared covariates"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        within_word,
        7.4
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "D. Propensity-score balance"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        balance_word,
        7.2
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
                "Notes: This is an exploratory comparator specificity analysis, not a primary causal comparative-effectiveness analysis. ",
                "The primary exploratory model (M1) used propensity-score overlap weights fixed in the outcome-blind B1-14 workflow. ",
                "M2 additionally adjusted for sex, frozen primary subsite, Stage IV-to-biologic interval, prior exposure to 0/1/2 major ",
                "cytotoxic classes, and active irinotecan/oxaliplatin context. Robust sandwich variance was used for weighted Cox models. ",
                "The interaction HR is the ratio of the anti-EGFR-versus-bevacizumab hazard ratio in bypass-positive patients to the ",
                "corresponding treatment HR in bypass-negative patients; values >1 are directionally consistent with a less favorable ",
                "anti-EGFR association among bypass-positive patients. Maximum absolute weighted SMD before outcome access was ",
                sprintf(
                    "%.3f",
                    max_weighted_smd
                ),
                ". Gate classification: ",
                interaction_gate,
                ". No causal or established predictive-biomarker inference is intended."
            ),
            fp_text(
                font.family = "Arial",
                font.size = 7.2
            )
        )
    )
)

print(
    doc,
    target = table_file
)


# ============================================================
# 17. Console summary
# ============================================================

cat("\nFOUR-CELL rwPFS SUMMARY\n")
cat("-----------------------\n")
print(
    four_cell
)

cat("\nCOMPARATOR MODEL CONTRASTS\n")
cat("--------------------------\n")
print(
    model_contrasts
)

cat("\nWITHIN-ARM SHARED-ADJUSTED BYPASS ASSOCIATIONS\n")
cat("---------------------------------------------\n")
print(
    within_arm
)

cat("\nINTERPRETATION\n")
cat("--------------\n")
cat(
    "PRIMARY M1 DECISION: ",
    interaction_gate,
    "\n",
    sep = ""
)
cat(
    interaction_wording,
    "\n"
)
cat(
    "M2 direction consistent with M1: ",
    m2_direction_consistent,
    "\n",
    sep = ""
)

cat("\n============================================================\n")
cat("B1-14B COMPLETE\n")
cat("============================================================\n")

cat(
    "\nAudit:\n",
    audit_file,
    "\n"
)

cat(
    "\nTable S8:\n",
    table_file,
    "\n"
)

cat(
    "\nRDS:\n",
    result_file,
    "\n"
)
