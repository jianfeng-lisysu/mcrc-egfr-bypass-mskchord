# ============================================================
# A7 / B-line
# 13C_B1_13c_specimen_timing_diagnostic.R
#
# POST-HOC DIAGNOSTIC ONLY — NO NEW GATE
#
# PURPOSE
#   Diagnose whether the smaller / less precise estimate under the
#   specimen <=730-day restriction reflects a clear specimen-age
#   gradient or simply a smaller, different subset.
#
# IMPORTANT
#   - Frozen high-confidence genomic exposure is unchanged.
#   - Primary rwPFS endpoint is unchanged.
#   - Model A remains the locked primary model.
#   - B1-13 Model E remains the predeclared wording ceiling.
#   - B1-13c DOES NOT create a new wording or significance gate.
#   - All prespecified cutoffs below are reported regardless of result.
#   - No cutoff is selected as "best".
#
# DIAGNOSTICS
#   1) Continuous bypass x specimen-age interaction:
#      specimen age is log1p(days), standardized to 1 SD.
#      The interaction HR is the ratio of the bypass HR per 1-SD older
#      log specimen age.
#
#      Two versions:
#        a) Expanded Model-C architecture.
#        b) Fixed PH-addressed Model-C architecture inherited from
#           B1-13A (same nuisance-term handling; no new PH selection).
#
#   2) Prespecified specimen-age restrictions:
#        <=365, <=540, <=730, <=1095 days, plus unrestricted.
#
#      For each cutoff, report:
#        N, events, bypass+, standard expanded HR/CI,
#        and the same fixed PH-addressed architecture.
#
#   3) Describe patients excluded by the <=730-day restriction,
#      especially bypass+ patients:
#        specimen age, rwPFS, event status.
#
# INTERPRETATION RED LINE
#   - A non-significant interaction does NOT prove absence of effect
#     modification because exposed N is small.
#   - If no clear monotonic pattern is seen, manuscript wording may say:
#       "No clear monotonic specimen-age gradient was observed,
#        although interaction analyses were underpowered."
#   - Do NOT say attenuation was "proved" to be only reduced precision.
#
# INPUT
#   04_results/B1_13_additional_sensitivity_adjustment.rds
#
# OUTPUTS
#   04_results/B1_13c_specimen_timing_diagnostic.rds
#   06_logs_and_audit/B1_13c_specimen_timing_diagnostic.txt
#   07_tables/TableS6_specimen_timing_diagnostic.docx
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

input_rds <- file.path(
    project_root,
    "04_results",
    "B1_13_additional_sensitivity_adjustment.rds"
)

output_rds <- file.path(
    project_root,
    "04_results",
    "B1_13c_specimen_timing_diagnostic.rds"
)

audit_file <- file.path(
    project_root,
    "06_logs_and_audit",
    "B1_13c_specimen_timing_diagnostic.txt"
)

table_file <- file.path(
    project_root,
    "07_tables",
    "TableS6_specimen_timing_diagnostic.docx"
)

if (!file.exists(input_rds)) {
    stop(
        "Missing B1-13 RDS:\n",
        input_rds
    )
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
cat("B1-13c SPECIMEN-TIMING DIAGNOSTIC\n")
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

extract_term <- function(
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
                SE = NA_real_,
                STATUS =
                    "NOT ESTIMABLE"
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
                SE = NA_real_,
                STATUS =
                    "TERM NOT ESTIMABLE"
            )
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
        SE = se,
        STATUS = "OK"
    )
}

fit_cox_safe <- function(
    formula,
    data,
    tt_fun = NULL
) {

    # coxph does not support model=TRUE for formulas containing tt().
    # Keep model=TRUE for ordinary Cox models, but turn it off whenever
    # a time-varying coefficient function is supplied.
    tryCatch(
        suppressWarnings(
            coxph(
                formula,
                data = data,
                ties = "efron",
                na.action = na.omit,
                x = TRUE,
                y = TRUE,
                model = is.null(tt_fun),
                tt = tt_fun
            )
        ),
        error = function(e) {
            message(
                "Cox model failed: ",
                conditionMessage(e)
            )
            NULL
        }
    )
}

model_complete_counts <- function(
    formula,
    data
) {

    mf <- tryCatch(
        model.frame(
            formula,
            data = data,
            na.action = na.omit
        ),
        error = function(e) NULL
    )

    if (is.null(mf)) {

        return(
            data.table(
                N = NA_integer_,
                EVENTS = NA_integer_,
                BYPASS_POS = NA_integer_
            )
        )
    }

    y <- model.response(
        mf
    )

    status <- if (
        inherits(
            y,
            "Surv"
        )
    ) {
        y[
            ,
            ncol(y)
        ]
    } else {
        rep(
            NA_real_,
            nrow(mf)
        )
    }

    data.table(
        N = nrow(mf),
        EVENTS =
            sum(
                status,
                na.rm = TRUE
            ),
        BYPASS_POS =
            if (
                "BYPASS_HIGH_CONFIDENCE_R" %in%
                    names(mf)
            ) {
                sum(
                    mf$BYPASS_HIGH_CONFIDENCE_R == 1L,
                    na.rm = TRUE
                )
            } else {
                NA_integer_
            }
    )
}


# ============================================================
# 2. Read object and hard anchors
# ============================================================

obj <- readRDS(
    input_rds
)

if (
    !"patient_data" %in%
        names(obj)
) {
    stop(
        "patient_data missing from B1-13 RDS."
    )
}

d <- as.data.table(
    obj$patient_data
)

if (
    nrow(d) != 191L ||
    d[
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        .N
    ] != 23L ||
    sum(
        d$RWPFS_EVENT_B13_R
    ) != 183L
) {
    stop(
        "Locked B1-13 anchors failed."
    )
}

required_cols <- c(
    "BYPASS_HIGH_CONFIDENCE_R",
    "RWPFS_DAYS_B13_R",
    "RWPFS_EVENT_B13_R",
    "GENDER_F_R",
    "SUBSITE_F_R",
    "ANTI_EGFR_AGENT_F_R",
    "LOG_STAGE4_TO_ANTI_EGFR_R",
    "PRIOR_MAJOR_N_F_B13_R",
    "ACTIVE_IRI_F_B13_R",
    "ACTIVE_OX_F_B13_R",
    "PANEL_OLD_F_B13_R",
    "SAMPLE_METASTATIC_F_B13_R",
    "SPECIMEN_AGE_DAYS_B13_R",
    "MET_SITE_COUNT_B13_R"
)

missing_cols <- setdiff(
    required_cols,
    names(d)
)

if (
    length(
        missing_cols
    ) > 0L
) {
    stop(
        "Missing required B1-13 variables:\n",
        paste(
            missing_cols,
            collapse = "\n"
        )
    )
}


# ============================================================
# 3. Specimen-age variables
# ============================================================

if (
    any(
        !is.finite(
            d$SPECIMEN_AGE_DAYS_B13_R
        ) |
            d$SPECIMEN_AGE_DAYS_B13_R < 0
    )
) {
    stop(
        "Invalid specimen-age values found."
    )
}

d[
    ,
    LOG_SPECIMEN_AGE_B13C :=
        log1p(
            SPECIMEN_AGE_DAYS_B13_R
        )
]

log_age_mean <- mean(
    d$LOG_SPECIMEN_AGE_B13C,
    na.rm = TRUE
)

log_age_sd <- sd(
    d$LOG_SPECIMEN_AGE_B13C,
    na.rm = TRUE
)

if (
    !is.finite(
        log_age_sd
    ) ||
    log_age_sd <= 0
) {
    stop(
        "Could not standardize specimen age."
    )
}

d[
    ,
    LOG_SPECIMEN_AGE_Z_B13C :=
        (
            LOG_SPECIMEN_AGE_B13C -
                log_age_mean
        ) /
            log_age_sd
]


# ============================================================
# 4. Fixed Model-C architecture
# ============================================================

surv_lhs <-
    "Surv(RWPFS_DAYS_B13_R, RWPFS_EVENT_B13_R)"

base_C_covars <- c(
    "BYPASS_HIGH_CONFIDENCE_R",
    "GENDER_F_R",
    "SUBSITE_F_R",
    "ANTI_EGFR_AGENT_F_R",
    "LOG_STAGE4_TO_ANTI_EGFR_R",
    "PRIOR_MAJOR_N_F_B13_R",
    "ACTIVE_IRI_F_B13_R",
    "ACTIVE_OX_F_B13_R",
    "PANEL_OLD_F_B13_R",
    "SAMPLE_METASTATIC_F_B13_R",
    "LOG_SPECIMEN_AGE_B13C",
    "MET_SITE_COUNT_B13_R"
)

formula_C_standard <- as.formula(
    paste0(
        surv_lhs,
        " ~ ",
        paste(
            base_C_covars,
            collapse = " + "
        )
    )
)


# ============================================================
# 5. Recover the FIXED PH-addressing architecture from B1-13A
# ============================================================

if (
    !"ph_safe_specs" %in%
        names(obj) ||
    is.null(
        obj$ph_safe_specs$C
    )
) {
    stop(
        "B1-13A PH-addressed/PH-addressed specification for Model C is missing. ",
        "Run 13A_B1_13_PH_patch_and_pooled_strata.R first."
    )
}

C_spec <- obj$ph_safe_specs$C

factor_viol_C <- C_spec$factor_viol
continuous_viol_C <- C_spec$continuous_viol

# Replace the old specimen-age variable name in the fixed architecture
# if it ever appeared among the continuous violators. In the current
# B1-13A run it was not a violator, but this keeps the script robust.
continuous_viol_C[
    continuous_viol_C ==
        "LOG_SPECIMEN_AGE_B13_R"
] <- "LOG_SPECIMEN_AGE_B13C"

base_C_covars_current <- base_C_covars

factor_viol_C <- intersect(
    factor_viol_C,
    base_C_covars_current
)

continuous_viol_C <- intersect(
    continuous_viol_C,
    base_C_covars_current
)

build_fixed_ph_formula <- function(
    covars,
    factor_viol,
    continuous_viol
) {

    retained <- setdiff(
        covars,
        factor_viol
    )

    rhs <- retained

    if (
        length(
            factor_viol
        ) > 0L
    ) {

        rhs <- c(
            rhs,
            paste0(
                "strata(",
                factor_viol,
                ")"
            )
        )
    }

    if (
        length(
            continuous_viol
        ) > 0L
    ) {

        rhs <- c(
            rhs,
            paste0(
                "tt(",
                continuous_viol,
                ")"
            )
        )
    }

    as.formula(
        paste0(
            surv_lhs,
            " ~ ",
            paste(
                rhs,
                collapse = " + "
            )
        )
    )
}

formula_C_PH_fixed <- build_fixed_ph_formula(
    base_C_covars_current,
    factor_viol_C,
    continuous_viol_C
)

tt_fun <- if (
    length(
        continuous_viol_C
    ) > 0L
) {

    function(
        x,
        t,
        ...
    ) {
        x *
            log(
                pmax(
                    t,
                    1
                )
            )
    }

} else {

    NULL
}


# ============================================================
# 6. Continuous bypass x specimen-age interaction
# ============================================================

interaction_covars <- setdiff(
    base_C_covars,
    c(
        "BYPASS_HIGH_CONFIDENCE_R",
        "LOG_SPECIMEN_AGE_B13C"
    )
)

formula_interaction_standard <- as.formula(
    paste0(
        surv_lhs,
        " ~ BYPASS_HIGH_CONFIDENCE_R * ",
        "LOG_SPECIMEN_AGE_Z_B13C + ",
        paste(
            interaction_covars,
            collapse = " + "
        )
    )
)

fit_interaction_standard <- fit_cox_safe(
    formula_interaction_standard,
    d
)

interaction_term <-
    "BYPASS_HIGH_CONFIDENCE_R:LOG_SPECIMEN_AGE_Z_B13C"

interaction_standard <- extract_term(
    fit_interaction_standard,
    interaction_term
)

interaction_standard[
    ,
    MODEL :=
        "Expanded Model-C interaction"
]

interaction_standard[
    ,
    N :=
        if (
            is.null(
                fit_interaction_standard
            )
        ) {
            NA_integer_
        } else {
            fit_interaction_standard$n
        }
]

interaction_standard[
    ,
    EVENTS :=
        if (
            is.null(
                fit_interaction_standard
            )
        ) {
            NA_integer_
        } else {
            fit_interaction_standard$nevent
        }
]

# Fixed PH-addressed interaction architecture.
interaction_base_covars <- c(
    "BYPASS_HIGH_CONFIDENCE_R",
    "LOG_SPECIMEN_AGE_Z_B13C",
    interaction_covars
)

interaction_factor_viol <- intersect(
    factor_viol_C,
    interaction_base_covars
)

interaction_cont_viol <- intersect(
    continuous_viol_C,
    interaction_base_covars
)

interaction_retained <- setdiff(
    interaction_base_covars,
    interaction_factor_viol
)

interaction_rhs <- c(
    "BYPASS_HIGH_CONFIDENCE_R",
    "LOG_SPECIMEN_AGE_Z_B13C",
    "BYPASS_HIGH_CONFIDENCE_R:LOG_SPECIMEN_AGE_Z_B13C",
    setdiff(
        interaction_retained,
        c(
            "BYPASS_HIGH_CONFIDENCE_R",
            "LOG_SPECIMEN_AGE_Z_B13C"
        )
    )
)

if (
    length(
        interaction_factor_viol
    ) > 0L
) {
    interaction_rhs <- c(
        interaction_rhs,
        paste0(
            "strata(",
            interaction_factor_viol,
            ")"
        )
    )
}

if (
    length(
        interaction_cont_viol
    ) > 0L
) {
    interaction_rhs <- c(
        interaction_rhs,
        paste0(
            "tt(",
            interaction_cont_viol,
            ")"
        )
    )
}

formula_interaction_PH <- as.formula(
    paste0(
        surv_lhs,
        " ~ ",
        paste(
            unique(
                interaction_rhs
            ),
            collapse = " + "
        )
    )
)

fit_interaction_PH <- fit_cox_safe(
    formula_interaction_PH,
    d,
    tt_fun = tt_fun
)

interaction_PH <- extract_term(
    fit_interaction_PH,
    interaction_term
)

interaction_PH[
    ,
    MODEL :=
        "Fixed PH-addressed Model-C interaction"
]

interaction_PH[
    ,
    N :=
        if (
            is.null(
                fit_interaction_PH
            )
        ) {
            NA_integer_
        } else {
            fit_interaction_PH$n
        }
]

interaction_PH[
    ,
    EVENTS :=
        if (
            is.null(
                fit_interaction_PH
            )
        ) {
            NA_integer_
        } else {
            fit_interaction_PH$nevent
        }
]

interaction_results <- rbindlist(
    list(
        interaction_standard,
        interaction_PH
    ),
    fill = TRUE
)

setcolorder(
    interaction_results,
    c(
        "MODEL",
        "N",
        "EVENTS",
        "HR",
        "LCL95",
        "UCL95",
        "P",
        "SE",
        "STATUS"
    )
)


# ============================================================
# 7. Cutoff analyses
# ============================================================

cutoffs <- c(
    Inf,
    365,
    540,
    730,
    1095
)

cutoff_labels <- c(
    "Unrestricted",
    "<=365 days",
    "<=540 days",
    "<=730 days",
    "<=1095 days"
)

fit_cutoff <- function(
    cutoff,
    label
) {

    ds <- if (
        is.infinite(
            cutoff
        )
    ) {
        copy(
            d
        )
    } else {
        d[
            SPECIMEN_AGE_DAYS_B13_R <=
                cutoff
        ]
    }

    restriction_n <- nrow(
        ds
    )

    restriction_events <- sum(
        ds$RWPFS_EVENT_B13_R
    )

    restriction_pos <- ds[
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        .N
    ]

    standard_fit <- fit_cox_safe(
        formula_C_standard,
        ds
    )

    standard_res <- extract_term(
        standard_fit,
        "BYPASS_HIGH_CONFIDENCE_R"
    )

    ph_fit <- fit_cox_safe(
        formula_C_PH_fixed,
        ds,
        tt_fun = tt_fun
    )

    ph_res <- extract_term(
        ph_fit,
        "BYPASS_HIGH_CONFIDENCE_R"
    )

    data.table(
        CUTOFF = label,
        MAX_SPECIMEN_AGE_DAYS =
            if (
                is.infinite(
                    cutoff
                )
            ) {
                NA_real_
            } else {
                cutoff
            },
        RESTRICTION_N =
            restriction_n,
        RESTRICTION_EVENTS =
            restriction_events,
        RESTRICTION_BYPASS_POS =
            restriction_pos,
        STANDARD_N =
            if (
                is.null(
                    standard_fit
                )
            ) {
                NA_integer_
            } else {
                standard_fit$n
            },
        STANDARD_EVENTS =
            if (
                is.null(
                    standard_fit
                )
            ) {
                NA_integer_
            } else {
                standard_fit$nevent
            },
        STANDARD_HR =
            standard_res$HR,
        STANDARD_LCL95 =
            standard_res$LCL95,
        STANDARD_UCL95 =
            standard_res$UCL95,
        STANDARD_P =
            standard_res$P,
        PH_N =
            if (
                is.null(
                    ph_fit
                )
            ) {
                NA_integer_
            } else {
                ph_fit$n
            },
        PH_EVENTS =
            if (
                is.null(
                    ph_fit
                )
            ) {
                NA_integer_
            } else {
                ph_fit$nevent
            },
        PH_HR =
            ph_res$HR,
        PH_LCL95 =
            ph_res$LCL95,
        PH_UCL95 =
            ph_res$UCL95,
        PH_P =
            ph_res$P
    )
}

cutoff_results <- rbindlist(
    Map(
        fit_cutoff,
        cutoffs,
        cutoff_labels
    ),
    fill = TRUE
)


# ============================================================
# 8. <=730-day included vs excluded descriptive audit
# ============================================================

d[
    ,
    WITHIN_730_B13C :=
        SPECIMEN_AGE_DAYS_B13_R <=
            730
]

included_730 <- d[
    WITHIN_730_B13C == TRUE
]

excluded_730 <- d[
    WITHIN_730_B13C == FALSE
]

summary_730 <- rbindlist(
    list(
        included_730[
            ,
            .(
                GROUP =
                    "Included <=730 days",
                N = .N,
                EVENTS =
                    sum(
                        RWPFS_EVENT_B13_R
                    ),
                BYPASS_POS =
                    sum(
                        BYPASS_HIGH_CONFIDENCE_R == 1L
                    ),
                SPECIMEN_AGE_DAYS =
                    median_iqr(
                        SPECIMEN_AGE_DAYS_B13_R
                    ),
                RWPFS_DAYS =
                    median_iqr(
                        RWPFS_DAYS_B13_R
                    )
            )
        ],

        excluded_730[
            ,
            .(
                GROUP =
                    "Excluded >730 days",
                N = .N,
                EVENTS =
                    sum(
                        RWPFS_EVENT_B13_R
                    ),
                BYPASS_POS =
                    sum(
                        BYPASS_HIGH_CONFIDENCE_R == 1L
                    ),
                SPECIMEN_AGE_DAYS =
                    median_iqr(
                        SPECIMEN_AGE_DAYS_B13_R
                    ),
                RWPFS_DAYS =
                    median_iqr(
                        RWPFS_DAYS_B13_R
                    )
            )
        ]
    ),
    fill = TRUE
)

excluded_by_bypass <- excluded_730[
    ,
    .(
        N = .N,
        EVENTS =
            sum(
                RWPFS_EVENT_B13_R
            ),
        MEDIAN_RWPFS_DAYS =
            median(
                RWPFS_DAYS_B13_R,
                na.rm = TRUE
            ),
        Q1_RWPFS_DAYS =
            quantile(
                RWPFS_DAYS_B13_R,
                0.25,
                na.rm = TRUE
            ),
        Q3_RWPFS_DAYS =
            quantile(
                RWPFS_DAYS_B13_R,
                0.75,
                na.rm = TRUE
            ),
        MEDIAN_SPECIMEN_AGE_DAYS =
            median(
                SPECIMEN_AGE_DAYS_B13_R,
                na.rm = TRUE
            )
    ),
    by =
        BYPASS_HIGH_CONFIDENCE_R
]

excluded_bypass_pos <- excluded_730[
    BYPASS_HIGH_CONFIDENCE_R == 1L,
    .(
        PATIENT_ID,
        SPECIMEN_AGE_DAYS_B13_R,
        RWPFS_DAYS_B13_R,
        RWPFS_EVENT_B13_R,
        PRIOR_MAJOR_N_F_B13_R,
        ACTIVE_IRI_F_B13_R,
        ACTIVE_OX_F_B13_R,
        MET_SITE_COUNT_B13_R
    )
]

setorder(
    excluded_bypass_pos,
    RWPFS_DAYS_B13_R
)


# ============================================================
# 9. Save RDS
# ============================================================

saveRDS(
    list(
        analysis_status =
            "POST-HOC DIAGNOSTIC ONLY; NO NEW GATE",

        fixed_C_PH_architecture = list(
            factor_strata =
                factor_viol_C,
            continuous_time_interactions =
                continuous_viol_C
        ),

        specimen_age_standardization = list(
            log1p_mean =
                log_age_mean,
            log1p_sd =
                log_age_sd
        ),

        interaction_results =
            interaction_results,

        cutoff_results =
            cutoff_results,

        included_excluded_730_summary =
            summary_730,

        excluded_730_by_bypass =
            excluded_by_bypass,

        excluded_730_bypass_positive_patients =
            excluded_bypass_pos
    ),
    output_rds
)


# ============================================================
# 10. Audit TXT
# ============================================================

audit_lines <- c(
    "B1-13c SPECIMEN-TIMING DIAGNOSTIC",
    "=================================",
    "",
    "STATUS",
    "------",
    "POST-HOC diagnostic only.",
    "Technical note: tt()-containing Cox models are fit with model=FALSE; this fixes the prior NA PH-addressed outputs without changing any standard model.",
    "No new significance or wording gate.",
    "All prespecified cutoffs are reported.",
    "Model E from B1-13 remains the wording ceiling.",
    "",
    "FIXED PH-ADDRESSED MODEL-C ARCHITECTURE",
    "---------------------------------------",
    paste0(
        "Categorical nuisance terms inherited for strata(): ",
        if (
            length(
                factor_viol_C
            ) == 0L
        ) {
            "None"
        } else {
            paste(
                factor_viol_C,
                collapse = ", "
            )
        }
    ),
    paste0(
        "Continuous nuisance terms inherited for time interaction: ",
        if (
            length(
                continuous_viol_C
            ) == 0L
        ) {
            "None"
        } else {
            paste(
                continuous_viol_C,
                collapse = ", "
            )
        }
    ),
    "",
    "CONTINUOUS BYPASS x SPECIMEN-AGE INTERACTION",
    "--------------------------------------------",
    "Specimen age = standardized log1p(days); interaction HR is the ratio",
    "of the bypass HR per 1-SD older log specimen age.",
    capture.output(
        print(
            interaction_results
        )
    ),
    "",
    "PRESPECIFIED SPECIMEN-AGE CUTOFFS",
    "---------------------------------",
    capture.output(
        print(
            cutoff_results
        )
    ),
    "",
    "<=730-DAY INCLUDED VS EXCLUDED",
    "------------------------------",
    capture.output(
        print(
            summary_730
        )
    ),
    "",
    "EXCLUDED >730 DAYS BY BYPASS STATUS",
    "-----------------------------------",
    capture.output(
        print(
            excluded_by_bypass
        )
    ),
    "",
    "EXCLUDED >730-DAY BYPASS-POSITIVE PATIENTS",
    "------------------------------------------",
    "Internal diagnostic listing only; do not put patient IDs in manuscript.",
    capture.output(
        print(
            excluded_bypass_pos
        )
    ),
    "",
    "INTERPRETATION GUARDRAILS",
    "-------------------------",
    "- Non-significant interaction does not prove no specimen-age modification.",
    "- Small exposed counts limit interaction power.",
    "- Do not select a cutoff after seeing results.",
    "- Do not claim attenuation is proven to be due only to reduced precision.",
    "- If no clear monotonic pattern is observed, manuscript wording may be:",
    "  'No clear monotonic specimen-age gradient was observed, although interaction analyses were underpowered.'",
    "",
    paste0(
        "Saved RDS: ",
        output_rds
    ),
    paste0(
        "Saved table: ",
        table_file
    )
)

writeLines(
    audit_lines,
    audit_file
)


# ============================================================
# 11. Word Supplementary Table S6
# ============================================================

interaction_word <- interaction_results[
    ,
    .(
        Model = MODEL,
        N,
        Events = EVENTS,
        `Interaction HR (95% CI)` =
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

cutoff_word <- cutoff_results[
    ,
    .(
        `Specimen-age restriction` =
            CUTOFF,
        `Restriction N` =
            RESTRICTION_N,
        `Bypass+, n` =
            RESTRICTION_BYPASS_POS,
        `Standard model N` =
            STANDARD_N,
        `Standard HR (95% CI)` =
            mapply(
                fmt_hr,
                STANDARD_HR,
                STANDARD_LCL95,
                STANDARD_UCL95
            ),
        `Fixed PH-addressed N` =
            PH_N,
        `PH-addressed HR (95% CI)` =
            mapply(
                fmt_hr,
                PH_HR,
                PH_LCL95,
                PH_UCL95
            )
    )
]

summary_word <- summary_730[
    ,
    .(
        Group = GROUP,
        N,
        Events = EVENTS,
        `Bypass+, n` =
            BYPASS_POS,
        `Specimen age, days median [IQR]` =
            SPECIMEN_AGE_DAYS,
        `rwPFS, days median [IQR]` =
            RWPFS_DAYS
    )
]

ft_three_line <- function(
    df,
    font_size = 7.7
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
            "Supplementary Table S6. Post-hoc specimen-timing diagnostic analyses",
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
    "Continuous specimen-age interaction"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        interaction_word,
        7.7
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "Prespecified specimen-age restrictions"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        cutoff_word,
        7.3
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "Description of the <=730-day restriction"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        summary_word,
        7.7
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
                "Notes: These analyses are post-hoc diagnostics and do not create a new inferential or wording gate. ",
                "Specimen age was modeled continuously as standardized log1p(days); the interaction HR represents the ",
                "ratio of the bypass HR per 1-SD older log specimen age. The fixed PH-addressed models use the same ",
                "nuisance-term handling identified for the full-cohort expanded Model C in B1-13A and do not reselect ",
                "PH adjustments separately within each cutoff. All cutoffs were specified before B1-13c was run and ",
                "are reported regardless of direction. A non-significant interaction should not be interpreted as proof ",
                "of no specimen-age effect modification because the bypass-positive subgroup is small."
            ),
            fp_text(
                font.family = "Arial",
                font.size = 7.5
            )
        )
    )
)

print(
    doc,
    target = table_file
)


# ============================================================
# 12. Console output
# ============================================================

cat("\nCONTINUOUS INTERACTION\n")
cat("----------------------\n")

print(
    interaction_results
)

cat("\nCUTOFF RESULTS\n")
cat("--------------\n")

print(
    cutoff_results
)

cat("\n<=730-DAY INCLUDED VS EXCLUDED\n")
cat("-----------------------------\n")

print(
    summary_730
)

cat("\nEXCLUDED >730-DAY BYPASS-POSITIVE PATIENTS\n")
cat("-----------------------------------------\n")

print(
    excluded_bypass_pos
)

cat("\n============================================================\n")
cat("B1-13c COMPLETE\n")
cat("============================================================\n")

cat(
    "\nAudit:\n",
    audit_file,
    "\n"
)

cat(
    "\nTable:\n",
    table_file,
    "\n"
)

cat(
    "\nRDS:\n",
    output_rds,
    "\n"
)
