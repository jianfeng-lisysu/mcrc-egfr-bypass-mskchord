# ============================================================
# A7 / B-line
# 13A_B1_13_PH_patch_and_pooled_strata.R
#
# PURPOSE
#   Final internal-statistics patch for B1-13 before freezing the
#   MSK discovery analysis.
#
# THIS SCRIPT DOES FOUR THINGS ONLY
#   1) Term-level PH diagnostics for Models B, C, and E.
#   2) Adds PH-addressed counterparts without replacing the standard models.
#   3) Adds one pooled bypass estimate stratified by prior major
#      cytotoxic exposure count (0/1/2).
#   4) Fixes manuscript-facing sparse-stratum and Model-E missingness
#      reporting.
#
# IMPORTANT
#   - Original Model A remains the locked primary adjusted analysis.
#   - Standard Models B/C/E remain visible even if PH-addressed versions differ.
#   - PH-addressed construction is rule-based and declared here before running:
#
#       a) Exposure BYPASS_HIGH_CONFIDENCE_R is NEVER stratified.
#       b) A nuisance TERM with cox.zph term-level P < 0.05 is flagged.
#       c) Flagged categorical nuisance terms are moved to strata().
#       d) Flagged continuous nuisance terms retain their main effect and
#          receive a time-varying coefficient tt(x)=x*log(time).
#       e) No nuisance term is modified because of the bypass HR/P value.
#
#   - If a PH-addressed model contains tt(), cox.zph is not applicable to that
#     model; this is explicitly reported rather than hidden.
#
# SPARSE-STRATUM DISPLAY RULE
#   A within-stratum adjusted estimate is manuscript-facing "Not estimable"
#   if:
#       SE > 5
#       OR UCL95 / LCL95 > 50
#       OR CI bounds are non-finite/non-positive.
#   The raw numeric estimate remains preserved in the RDS/audit.
#
# OUTPUTS
#   OVERWRITES:
#     04_results/B1_13_additional_sensitivity_adjustment.rds
#     07_tables/TableS5_additional_sensitivity_analyses.docx
#
#   CREATES/OVERWRITES:
#     06_logs_and_audit/B1_13_PH_patch_audit.txt
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

rds_file <- file.path(
    project_root,
    "04_results",
    "B1_13_additional_sensitivity_adjustment.rds"
)

table_file <- file.path(
    project_root,
    "07_tables",
    "TableS5_additional_sensitivity_analyses.docx"
)

audit_file <- file.path(
    project_root,
    "06_logs_and_audit",
    "B1_13_PH_patch_audit.txt"
)

if (!file.exists(rds_file)) {
    stop(
        "Missing B1-13 RDS:\n",
        rds_file
    )
}

required_packages <- c(
    "data.table",
    "survival",
    "officer",
    "flextable"
)

for (pkg in required_packages) {
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
cat("B1-13A PH PATCH + POOLED PRIOR-STRATIFIED ESTIMATE\n")
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

extract_exposure <- function(
    fit,
    model_id,
    model_label,
    status = "OK",
    allow_zph = TRUE
) {

    if (is.null(fit)) {

        return(
            data.table(
                MODEL_ID = model_id,
                MODEL = model_label,
                N = NA_integer_,
                EVENTS = NA_integer_,
                BYPASS_POS = NA_integer_,
                HR = NA_real_,
                LCL95 = NA_real_,
                UCL95 = NA_real_,
                P = NA_real_,
                SE = NA_real_,
                PH_BYPASS_P = NA_real_,
                PH_GLOBAL_P = NA_real_,
                STATUS = status
            )
        )
    }

    s <- summary(fit)

    term <- "BYPASS_HIGH_CONFIDENCE_R"

    if (
        !term %in%
            rownames(
                s$coefficients
            )
    ) {

        return(
            data.table(
                MODEL_ID = model_id,
                MODEL = model_label,
                N = fit$n,
                EVENTS = fit$nevent,
                BYPASS_POS = NA_integer_,
                HR = NA_real_,
                LCL95 = NA_real_,
                UCL95 = NA_real_,
                P = NA_real_,
                SE = NA_real_,
                PH_BYPASS_P = NA_real_,
                PH_GLOBAL_P = NA_real_,
                STATUS =
                    paste0(
                        status,
                        "; bypass term not estimable"
                    )
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

    ph_bypass <- NA_real_
    ph_global <- NA_real_

    if (allow_zph) {

        z <- tryCatch(
            cox.zph(
                fit,
                transform = "km",
                terms = TRUE,
                singledf = TRUE
            ),
            error = function(e) NULL
        )

        if (!is.null(z)) {

            ztab <- as.data.frame(
                z$table
            )

            if (
                term %in%
                    rownames(
                        ztab
                    )
            ) {
                ph_bypass <-
                    ztab[
                        term,
                        "p"
                    ]
            }

            if (
                "GLOBAL" %in%
                    rownames(
                        ztab
                    )
            ) {
                ph_global <-
                    ztab[
                        "GLOBAL",
                        "p"
                    ]
            }
        }
    }

    data.table(
        MODEL_ID = model_id,
        MODEL = model_label,
        N = fit$n,
        EVENTS = fit$nevent,
        BYPASS_POS = NA_integer_,
        HR = exp(b),
        LCL95 = exp(
            b - 1.96 * se
        ),
        UCL95 = exp(
            b + 1.96 * se
        ),
        P = p,
        SE = se,
        PH_BYPASS_P = ph_bypass,
        PH_GLOBAL_P = ph_global,
        STATUS = status
    )
}

fit_basic <- function(
    formula,
    data,
    model_id,
    model_label,
    allow_zph = TRUE
) {

    mf <- tryCatch(
        model.frame(
            formula,
            data = data,
            na.action = na.omit
        ),
        error = function(e) NULL
    )

    bypass_pos <- if (
        !is.null(mf) &&
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

    fit <- tryCatch(
        suppressWarnings(
            coxph(
                formula,
                data = data,
                ties = "efron",
                na.action = na.omit,
                x = TRUE,
                y = TRUE,
                model = TRUE
            )
        ),
        error = function(e) NULL
    )

    out <- extract_exposure(
        fit,
        model_id,
        model_label,
        status =
            if (
                is.null(fit)
            ) {
                "NOT ESTIMABLE / Cox model failed"
            } else {
                "OK"
            },
        allow_zph = allow_zph
    )

    out[
        ,
        BYPASS_POS :=
            bypass_pos
    ]

    list(
        fit = fit,
        result = out
    )
}

zph_term_table <- function(
    fit,
    model_id
) {

    if (is.null(fit)) {
        return(
            data.table(
                MODEL_ID = model_id,
                TERM = NA_character_,
                CHISQ = NA_real_,
                DF = NA_real_,
                P = NA_real_
            )
        )
    }

    z <- cox.zph(
        fit,
        transform = "km",
        terms = TRUE,
        singledf = TRUE
    )

    zz <- as.data.table(
        as.data.frame(
            z$table
        ),
        keep.rownames = "TERM"
    )

    setnames(
        zz,
        old = names(zz),
        new = toupper(
            names(zz)
        )
    )

    if (
        !"CHISQ" %in%
            names(zz)
    ) {
        stop(
            "Unexpected cox.zph table structure."
        )
    }

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

    if (
        !"P" %in%
            names(zz)
    ) {
        stop(
            "Unexpected cox.zph table: P column missing."
        )
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

is_sparse_unstable <- function(
    hr,
    lo,
    hi,
    p = NA_real_
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
        ) ||
        lo <= 0 ||
        hi <= 0
    ) {
        return(TRUE)
    }

    se_from_ci <-
        (
            log(hi) -
                log(lo)
        ) /
        (
            2 *
                1.96
        )

    ci_ratio <-
        hi /
            lo

    isTRUE(
        se_from_ci > 5 ||
            ci_ratio > 50
    )
}


# ============================================================
# 2. Read B1-13 object and verify anchors
# ============================================================

obj <- readRDS(
    rds_file
)

needed_slots <- c(
    "patient_data",
    "models",
    "formulas",
    "model_results",
    "within_prior_depth_strata",
    "prior_depth_distribution"
)

missing_slots <- setdiff(
    needed_slots,
    names(obj)
)

if (
    length(
        missing_slots
    ) > 0L
) {
    stop(
        "B1-13 RDS missing slots:\n",
        paste(
            missing_slots,
            collapse = "\n"
        )
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
        "B1-13 locked anchors failed."
    )
}

if (
    is.null(
        obj$models$B
    ) ||
    is.null(
        obj$models$C
    ) ||
    is.null(
        obj$models$E
    )
) {
    stop(
        "Standard Models B/C/E are missing from B1-13 object."
    )
}


# ============================================================
# 3. Term-level PH diagnostics for standard B/C/E
# ============================================================

ph_B <- zph_term_table(
    obj$models$B,
    "B"
)

ph_C <- zph_term_table(
    obj$models$C,
    "C"
)

ph_E <- zph_term_table(
    obj$models$E,
    "E"
)

ph_terms <- rbindlist(
    list(
        ph_B,
        ph_C,
        ph_E
    ),
    fill = TRUE
)

factor_nuisance <- c(
    "GENDER_F_R",
    "SUBSITE_F_R",
    "ANTI_EGFR_AGENT_F_R",
    "PRIOR_MAJOR_N_F_B13_R",
    "ACTIVE_IRI_F_B13_R",
    "ACTIVE_OX_F_B13_R",
    "PANEL_OLD_F_B13_R",
    "SAMPLE_METASTATIC_F_B13_R"
)

continuous_nuisance <- c(
    "LOG_STAGE4_TO_ANTI_EGFR_R",
    "LOG_SPECIMEN_AGE_B13_R",
    "MET_SITE_COUNT_B13_R"
)

allowed_nuisance <- c(
    factor_nuisance,
    continuous_nuisance
)

get_violators <- function(
    ph_dt
) {

    ph_dt[
        TERM != "GLOBAL" &
            TERM !=
                "BYPASS_HIGH_CONFIDENCE_R" &
            TERM %chin%
                allowed_nuisance &
            is.finite(P) &
            P < 0.05,
        TERM
    ]
}

viol_B <- get_violators(
    ph_B
)

viol_C <- get_violators(
    ph_C
)

viol_E <- get_violators(
    ph_E
)


# ============================================================
# 4. Rule-based PH-addressed builder
# ============================================================

surv_lhs <-
    "Surv(RWPFS_DAYS_B13_R, RWPFS_EVENT_B13_R)"

covars_B <- c(
    "BYPASS_HIGH_CONFIDENCE_R",
    "GENDER_F_R",
    "SUBSITE_F_R",
    "ANTI_EGFR_AGENT_F_R",
    "LOG_STAGE4_TO_ANTI_EGFR_R",
    "PRIOR_MAJOR_N_F_B13_R",
    "ACTIVE_IRI_F_B13_R",
    "ACTIVE_OX_F_B13_R"
)

covars_C <- c(
    covars_B,
    "PANEL_OLD_F_B13_R",
    "SAMPLE_METASTATIC_F_B13_R",
    "LOG_SPECIMEN_AGE_B13_R",
    "MET_SITE_COUNT_B13_R"
)

build_ph_safe_formula <- function(
    base_covars,
    violators
) {

    factor_viol <- intersect(
        violators,
        factor_nuisance
    )

    continuous_viol <- intersect(
        violators,
        continuous_nuisance
    )

    retained_covars <- setdiff(
        base_covars,
        factor_viol
    )

    rhs_terms <- retained_covars

    if (
        length(
            factor_viol
        ) > 0L
    ) {

        rhs_terms <- c(
            rhs_terms,
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

        rhs_terms <- c(
            rhs_terms,
            paste0(
                "tt(",
                continuous_viol,
                ")"
            )
        )
    }

    f <- as.formula(
        paste0(
            surv_lhs,
            " ~ ",
            paste(
                rhs_terms,
                collapse = " + "
            )
        )
    )

    list(
        formula = f,
        factor_viol = factor_viol,
        continuous_viol = continuous_viol
    )
}

fit_ph_safe <- function(
    base_covars,
    violators,
    data,
    model_id,
    model_label
) {

    spec <- build_ph_safe_formula(
        base_covars,
        violators
    )

    has_tt <-
        length(
            spec$continuous_viol
        ) > 0L

    fit <- tryCatch(
        suppressWarnings(
            coxph(
                spec$formula,
                data = data,
                ties = "efron",
                na.action = na.omit,
                x = TRUE,
                y = TRUE,
                tt =
                    if (has_tt) {
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
            )
        ),
        error = function(e) NULL
    )

    status_parts <- c(
        "PH-addressed"
    )

    if (
        length(
            spec$factor_viol
        ) > 0L
    ) {
        status_parts <- c(
            status_parts,
            paste0(
                "stratified: ",
                paste(
                    spec$factor_viol,
                    collapse = ", "
                )
            )
        )
    }

    if (
        length(
            spec$continuous_viol
        ) > 0L
    ) {
        status_parts <- c(
            status_parts,
            paste0(
                "time interaction: ",
                paste(
                    spec$continuous_viol,
                    collapse = ", "
                )
            )
        )
    }

    if (
        length(
            violators
        ) == 0L
    ) {
        status_parts <- c(
            status_parts,
            "no nuisance term required modification"
        )
    }

    out <- extract_exposure(
        fit,
        model_id,
        model_label,
        status =
            if (
                is.null(fit)
            ) {
                "NOT ESTIMABLE / PH-addressed model failed"
            } else {
                paste(
                    status_parts,
                    collapse = "; "
                )
            },
        allow_zph =
            !has_tt
    )

    mf <- tryCatch(
        model.frame(
            spec$formula,
            data = data,
            na.action = na.omit
        ),
        error = function(e) NULL
    )

    if (
        !is.null(mf) &&
        "BYPASS_HIGH_CONFIDENCE_R" %in%
            names(mf)
    ) {
        out[
            ,
            BYPASS_POS :=
                sum(
                    mf$BYPASS_HIGH_CONFIDENCE_R == 1L,
                    na.rm = TRUE
                )
        ]
    }

    list(
        fit = fit,
        result = out,
        spec = spec
    )
}

phsafe_B <- fit_ph_safe(
    covars_B,
    viol_B,
    d,
    "B-PH",
    "Treatment-context adjustment, PH-addressed"
)

phsafe_C <- fit_ph_safe(
    covars_C,
    viol_C,
    d,
    "C-PH",
    "Expanded clinical/genomic-context adjustment, PH-addressed"
)

d_E <- d[
    SPECIMEN_WITHIN_730_B13_R == 1L
]

phsafe_E <- fit_ph_safe(
    covars_C,
    viol_E,
    d_E,
    "E-PH",
    "Specimen <=730-day expanded model, PH-addressed"
)


# ============================================================
# 5. Pooled estimate with strata(prior 0/1/2)
# ============================================================

formula_prior_strat <- as.formula(
    paste0(
        surv_lhs,
        " ~ BYPASS_HIGH_CONFIDENCE_R + ",
        "GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
        "LOG_STAGE4_TO_ANTI_EGFR_R + ",
        "ACTIVE_IRI_F_B13_R + ACTIVE_OX_F_B13_R + ",
        "strata(PRIOR_MAJOR_N_F_B13_R)"
    )
)

prior_strat_fit <- fit_basic(
    formula_prior_strat,
    d,
    "P",
    "Pooled estimate stratified by prior major cytotoxic exposure count"
)


# ============================================================
# 6. Model-E restriction vs complete-case missingness audit
# ============================================================

model_C_covars <- c(
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
    "LOG_SPECIMEN_AGE_B13_R",
    "MET_SITE_COUNT_B13_R"
)

missing_for_row <- function(
    row_index
) {

    miss <- character(0)

    for (v in model_C_covars) {

        x <- d_E[[v]][row_index]

        is_missing <- if (
            is.numeric(
                d_E[[v]]
            ) ||
            is.integer(
                d_E[[v]]
            )
        ) {
            is.na(x) ||
                !is.finite(
                    as.numeric(x)
                )
        } else {
            is.na(x) ||
                !nzchar(
                    trimws(
                        as.character(x)
                    )
                )
        }

        if (is_missing) {
            miss <- c(
                miss,
                v
            )
        }
    }

    paste(
        miss,
        collapse = "; "
    )
}

d_E[
    ,
    MISSING_C_MODEL_VARS_B13A :=
        vapply(
            seq_len(.N),
            missing_for_row,
            character(1)
        )
]

e_missing <- d_E[
    nzchar(
        MISSING_C_MODEL_VARS_B13A
    ),
    .(
        PATIENT_ID,
        BYPASS_HIGH_CONFIDENCE_R,
        RWPFS_EVENT_B13_R,
        RWPFS_DAYS_B13_R,
        MISSING_C_MODEL_VARS_B13A
    )
]

restriction_n <- nrow(
    d_E
)

restriction_events <- sum(
    d_E$RWPFS_EVENT_B13_R
)

restriction_bypass_pos <- d_E[
    BYPASS_HIGH_CONFIDENCE_R == 1L,
    .N
]

complete_n_E <- obj$model_results[
    MODEL_ID == "E",
    N
]

complete_events_E <- obj$model_results[
    MODEL_ID == "E",
    EVENTS
]


# ============================================================
# 7. Sparse-stratum display cleanup
# ============================================================

stratum_raw <- as.data.table(
    obj$within_prior_depth_strata
)

stratum_display <- copy(
    stratum_raw
)

stratum_display[
    ,
    UNSTABLE_B13A :=
        mapply(
            is_sparse_unstable,
            HR,
            LCL95,
            UCL95,
            P
        )
]

stratum_display[
    UNSTABLE_B13A == TRUE &
        grepl(
            "core-adjusted",
            ANALYSIS,
            fixed = TRUE
        ),
    `:=`(
        DISPLAY_HR =
            "Not estimable",

        DISPLAY_STATUS =
            paste0(
                "Not estimable—sparse/unstable stratum (N=",
                N,
                "; bypass+=",
                BYPASS_POS,
                ")"
            )
    )
]

stratum_display[
    is.na(
        DISPLAY_HR
    ),
    DISPLAY_HR :=
        mapply(
            fmt_hr,
            HR,
            LCL95,
            UCL95
        )
]

stratum_display[
    is.na(
        DISPLAY_STATUS
    ),
    DISPLAY_STATUS :=
        fifelse(
            STATUS == "OK",
            "OK",
            STATUS
        )
]


# ============================================================
# 8. Assemble manuscript-facing model table
# ============================================================

std <- as.data.table(
    obj$model_results
)

std_keep <- std[
    MODEL_ID %chin%
        c(
            "A",
            "B",
            "C",
            "D",
            "E"
        )
]

if (
    !"SE" %in%
        names(std_keep)
) {
    std_keep[
        ,
        SE :=
            (
                log(
                    UCL95
                ) -
                log(
                    LCL95
                )
            ) /
            (
                2 *
                    1.96
            )
    ]
}

all_models <- rbindlist(
    list(
        std_keep,
        phsafe_B$result,
        phsafe_C$result,
        phsafe_E$result,
        prior_strat_fit$result
    ),
    fill = TRUE
)

model_order <- c(
    "A",
    "B",
    "B-PH",
    "C",
    "C-PH",
    "D",
    "E",
    "E-PH",
    "P"
)

all_models[
    ,
    MODEL_ID :=
        factor(
            MODEL_ID,
            levels = model_order
        )
]

setorder(
    all_models,
    MODEL_ID
)

all_models[
    ,
    MODEL_ID :=
        as.character(
            MODEL_ID
        )
]


# ============================================================
# 9. Update RDS, preserving originals
# ============================================================

obj$term_ph_diagnostics <- ph_terms

obj$ph_safe_models <- list(
    B = phsafe_B$fit,
    C = phsafe_C$fit,
    E = phsafe_E$fit
)

obj$ph_safe_specs <- list(
    B = phsafe_B$spec,
    C = phsafe_C$spec,
    E = phsafe_E$spec
)

obj$ph_safe_results <- rbindlist(
    list(
        phsafe_B$result,
        phsafe_C$result,
        phsafe_E$result
    ),
    fill = TRUE
)

obj$prior_stratified_pooled_model <-
    prior_strat_fit$fit

obj$prior_stratified_pooled_result <-
    prior_strat_fit$result

obj$model_e_missingness <- list(
    restriction_N =
        restriction_n,
    restriction_events =
        restriction_events,
    restriction_bypass_pos =
        restriction_bypass_pos,
    complete_case_model_N =
        complete_n_E,
    complete_case_model_events =
        complete_events_E,
    excluded_due_model_covariates =
        e_missing
)

obj$within_prior_depth_strata_display <-
    stratum_display

obj$B1_13A_all_models_display <-
    all_models

saveRDS(
    obj,
    rds_file
)


# ============================================================
# 10. Audit
# ============================================================

audit_lines <- c(
    "B1-13A PH PATCH + POOLED PRIOR-STRATIFIED ESTIMATE",
    "=================================================",
    "",
    "ANALYSIS HIERARCHY",
    "------------------",
    "Model A remains the original locked primary adjusted model.",
    "Standard Models B/C/E remain reported.",
    "B-PH/C-PH/E-PH are post-hoc PH-addressed counterparts.",
    "P is a pooled bypass estimate stratified by prior major cytotoxic exposure count (0/1/2).",
    "No genomic exposure, endpoint, or primary model was redefined.",
    "",
    "TERM-LEVEL PH DIAGNOSTICS: MODEL B",
    "----------------------------------",
    capture.output(
        print(
            ph_B
        )
    ),
    "",
    "Flagged nuisance terms for B (P<0.05):",
    if (
        length(
            viol_B
        ) == 0L
    ) {
        "None"
    } else {
        paste(
            viol_B,
            collapse = ", "
        )
    },
    "",
    "TERM-LEVEL PH DIAGNOSTICS: MODEL C",
    "----------------------------------",
    capture.output(
        print(
            ph_C
        )
    ),
    "",
    "Flagged nuisance terms for C (P<0.05):",
    if (
        length(
            viol_C
        ) == 0L
    ) {
        "None"
    } else {
        paste(
            viol_C,
            collapse = ", "
        )
    },
    "",
    "TERM-LEVEL PH DIAGNOSTICS: MODEL E",
    "----------------------------------",
    capture.output(
        print(
            ph_E
        )
    ),
    "",
    "Flagged nuisance terms for E (P<0.05):",
    if (
        length(
            viol_E
        ) == 0L
    ) {
        "None"
    } else {
        paste(
            viol_E,
            collapse = ", "
        )
    },
    "",
    "STANDARD + PH-ADDRESSED + PRIOR-STRATIFIED RESULTS",
    "---------------------------------------------",
    capture.output(
        print(
            all_models
        )
    ),
    "",
    "PH-ADDRESSED CONSTRUCTION",
    "--------------------",
    paste0(
        "B-PH factor-stratified terms: ",
        if (
            length(
                phsafe_B$spec$factor_viol
            ) == 0L
        ) {
            "None"
        } else {
            paste(
                phsafe_B$spec$factor_viol,
                collapse = ", "
            )
        }
    ),
    paste0(
        "B-PH continuous time-interaction terms: ",
        if (
            length(
                phsafe_B$spec$continuous_viol
            ) == 0L
        ) {
            "None"
        } else {
            paste(
                phsafe_B$spec$continuous_viol,
                collapse = ", "
            )
        }
    ),
    paste0(
        "C-PH factor-stratified terms: ",
        if (
            length(
                phsafe_C$spec$factor_viol
            ) == 0L
        ) {
            "None"
        } else {
            paste(
                phsafe_C$spec$factor_viol,
                collapse = ", "
            )
        }
    ),
    paste0(
        "C-PH continuous time-interaction terms: ",
        if (
            length(
                phsafe_C$spec$continuous_viol
            ) == 0L
        ) {
            "None"
        } else {
            paste(
                phsafe_C$spec$continuous_viol,
                collapse = ", "
            )
        }
    ),
    paste0(
        "E-PH factor-stratified terms: ",
        if (
            length(
                phsafe_E$spec$factor_viol
            ) == 0L
        ) {
            "None"
        } else {
            paste(
                phsafe_E$spec$factor_viol,
                collapse = ", "
            )
        }
    ),
    paste0(
        "E-PH continuous time-interaction terms: ",
        if (
            length(
                phsafe_E$spec$continuous_viol
            ) == 0L
        ) {
            "None"
        } else {
            paste(
                phsafe_E$spec$continuous_viol,
                collapse = ", "
            )
        }
    ),
    "",
    "MODEL E MISSINGNESS AUDIT",
    "-------------------------",
    paste0(
        "Specimen <=730-day restriction cohort: N=",
        restriction_n,
        ", events=",
        restriction_events,
        ", bypass+=",
        restriction_bypass_pos
    ),
    paste0(
        "Standard Model E complete-case analysis: N=",
        complete_n_E,
        ", events=",
        complete_events_E
    ),
    paste0(
        "Patients excluded by Model-E covariate complete-case handling: ",
        nrow(
            e_missing
        )
    ),
    capture.output(
        print(
            e_missing
        )
    ),
    "",
    "SPARSE-STRATUM MANUSCRIPT DISPLAY",
    "---------------------------------",
    "Rule: adjusted stratum estimate is displayed as Not estimable when SE>5 or UCL/LCL>50 or CI is invalid.",
    capture.output(
        print(
            stratum_display[
                ,
                .(
                    PRIOR_MAJOR_N,
                    ANALYSIS,
                    N,
                    EVENTS,
                    BYPASS_POS,
                    BYPASS_NEG,
                    HR,
                    LCL95,
                    UCL95,
                    DISPLAY_HR,
                    DISPLAY_STATUS
                )
            ]
        )
    ),
    "",
    "INTERPRETATION GUARDRAILS",
    "-------------------------",
    "- PH-addressed models do not replace standard Models B/C/E.",
    "- A nuisance-term PH violation is handled without modifying the bypass exposure.",
    "- Pooled strata(prior012) estimate addresses treatment-depth baseline-hazard heterogeneity; it is not a formal line-of-therapy model.",
    "- Sparse within-stratum estimates are diagnostic and not independent validation.",
    "- Model E remains the predeclared wording ceiling from B1-13.",
    "- Do not use predictive/causal language without a valid comparator interaction.",
    "",
    paste0(
        "Updated RDS: ",
        rds_file
    ),
    paste0(
        "Updated Table S5: ",
        table_file
    )
)

writeLines(
    audit_lines,
    audit_file
)


# ============================================================
# 11. Rebuild manuscript-facing Table S5
# ============================================================

model_table_word <- all_models[
    ,
    .(
        Model = MODEL_ID,
        Description = MODEL,
        N,
        Events = EVENTS,
        `Bypass+, n` =
            BYPASS_POS,
        `HR (95% CI)` =
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
            ),
        `PH global P` =
            vapply(
                PH_GLOBAL_P,
                fmt_p,
                character(1)
            ),
        Status = STATUS
    )
]

stratum_table_word <- stratum_display[
    ,
    .(
        `Prior major classes` =
            PRIOR_MAJOR_N,
        Analysis = ANALYSIS,
        N,
        Events = EVENTS,
        `Bypass+, n` =
            BYPASS_POS,
        `Bypass-, n` =
            BYPASS_NEG,
        `HR (95% CI)` =
            DISPLAY_HR,
        Status =
            DISPLAY_STATUS
    )
]

ft_three_line <- function(
    df,
    font_size = 7.5
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
        j = c(
            1,
            2
        ),
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
            "Supplementary Table S5. Additional treatment-context and proportional-hazards sensitivity analyses",
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
    "A-E standard models, PH-addressed counterparts, and pooled prior-treatment-depth stratification"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        model_table_word,
        7.2
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "Within-stratum diagnostics by prior major cytotoxic exposure count"
)

doc <- body_add_flextable(
    doc,
    ft_three_line(
        stratum_table_word,
        7.2
    )
)

doc <- body_add_par(
    doc,
    ""
)

e_note <- if (
    nrow(
        e_missing
    ) == 1L
) {
    paste0(
        "The specimen <=730-day restriction yielded ",
        restriction_n,
        " patients; standard Model E included ",
        complete_n_E,
        " patients after complete-case handling because one patient lacked ",
        e_missing$MISSING_C_MODEL_VARS_B13A[1],
        "."
    )
} else if (
    nrow(
        e_missing
    ) > 1L
) {
    paste0(
        "The specimen <=730-day restriction yielded ",
        restriction_n,
        " patients; standard Model E included ",
        complete_n_E,
        " patients after complete-case handling of model covariates."
    )
} else {
    paste0(
        "The specimen <=730-day restriction and Model E both included ",
        complete_n_E,
        " patients."
    )
}

doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            paste0(
                "Notes: Model A is the original locked primary adjusted rwPFS model. ",
                "Models B-E are post-hoc sensitivity analyses. ",
                "B-PH, C-PH, and E-PH retain the corresponding standard-model population ",
                "and handle term-level proportional-hazards violations using prespecified ",
                "stratification for categorical nuisance terms and time-varying coefficients ",
                "for continuous nuisance terms; standard models remain shown. ",
                "Model P estimates the pooled bypass association while stratifying baseline ",
                "hazards by prior exposure to 0, 1, or 2 major cytotoxic classes. ",
                "The 0/1/2 measure is a treatment-depth proxy, not a formal line-of-therapy number. ",
                e_note,
                " Sparse within-stratum adjusted estimates meeting the predefined instability ",
                "criterion are displayed as not estimable, while raw values remain archived."
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
# 12. Console summary
# ============================================================

cat("\nTERM-LEVEL PH VIOLATORS\n")
cat("-----------------------\n")

cat(
    "Model B: ",
    if (
        length(
            viol_B
        ) == 0L
    ) {
        "None"
    } else {
        paste(
            viol_B,
            collapse = ", "
        )
    },
    "\n"
)

cat(
    "Model C: ",
    if (
        length(
            viol_C
        ) == 0L
    ) {
        "None"
    } else {
        paste(
            viol_C,
            collapse = ", "
        )
    },
    "\n"
)

cat(
    "Model E: ",
    if (
        length(
            viol_E
        ) == 0L
    ) {
        "None"
    } else {
        paste(
            viol_E,
            collapse = ", "
        )
    },
    "\n\n"
)

cat("STANDARD + PH-ADDRESSED + PRIOR-STRATIFIED RESULTS\n")
cat("---------------------------------------------\n")

print(
    all_models[
        ,
        .(
            MODEL_ID,
            MODEL,
            N,
            EVENTS,
            BYPASS_POS,
            HR,
            LCL95,
            UCL95,
            P,
            PH_BYPASS_P,
            PH_GLOBAL_P,
            STATUS
        )
    ]
)

cat("\nMODEL E MISSINGNESS\n")
cat("-------------------\n")

cat(
    "Restriction cohort: N=",
    restriction_n,
    ", events=",
    restriction_events,
    ", bypass+=",
    restriction_bypass_pos,
    "\n",
    sep = ""
)

cat(
    "Complete-case Model E: N=",
    complete_n_E,
    ", events=",
    complete_events_E,
    "\n",
    sep = ""
)

print(
    e_missing
)

cat("\nSPARSE-STRATUM DISPLAY CHECK\n")
cat("----------------------------\n")

print(
    stratum_display[
        ,
        .(
            PRIOR_MAJOR_N,
            ANALYSIS,
            N,
            BYPASS_POS,
            DISPLAY_HR,
            DISPLAY_STATUS
        )
    ]
)

cat("\n============================================================\n")
cat("B1-13A COMPLETE\n")
cat("============================================================\n")

cat(
    "\nAudit:\n",
    audit_file,
    "\n"
)

cat(
    "\nUpdated Table S5:\n",
    table_file,
    "\n"
)

cat(
    "\nUpdated RDS:\n",
    rds_file,
    "\n"
)
