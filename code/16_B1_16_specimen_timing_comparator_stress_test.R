# ============================================================
# A7 / B-line
# 16_B1_16_specimen_timing_comparator_stress_test.R
#
# POST-HOC DIAGNOSTIC ONLY (AFTER OUTCOME ACCESS)
# NO NEW PRIMARY ANALYSIS / NO NEW EXPOSURE / NO CUTOFF SHOPPING
#
# PURPOSE
#   Address one specific question:
#   Could the exploratory anti-EGFR vs bevacizumab treatment-by-bypass
#   interaction be materially explained by between-arm differences in
#   pretreatment specimen age?
#
# FROZEN INPUTS
#   04_results/B1_14_bevacizumab_comparator_feasibility.rds
#   04_results/B1_14B_exploratory_comparator_rwPFS.rds
#   04_results/B1_13c_specimen_timing_diagnostic.rds
#
# ANALYSES ALLOWED IN B1-16
#   A) SHARED <=730-DAY SPECIMEN RESTRICTION
#      - Apply the SAME <=730-day specimen-to-biologic restriction to
#        anti-EGFR and bevacizumab.
#      - Refit the EXACT frozen B1-14 logistic propensity model.
#      - Recalculate overlap weights from that refitted PS.
#      - Refit the EXACT frozen B1-14B M1 and M2 models.
#      - Robust sandwich variance; patient clustering unchanged.
#      - No trimming to empirical common support.
#
#   B) BEVACIZUMAB BYPASS+ LEAVE-ONE-OUT INFLUENCE DIAGNOSTIC
#      - Start from the frozen 319-patient PS-complete B1-14B population.
#      - There must be exactly 14 bevacizumab bypass+ patients.
#      - Delete one such patient at a time.
#      - In EACH replicate: refit frozen PS -> recalculate overlap weights
#        -> refit frozen M1.
#      - Report all 14 interaction HRs; do not select a favorable replicate.
#
#   C) EXISTING OUTCOME-BLIND COMPONENT COMPOSITION
#      - Read directly from B1-14 outcome-blind component_summary.
#      - No outcome-based gene/component analysis is performed.
#
#   D) PREEXISTING SPECIMEN-TIMING CONTEXT
#      - Read B1-13c results only.
#      - No new 365/540/1095 comparator models.
#      - No non-nested interval HRs.
#      - No treatment x bypass x specimen-age three-way interaction.
#      - No Stage-IV-to-specimen covariate added post hoc.
#      - No <=365 patient-level outcome excavation.
#
# ------------------------------------------------------------
# MANUSCRIPT FRAMING RULES — FROZEN BEFORE B1-16 RESULTS
# ------------------------------------------------------------
# These rules govern PRESENTATION only. They are not new hypothesis tests
# and never replace the frozen primary anti-EGFR rwPFS result.
#
# TIER A — "MAIN FIGURE / TITLE ELIGIBLE"
#   Restricted M1 interaction HR >= 2.50 AND 95% CI lower bound > 1.
#   -> Figure 4 may remain a main figure.
#   -> A title mentioning exploratory treatment heterogeneity is eligible,
#      but not mandatory.
#   -> Comparator remains observational, exploratory, and noncausal.
#
# TIER B — "MAIN FIGURE / TITLE NO HETEROGENEITY"
#   Restricted M1 HR > 1, but Tier A is not met, AND
#   (restricted M1 HR >= 1.50 OR 95% CI lower bound > 1).
#   -> Figure 4 may remain in the main manuscript.
#   -> Remove treatment heterogeneity from the title.
#   -> Comparator is a secondary exploratory finding.
#
# TIER C — "DOWNGRADE COMPARATOR"
#   Restricted M1 HR <= 1, OR
#   (restricted M1 HR < 1.50 AND 95% CI lower bound <= 1).
#   -> Remove treatment heterogeneity from the title.
#   -> Figure 4 should be considered for Supplementary placement.
#   -> Comparator becomes exploratory/uncertain context, not a central claim.
#
# LOO INFLUENCE RULE
#   LOO does NOT automatically change Tier A/B/C.
#   It is reported separately:
#     - if all 14 successful LOO interaction HRs remain >1:
#         "no single-patient directional reversal observed"
#     - if any successful LOO interaction HR <=1:
#         "single-patient directional reversal observed"
#     - any failed replicate is explicitly reported.
#
# SPECIMEN-TIMING WORDING — REGARDLESS OF B1-16 RESULT
#   "The magnitude of the anti-EGFR association was sensitive to specimen
#    recency, with smaller and less precise estimates under more restrictive
#    specimen-age windows; continuous interaction analyses did not establish
#    effect modification."
#
# NEVER SAY
#   - predictive biomarker (as established conclusion)
#   - causal interaction
#   - treatment-selection rule
#   - independent validation
#   - specimen timing did not matter
#   - effect modification was proved/disproved
#
# OUTPUTS — SAME FILENAMES OVERWRITTEN ON RERUN
#   04_results/B1_16_specimen_timing_comparator_stress_test.rds
#   06_logs_and_audit/B1_16_specimen_timing_comparator_stress_test.txt
#   07_tables/B1_16_specimen_timing_comparator_stress_test.docx
#
# IMPORTANT
#   This diagnostic Word file is NOT assigned a Supplementary Table number.
#   After the result is interpreted under the frozen rules, selected panels
#   can be merged into the final supplementary package without renumbering
#   the current S1-S9 prematurely.
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

results_dir <- file.path(project_root, "04_results")
audit_dir   <- file.path(project_root, "06_logs_and_audit")
table_dir   <- file.path(project_root, "07_tables")

b114_file <- file.path(
    results_dir,
    "B1_14_bevacizumab_comparator_feasibility.rds"
)

b114b_file <- file.path(
    results_dir,
    "B1_14B_exploratory_comparator_rwPFS.rds"
)

b113c_file <- file.path(
    results_dir,
    "B1_13c_specimen_timing_diagnostic.rds"
)

out_rds <- file.path(
    results_dir,
    "B1_16_specimen_timing_comparator_stress_test.rds"
)

out_txt <- file.path(
    audit_dir,
    "B1_16_specimen_timing_comparator_stress_test.txt"
)

out_doc <- file.path(
    table_dir,
    "B1_16_specimen_timing_comparator_stress_test.docx"
)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
    "data.table",
    "survival",
    "officer",
    "flextable"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]

if (length(missing_packages) > 0L) {
    stop(
        "Missing required packages: ",
        paste(missing_packages, collapse = ", ")
    )
}

library(data.table)
library(survival)
library(officer)
library(flextable)

for (f in c(b114_file, b114b_file, b113c_file)) {
    if (!file.exists(f)) {
        stop("Missing required frozen input:\n", f)
    }
}

# ============================================================
# 1. Helpers
# ============================================================

fmt_num <- function(x, digits = 2L) {
    ifelse(
        is.finite(x),
        formatC(x, format = "f", digits = digits),
        NA_character_
    )
}

fmt_p <- function(p) {
    if (!is.finite(p)) {
        return(NA_character_)
    }
    if (p < 0.001) {
        return("<.001")
    }
    sub("^0", "", sprintf("%.3f", p))
}

fmt_hr_ci <- function(hr, lo, hi, digits = 2L) {
    if (
        !is.finite(hr) ||
        !is.finite(lo) ||
        !is.finite(hi)
    ) {
        return(NA_character_)
    }

    paste0(
        fmt_num(hr, digits),
        " (",
        fmt_num(lo, digits),
        "–",
        fmt_num(hi, digits),
        ")"
    )
}

weighted_mean_safe <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) {
        return(NA_real_)
    }
    sum(w[ok] * x[ok]) / sum(w[ok])
}

weighted_var_safe <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (sum(ok) < 2L) {
        return(NA_real_)
    }
    mu <- weighted_mean_safe(x[ok], w[ok])
    sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok])
}

smd_one <- function(x, z, w = NULL) {
    if (is.null(w)) {
        w <- rep(1, length(x))
    }

    ok <- is.finite(x) & !is.na(z) & is.finite(w) & w > 0
    x <- x[ok]
    z <- z[ok]
    w <- w[ok]

    if (
        length(x) == 0L ||
        !all(z %in% c(0, 1)) ||
        !any(z == 0) ||
        !any(z == 1)
    ) {
        return(NA_real_)
    }

    m1 <- weighted_mean_safe(x[z == 1], w[z == 1])
    m0 <- weighted_mean_safe(x[z == 0], w[z == 0])

    v1 <- weighted_var_safe(x[z == 1], w[z == 1])
    v0 <- weighted_var_safe(x[z == 0], w[z == 0])

    denom <- sqrt((v1 + v0) / 2)

    if (!is.finite(denom) || denom < .Machine$double.eps) {
        if (isTRUE(all.equal(m1, m0, tolerance = 1e-12))) {
            return(0)
        }
        return(NA_real_)
    }

    (m1 - m0) / denom
}

extract_coef <- function(fit, term) {
    b <- coef(fit)

    if (!(term %in% names(b))) {
        stop("Term not found in Cox model: ", term)
    }

    V <- vcov(fit)
    beta <- unname(b[term])
    se <- sqrt(unname(V[term, term]))

    data.table(
        TERM = term,
        BETA = beta,
        SE = se,
        HR = exp(beta),
        LCL95 = exp(beta - 1.96 * se),
        UCL95 = exp(beta + 1.96 * se),
        P = 2 * pnorm(-abs(beta / se))
    )
}

linear_combo <- function(fit, terms) {
    b <- coef(fit)
    V <- vcov(fit)

    missing_terms <- setdiff(terms, names(b))
    if (length(missing_terms) > 0L) {
        stop(
            "Linear-combination term(s) absent: ",
            paste(missing_terms, collapse = ", ")
        )
    }

    L <- rep(0, length(b))
    names(L) <- names(b)
    L[terms] <- 1

    beta <- sum(L * b)
    se <- sqrt(
        as.numeric(
            t(L) %*% V %*% L
        )
    )

    data.table(
        BETA = beta,
        SE = se,
        HR = exp(beta),
        LCL95 = exp(beta - 1.96 * se),
        UCL95 = exp(beta + 1.96 * se),
        P = 2 * pnorm(-abs(beta / se))
    )
}

safe_cox_zph <- function(fit, model_id) {
    z <- tryCatch(
        cox.zph(fit),
        error = function(e) e
    )

    if (inherits(z, "error")) {
        return(
            data.table(
                MODEL_ID = model_id,
                TERM = "COX.ZPH_FAILED",
                CHISQ = NA_real_,
                DF = NA_real_,
                P = NA_real_,
                MESSAGE = conditionMessage(z)
            )
        )
    }

    tab <- as.data.table(
        as.data.frame(z$table),
        keep.rownames = "TERM"
    )

    setnames(
        tab,
        old = intersect(
            c("chisq", "df", "p"),
            names(tab)
        ),
        new = c(
            "CHISQ",
            "DF",
            "P"
        )[seq_along(
            intersect(
                c("chisq", "df", "p"),
                names(tab)
            )
        )]
    )

    if (!("CHISQ" %in% names(tab))) tab[, CHISQ := NA_real_]
    if (!("DF" %in% names(tab)))    tab[, DF := NA_real_]
    if (!("P" %in% names(tab)))     tab[, P := NA_real_]

    tab[
        ,
        `:=`(
            MODEL_ID = model_id,
            MESSAGE = NA_character_
        )
    ]

    setcolorder(
        tab,
        c(
            "MODEL_ID",
            "TERM",
            "CHISQ",
            "DF",
            "P",
            "MESSAGE"
        )
    )

    tab[]
}

three_line <- function(x, size = 7.4) {
    ft <- flextable(x)
    ft <- border_remove(ft)

    border_top <- fp_border(
        color = "black",
        width = 1
    )

    border_mid <- fp_border(
        color = "black",
        width = 0.6
    )

    ft <- hline_top(
        ft,
        border = border_top,
        part = "header"
    )

    ft <- hline_bottom(
        ft,
        border = border_mid,
        part = "header"
    )

    ft <- hline_bottom(
        ft,
        border = border_top,
        part = "body"
    )

    ft <- font(
        ft,
        fontname = "Arial",
        part = "all"
    )

    ft <- fontsize(
        ft,
        size = size,
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

    ft <- autofit(ft)
    ft
}

# ============================================================
# 2. Load frozen objects and hard-anchor them
# ============================================================

b114 <- readRDS(b114_file)
b114b <- readRDS(b114b_file)
b113c <- readRDS(b113c_file)

required_b114 <- c(
    "patient_audit",
    "ps_formula",
    "gate_decision",
    "component_summary"
)

missing_b114 <- setdiff(
    required_b114,
    names(b114)
)

if (length(missing_b114) > 0L) {
    stop(
        "B1-14 object missing slot(s): ",
        paste(missing_b114, collapse = ", ")
    )
}

if (!identical(as.character(b114$gate_decision), "GO")) {
    stop("B1-14 gate_decision is not GO.")
}

required_b114b <- c(
    "analysis_data",
    "model_contrasts",
    "interaction_gate"
)

missing_b114b <- setdiff(
    required_b114b,
    names(b114b)
)

if (length(missing_b114b) > 0L) {
    stop(
        "B1-14B object missing slot(s): ",
        paste(missing_b114b, collapse = ", ")
    )
}

d0 <- as.data.table(
    copy(
        b114b$analysis_data
    )
)

required_vars <- c(
    "PATIENT_ID",
    "FIRST_BIOLOGIC_R",
    "BYPASS_HIGH_CONFIDENCE_R",
    "SPECIMEN_TO_BIOLOGIC_DAYS_R",
    "RWPFS_DAYS_B14B_R",
    "RWPFS_EVENT_B14B_R",
    "GENDER_F_SHARED_R",
    "SUBSITE_F_SHARED_R",
    "LOG_STAGE4_TO_BIOLOGIC_R",
    "PRIOR_MAJOR_N_F_SHARED_R",
    "ACTIVE_IRI_F_SHARED_R",
    "ACTIVE_OX_F_SHARED_R"
)

ps_vars <- all.vars(
    b114$ps_formula
)

required_vars <- unique(
    c(
        required_vars,
        ps_vars
    )
)

missing_vars <- setdiff(
    required_vars,
    names(d0)
)

if (length(missing_vars) > 0L) {
    stop(
        "Frozen B1-14B analysis_data missing variable(s):\n",
        paste(missing_vars, collapse = "\n")
    )
}

d0[
    ,
    ANTI_EGFR_TREAT_R :=
        as.integer(
            FIRST_BIOLOGIC_R ==
                "anti-EGFR"
        )
]

hard_anchor <- data.table(
    Metric = c(
        "PS-complete N",
        "rwPFS events",
        "Anti-EGFR N",
        "Bevacizumab N",
        "Anti-EGFR bypass+",
        "Bevacizumab bypass+"
    ),
    Observed = c(
        nrow(d0),
        sum(d0$RWPFS_EVENT_B14B_R),
        d0[
            FIRST_BIOLOGIC_R ==
                "anti-EGFR",
            .N
        ],
        d0[
            FIRST_BIOLOGIC_R ==
                "bevacizumab",
            .N
        ],
        d0[
            FIRST_BIOLOGIC_R ==
                "anti-EGFR" &
            BYPASS_HIGH_CONFIDENCE_R ==
                1L,
            .N
        ],
        d0[
            FIRST_BIOLOGIC_R ==
                "bevacizumab" &
            BYPASS_HIGH_CONFIDENCE_R ==
                1L,
            .N
        ]
    ),
    Expected = c(
        319L,
        305L,
        188L,
        131L,
        23L,
        14L
    )
)

hard_anchor[
    ,
    PASS :=
        Observed ==
            Expected
]

if (any(!hard_anchor$PASS)) {
    print(hard_anchor)
    stop("B1-16 hard anchors failed. Do not proceed.")
}

original_m1 <- as.data.table(
    copy(
        b114b$model_contrasts
    )
)[
    MODEL_ID == "M1"
]

if (
    nrow(original_m1) != 1L ||
    !is.finite(original_m1$INTERACTION_HR)
) {
    stop("Could not recover frozen original B1-14B M1 interaction.")
}

if (
    abs(
        original_m1$INTERACTION_HR -
            3.08
    ) >
        0.10
) {
    stop(
        "Frozen original M1 interaction anchor differs materially from 3.08."
    )
}

# ============================================================
# 3. Exact frozen propensity-score / overlap-weight refit helper
# ============================================================

ps_covars <- setdiff(
    all.vars(b114$ps_formula),
    "FIRST_BIOLOGIC_R"
)

fit_ps_ow <- function(din, label = "") {

    d <- as.data.table(
        copy(din)
    )

    d[
        ,
        ANTI_EGFR_TREAT_R :=
            as.integer(
                FIRST_BIOLOGIC_R ==
                    "anti-EGFR"
            )
    ]

    cc <- complete.cases(
        d[
            ,
            ..ps_covars
        ]
    )

    d <- d[
        cc
    ]

    if (
        !any(
            d$FIRST_BIOLOGIC_R ==
                "anti-EGFR"
        ) ||
        !any(
            d$FIRST_BIOLOGIC_R ==
                "bevacizumab"
        )
    ) {
        stop(
            "Both treatment arms are required in PS refit: ",
            label
        )
    }

    ps_fit <- glm(
        b114$ps_formula,
        data = d,
        family = binomial()
    )

    if (
        !isTRUE(ps_fit$converged) ||
        any(
            !is.finite(
                coef(ps_fit)
            )
        )
    ) {
        stop(
            "Propensity model failed: ",
            label
        )
    }

    d[
        ,
        PS_B116_R :=
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
                d$PS_B116_R
            )
        ) ||
        any(
            d$PS_B116_R <=
                0 |
            d$PS_B116_R >=
                1
        )
    ) {
        stop(
            "Invalid propensity score(s): ",
            label
        )
    }

    d[
        ,
        OVERLAP_WEIGHT_B116_R :=
            fifelse(
                ANTI_EGFR_TREAT_R ==
                    1L,
                1 -
                    PS_B116_R,
                PS_B116_R
            )
    ]

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
                x <- X[, v]

                data.table(
                    Covariate = v,
                    SMD_unweighted =
                        smd_one(
                            x,
                            d$ANTI_EGFR_TREAT_R
                        ),
                    SMD_overlap_weighted =
                        smd_one(
                            x,
                            d$ANTI_EGFR_TREAT_R,
                            d$OVERLAP_WEIGHT_B116_R
                        )
                )
            }
        )
    )

    balance[
        ,
        ABS_SMD_overlap_weighted :=
            abs(
                SMD_overlap_weighted
            )
    ]

    max_abs_smd <- max(
        balance$ABS_SMD_overlap_weighted,
        na.rm = TRUE
    )

    list(
        data = d,
        ps_model = ps_fit,
        balance = balance,
        max_abs_smd = max_abs_smd
    )
}

# ============================================================
# 4. Frozen M1/M2 model helper
# ============================================================

term_treatment <-
    "ANTI_EGFR_TREAT_R"

term_bypass <-
    "BYPASS_HIGH_CONFIDENCE_R"

term_interaction <-
    "ANTI_EGFR_TREAT_R:BYPASS_HIGH_CONFIDENCE_R"

fit_models <- function(dw) {

    fit_M1 <- coxph(
        Surv(
            RWPFS_DAYS_B14B_R,
            RWPFS_EVENT_B14B_R
        ) ~
            ANTI_EGFR_TREAT_R *
            BYPASS_HIGH_CONFIDENCE_R,
        data = dw,
        weights =
            OVERLAP_WEIGHT_B116_R,
        cluster =
            PATIENT_ID,
        robust = TRUE,
        ties = "efron",
        x = TRUE,
        y = TRUE
    )

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
            OVERLAP_WEIGHT_B116_R,
        cluster =
            PATIENT_ID,
        robust = TRUE,
        ties = "efron",
        x = TRUE,
        y = TRUE
    )

    extract_one <- function(fit, id, label) {

        inter <- extract_coef(
            fit,
            term_interaction
        )

        t_neg <- linear_combo(
            fit,
            term_treatment
        )

        t_pos <- linear_combo(
            fit,
            c(
                term_treatment,
                term_interaction
            )
        )

        b_bev <- linear_combo(
            fit,
            term_bypass
        )

        b_anti <- linear_combo(
            fit,
            c(
                term_bypass,
                term_interaction
            )
        )

        list(
            interaction =
                data.table(
                    MODEL_ID = id,
                    MODEL = label,
                    N = fit$n,
                    EVENTS = fit$nevent,
                    INTERACTION_HR = inter$HR,
                    LCL95 = inter$LCL95,
                    UCL95 = inter$UCL95,
                    P = inter$P
                ),
            contrasts =
                data.table(
                    MODEL_ID = id,
                    Contrast = c(
                        "Anti-EGFR vs bevacizumab HR in bypass−",
                        "Anti-EGFR vs bevacizumab HR in bypass+",
                        "Bypass HR within bevacizumab",
                        "Bypass HR within anti-EGFR"
                    ),
                    HR = c(
                        t_neg$HR,
                        t_pos$HR,
                        b_bev$HR,
                        b_anti$HR
                    ),
                    LCL95 = c(
                        t_neg$LCL95,
                        t_pos$LCL95,
                        b_bev$LCL95,
                        b_anti$LCL95
                    ),
                    UCL95 = c(
                        t_neg$UCL95,
                        t_pos$UCL95,
                        b_bev$UCL95,
                        b_anti$UCL95
                    ),
                    P = c(
                        t_neg$P,
                        t_pos$P,
                        b_bev$P,
                        b_anti$P
                    )
                )
        )
    }

    m1 <- extract_one(
        fit_M1,
        "M1",
        "Overlap-weighted interaction"
    )

    m2 <- extract_one(
        fit_M2,
        "M2",
        "Overlap-weighted doubly adjusted interaction"
    )

    list(
        M1 = fit_M1,
        M2 = fit_M2,
        interactions =
            rbindlist(
                list(
                    m1$interaction,
                    m2$interaction
                )
            ),
        contrasts =
            rbindlist(
                list(
                    m1$contrasts,
                    m2$contrasts
                )
            ),
        ph =
            rbindlist(
                list(
                    safe_cox_zph(
                        fit_M1,
                        "M1"
                    ),
                    safe_cox_zph(
                        fit_M2,
                        "M2"
                    )
                ),
                fill = TRUE
            )
    )
}

# ============================================================
# 5. A — shared <=730-day comparator restriction
# ============================================================

d730_source <- d0[
    is.finite(
        SPECIMEN_TO_BIOLOGIC_DAYS_R
    ) &
    SPECIMEN_TO_BIOLOGIC_DAYS_R >=
        0 &
    SPECIMEN_TO_BIOLOGIC_DAYS_R <=
        730
]

if (
    !any(
        d730_source$FIRST_BIOLOGIC_R ==
            "anti-EGFR"
    ) ||
    !any(
        d730_source$FIRST_BIOLOGIC_R ==
            "bevacizumab"
    )
) {
    stop("<=730-day restriction removed one treatment arm.")
}

ps730 <- fit_ps_ow(
    d730_source,
    "<=730-day shared restriction"
)

d730 <- ps730$data

if (
    !is.finite(
        ps730$max_abs_smd
    ) ||
    ps730$max_abs_smd >
        0.10
) {
    stop(
        "<=730-day refitted overlap-weight balance gate failed: max |SMD| = ",
        ps730$max_abs_smd
    )
}

models730 <- fit_models(
    d730
)

cell730 <- d730[
    ,
    .(
        N = .N,
        Events =
            sum(
                RWPFS_EVENT_B14B_R
            ),
        `Bypass+, n` =
            sum(
                BYPASS_HIGH_CONFIDENCE_R ==
                    1L
            ),
        `Specimen age, median [IQR]` =
            paste0(
                fmt_num(
                    median(
                        SPECIMEN_TO_BIOLOGIC_DAYS_R,
                        na.rm = TRUE
                    ),
                    1
                ),
                " [",
                fmt_num(
                    quantile(
                        SPECIMEN_TO_BIOLOGIC_DAYS_R,
                        0.25,
                        na.rm = TRUE,
                        names = FALSE
                    ),
                    1
                ),
                ", ",
                fmt_num(
                    quantile(
                        SPECIMEN_TO_BIOLOGIC_DAYS_R,
                        0.75,
                        na.rm = TRUE,
                        names = FALSE
                    ),
                    1
                ),
                "]"
            )
    ),
    by = FIRST_BIOLOGIC_R
]

setnames(
    cell730,
    "FIRST_BIOLOGIC_R",
    "Treatment"
)

m1_730 <- models730$interactions[
    MODEL_ID ==
        "M1"
]

if (
    nrow(m1_730) !=
        1L
) {
    stop("Restricted M1 extraction failed.")
}

# ============================================================
# 6. Frozen manuscript framing gate
# ============================================================

restricted_hr <- m1_730$INTERACTION_HR
restricted_lo <- m1_730$LCL95
restricted_hi <- m1_730$UCL95

if (
    is.finite(restricted_hr) &&
    is.finite(restricted_lo) &&
    restricted_hr >= 2.50 &&
    restricted_lo > 1
) {

    framing_tier <-
        "TIER A — MAIN FIGURE / TITLE ELIGIBLE"

    figure4_action <-
        "Figure 4 may remain a main figure."

    title_action <-
        "A title mentioning exploratory treatment heterogeneity is eligible but not mandatory."

    comparator_action <-
        "Comparator remains exploratory, observational, and noncausal."

} else if (
    is.finite(restricted_hr) &&
    is.finite(restricted_lo) &&
    restricted_hr > 1 &&
    (
        restricted_hr >= 1.50 ||
        restricted_lo > 1
    )
) {

    framing_tier <-
        "TIER B — MAIN FIGURE / TITLE NO HETEROGENEITY"

    figure4_action <-
        "Figure 4 may remain in the main manuscript."

    title_action <-
        "Remove treatment heterogeneity from the title."

    comparator_action <-
        "Comparator is a secondary exploratory finding."

} else {

    framing_tier <-
        "TIER C — DOWNGRADE COMPARATOR"

    figure4_action <-
        "Consider moving Figure 4 to Supplementary material."

    title_action <-
        "Remove treatment heterogeneity from the title."

    comparator_action <-
        "Comparator becomes exploratory/uncertain context and not a central claim."
}

framing_table <- data.table(
    Item = c(
        "Frozen framing tier",
        "Figure 4",
        "Title",
        "Comparator wording"
    ),
    Decision = c(
        framing_tier,
        figure4_action,
        title_action,
        comparator_action
    )
)

# ============================================================
# 7. B — bevacizumab bypass+ leave-one-out full PS->OW->M1
# ============================================================

bev_pos_ids <- sort(
    unique(
        as.character(
            d0[
                FIRST_BIOLOGIC_R ==
                    "bevacizumab" &
                BYPASS_HIGH_CONFIDENCE_R ==
                    1L,
                PATIENT_ID
            ]
        )
    )
)

if (
    length(
        bev_pos_ids
    ) !=
        14L
) {
    stop(
        "Expected exactly 14 bevacizumab bypass+ PS-complete patients; observed ",
        length(bev_pos_ids)
    )
}

loo_results <- rbindlist(
    lapply(
        seq_along(
            bev_pos_ids
        ),
        function(i) {

            id_i <- bev_pos_ids[i]

            result <- tryCatch(
                {
                    di <- d0[
                        as.character(
                            PATIENT_ID
                        ) !=
                            id_i
                    ]

                    psi <- fit_ps_ow(
                        di,
                        paste0(
                            "BEV bypass+ LOO ",
                            i
                        )
                    )

                    dwi <- psi$data

                    fit_i <- coxph(
                        Surv(
                            RWPFS_DAYS_B14B_R,
                            RWPFS_EVENT_B14B_R
                        ) ~
                            ANTI_EGFR_TREAT_R *
                            BYPASS_HIGH_CONFIDENCE_R,
                        data = dwi,
                        weights =
                            OVERLAP_WEIGHT_B116_R,
                        cluster =
                            PATIENT_ID,
                        robust = TRUE,
                        ties = "efron",
                        x = TRUE,
                        y = TRUE
                    )

                    inter_i <- extract_coef(
                        fit_i,
                        term_interaction
                    )

                    data.table(
                        LOO = sprintf(
                            "LOO_%02d",
                            i
                        ),
                        PATIENT_ID_INTERNAL = id_i,
                        STATUS = "SUCCESS",
                        N = fit_i$n,
                        EVENTS = fit_i$nevent,
                        INTERACTION_HR =
                            inter_i$HR,
                        LCL95 =
                            inter_i$LCL95,
                        UCL95 =
                            inter_i$UCL95,
                        P =
                            inter_i$P,
                        MAX_ABS_WEIGHTED_SMD =
                            psi$max_abs_smd,
                        MESSAGE =
                            NA_character_
                    )
                },
                error = function(e) {
                    data.table(
                        LOO = sprintf(
                            "LOO_%02d",
                            i
                        ),
                        PATIENT_ID_INTERNAL = id_i,
                        STATUS = "FAILED",
                        N = NA_integer_,
                        EVENTS = NA_integer_,
                        INTERACTION_HR =
                            NA_real_,
                        LCL95 =
                            NA_real_,
                        UCL95 =
                            NA_real_,
                        P =
                            NA_real_,
                        MAX_ABS_WEIGHTED_SMD =
                            NA_real_,
                        MESSAGE =
                            conditionMessage(e)
                    )
                }
            )

            result
        }
    ),
    fill = TRUE
)

loo_success <- loo_results[
    STATUS ==
        "SUCCESS"
]

loo_failure_n <- loo_results[
    STATUS !=
        "SUCCESS",
    .N
]

loo_directional_reversal <- any(
    loo_success$INTERACTION_HR <=
        1,
    na.rm = TRUE
)

if (
    loo_failure_n >
        0L
) {

    loo_interpretation <-
        paste0(
            "INCOMPLETE: ",
            loo_failure_n,
            " of 14 LOO replicates failed; failures are reported and no favorable subset is selected."
        )

} else if (
    loo_directional_reversal
) {

    loo_interpretation <-
        "Single-patient directional reversal observed: at least one successful LOO interaction HR was <=1."

} else {

    loo_interpretation <-
        "No single-patient directional reversal observed: all 14 successful LOO interaction HRs remained >1."
}

loo_summary <- data.table(
    Metric = c(
        "Planned LOO replicates",
        "Successful replicates",
        "Failed replicates",
        "Minimum interaction HR",
        "Maximum interaction HR",
        "Replicates with HR >1",
        "Replicates with 95% CI lower bound >1",
        "Interpretation"
    ),
    Value = c(
        "14",
        as.character(
            nrow(
                loo_success
            )
        ),
        as.character(
            loo_failure_n
        ),
        if (
            nrow(
                loo_success
            ) >
                0L
        ) {
            fmt_num(
                min(
                    loo_success$INTERACTION_HR,
                    na.rm = TRUE
                ),
                2
            )
        } else {
            NA_character_
        },
        if (
            nrow(
                loo_success
            ) >
                0L
        ) {
            fmt_num(
                max(
                    loo_success$INTERACTION_HR,
                    na.rm = TRUE
                ),
                2
            )
        } else {
            NA_character_
        },
        as.character(
            sum(
                loo_success$INTERACTION_HR >
                    1,
                na.rm = TRUE
            )
        ),
        as.character(
            sum(
                loo_success$LCL95 >
                    1,
                na.rm = TRUE
            )
        ),
        loo_interpretation
    )
)

# ============================================================
# 8. C — existing outcome-blind component composition
# ============================================================

component <- as.data.table(
    copy(
        b114$component_summary
    )
)

required_component_cols <- c(
    "Treatment",
    "Component",
    "N"
)

missing_component_cols <- setdiff(
    required_component_cols,
    names(component)
)

if (
    length(
        missing_component_cols
    ) >
        0L
) {
    stop(
        "B1-14 component_summary missing column(s): ",
        paste(
            missing_component_cols,
            collapse = ", "
        )
    )
}

component[
    ,
    HC_DENOM :=
        fifelse(
            Treatment ==
                "anti-EGFR",
            23,
            16
        )
]

component[
    ,
    `Percent of HC+ carriers` :=
        100 *
        N /
        HC_DENOM
]

component_word <- component[
    ,
    .(
        Treatment,
        Component,
        `Component carriers, n` =
            N,
        `Percent of HC+ carriers` =
            sprintf(
                "%.1f%%",
                `Percent of HC+ carriers`
            )
    )
]

# ============================================================
# 9. D — read existing B1-13c specimen-timing context only
# ============================================================

s6_existing <- list(
    interaction_results =
        if (
            "interaction_results" %in%
                names(b113c)
        ) {
            b113c$interaction_results
        } else {
            NULL
        },
    cutoff_results =
        if (
            "cutoff_results" %in%
                names(b113c)
        ) {
            b113c$cutoff_results
        } else {
            NULL
        },
    included_excluded_730_summary =
        if (
            "included_excluded_730_summary" %in%
                names(b113c)
        ) {
            b113c$included_excluded_730_summary
        } else {
            NULL
        }
)

# ============================================================
# 10. Publication/audit tables
# ============================================================

interaction_word <- copy(
    models730$interactions
)

interaction_word[
    ,
    `Interaction HR (95% CI)` :=
        mapply(
            fmt_hr_ci,
            INTERACTION_HR,
            LCL95,
            UCL95,
            MoreArgs = list(
                digits = 2L
            )
        )
]

interaction_word[
    ,
    `P value` :=
        vapply(
            P,
            fmt_p,
            character(1)
        )
]

interaction_word <- interaction_word[
    ,
    .(
        Model = MODEL,
        N,
        Events = EVENTS,
        `Interaction HR (95% CI)`,
        `P value`
    )
]

contrast_word <- copy(
    models730$contrasts
)

contrast_word[
    ,
    `HR (95% CI)` :=
        mapply(
            fmt_hr_ci,
            HR,
            LCL95,
            UCL95,
            MoreArgs = list(
                digits = 2L
            )
        )
]

contrast_word[
    ,
    `P value` :=
        vapply(
            P,
            fmt_p,
            character(1)
        )
]

contrast_word <- contrast_word[
    ,
    .(
        Model = MODEL_ID,
        Contrast,
        `HR (95% CI)`,
        `P value`
    )
]

loo_word <- copy(
    loo_results
)

loo_word[
    ,
    `Interaction HR (95% CI)` :=
        ifelse(
            STATUS ==
                "SUCCESS",
            mapply(
                fmt_hr_ci,
                INTERACTION_HR,
                LCL95,
                UCL95,
                MoreArgs = list(
                    digits = 2L
                )
            ),
            NA_character_
        )
]

loo_word[
    ,
    `P value` :=
        ifelse(
            STATUS ==
                "SUCCESS",
            vapply(
                P,
                fmt_p,
                character(1)
            ),
            NA_character_
        )
]

loo_word <- loo_word[
    ,
    .(
        LOO,
        Status = STATUS,
        N,
        Events = EVENTS,
        `Interaction HR (95% CI)`,
        `P value`,
        `Max |weighted SMD|` =
            ifelse(
                is.finite(
                    MAX_ABS_WEIGHTED_SMD
                ),
                sprintf(
                    "%.4f",
                    MAX_ABS_WEIGHTED_SMD
                ),
                NA_character_
            ),
        Message = MESSAGE
    )
]

balance_word <- copy(
    ps730$balance
)

balance_word[
    ,
    `Unweighted SMD` :=
        sprintf(
            "%.3f",
            SMD_unweighted
        )
]

balance_word[
    ,
    `Overlap-weighted SMD` :=
        sprintf(
            "%.3f",
            SMD_overlap_weighted
        )
]

balance_word <- balance_word[
    ,
    .(
        Covariate,
        `Unweighted SMD`,
        `Overlap-weighted SMD`
    )
]

# ============================================================
# 11. Save RDS
# ============================================================

saveRDS(
    list(
        analysis_status =
            "POST-HOC (POST-OUTCOME) DIAGNOSTIC ONLY",
        frozen_rules = c(
            "Only shared <=730-day comparator restriction is allowed.",
            "No <=540/365 comparator sensitivity is run.",
            "PS formula is inherited unchanged from B1-14.",
            "M1/M2 formulas are inherited unchanged from B1-14B.",
            "BEV bypass+ LOO refits PS and overlap weights in every replicate.",
            "Component composition is read from outcome-blind B1-14.",
            "B1-13c specimen-timing results are read only; not re-estimated.",
            "Tier A: restricted M1 HR >=2.50 and LCL95 >1.",
            "Tier B: restricted M1 HR >1 and (HR >=1.50 or LCL95 >1), but Tier A not met.",
            "Tier C: restricted M1 HR <=1 or (HR <1.50 and LCL95 <=1).",
            "LOO influence is reported separately and does not automatically change the tier."
        ),
        hard_anchor = hard_anchor,
        original_B1_14B_M1 = original_m1,
        restricted_730 = list(
            source_N = nrow(d730_source),
            analysis_N = nrow(d730),
            arm_summary = cell730,
            ps_formula = b114$ps_formula,
            ps_model = ps730$ps_model,
            balance = ps730$balance,
            max_abs_weighted_smd =
                ps730$max_abs_smd,
            data = d730,
            models = list(
                M1 = models730$M1,
                M2 = models730$M2
            ),
            interactions =
                models730$interactions,
            contrasts =
                models730$contrasts,
            ph_diagnostics =
                models730$ph
        ),
        manuscript_framing = list(
            tier = framing_tier,
            figure4_action = figure4_action,
            title_action = title_action,
            comparator_action = comparator_action
        ),
        bev_bypass_pos_loo = list(
            planned_ids_internal = bev_pos_ids,
            results = loo_results,
            summary = loo_summary,
            interpretation = loo_interpretation
        ),
        outcome_blind_component_composition =
            component,
        existing_B1_13c_context =
            s6_existing
    ),
    out_rds
)

# ============================================================
# 12. Audit TXT
# ============================================================

audit_lines <- c(
    "B1-16 SPECIMEN-TIMING COMPARATOR STRESS TEST",
    "============================================",
    "",
    "STATUS",
    "------",
    "Post-hoc (post-outcome) diagnostic only.",
    "No new primary analysis, no new genomic exposure, no cutoff shopping.",
    "",
    "FROZEN RULES",
    "------------",
    "- Shared specimen restriction: <=730 days in BOTH treatment arms.",
    "- Refit exact B1-14 PS after restriction.",
    "- Recalculate overlap weights after PS refit.",
    "- Refit exact B1-14B M1/M2.",
    "- BEV bypass+ LOO: 14 deletions; each refits PS -> OW -> M1.",
    "- No <=540/365 comparator sensitivity.",
    "- No three-way treatment x bypass x specimen-age interaction.",
    "- No new Stage-IV-to-specimen adjustment.",
    "- No patient-level <=365 outcome excavation.",
    "",
    "HARD ANCHORS",
    "------------",
    capture.output(
        print(
            hard_anchor
        )
    ),
    "",
    "ORIGINAL FROZEN B1-14B M1",
    "--------------------------",
    capture.output(
        print(
            original_m1
        )
    ),
    "",
    "SHARED <=730-DAY POPULATION",
    "---------------------------",
    capture.output(
        print(
            cell730
        )
    ),
    paste0(
        "Maximum absolute overlap-weighted SMD after <=730 PS refit = ",
        sprintf(
            "%.6f",
            ps730$max_abs_smd
        )
    ),
    "",
    "SHARED <=730-DAY INTERACTION MODELS",
    "-----------------------------------",
    capture.output(
        print(
            models730$interactions
        )
    ),
    "",
    "SHARED <=730-DAY DERIVED CONTRASTS",
    "----------------------------------",
    capture.output(
        print(
            models730$contrasts
        )
    ),
    "",
    "SHARED <=730-DAY PH DIAGNOSTICS",
    "--------------------------------",
    capture.output(
        print(
            models730$ph
        )
    ),
    "",
    "FROZEN MANUSCRIPT FRAMING DECISION",
    "----------------------------------",
    paste0(
        "Tier: ",
        framing_tier
    ),
    paste0(
        "Figure 4: ",
        figure4_action
    ),
    paste0(
        "Title: ",
        title_action
    ),
    paste0(
        "Comparator: ",
        comparator_action
    ),
    "",
    "BEVACIZUMAB BYPASS+ LEAVE-ONE-OUT",
    "---------------------------------",
    capture.output(
        print(
            loo_results
        )
    ),
    "",
    capture.output(
        print(
            loo_summary
        )
    ),
    "",
    "OUTCOME-BLIND COMPONENT COMPOSITION",
    "-----------------------------------",
    capture.output(
        print(
            component
        )
    ),
    "",
    "PREEXISTING B1-13c CONTEXT",
    "--------------------------",
    if (
        !is.null(
            s6_existing$interaction_results
        )
    ) {
        capture.output(
            print(
                s6_existing$interaction_results
            )
        )
    } else {
        "interaction_results slot unavailable"
    },
    "",
    if (
        !is.null(
            s6_existing$included_excluded_730_summary
        )
    ) {
        capture.output(
            print(
                s6_existing$included_excluded_730_summary
            )
        )
    } else {
        "included_excluded_730_summary slot unavailable"
    },
    "",
    "INTERPRETATION GUARDRAILS",
    "-------------------------",
    "- The frozen primary anti-EGFR rwPFS result is not redefined by B1-16.",
    "- Specimen recency remains an important limitation regardless of B1-16 tier.",
    "- Overlap-weighted exact mean balance of modeled covariates is expected by construction under logistic-PS overlap weighting and is not evidence of causal exchangeability.",
    "- Common-support coverage is diagnostic; the overlap-weighted analysis is not trimmed to the empirical common-support interval.",
    "- LOO is an influence diagnostic, not a search for a favorable subset.",
    "- Component composition is outcome-blind and descriptive.",
    "- No causal treatment-selection or established predictive-biomarker conclusion is permitted.",
    "",
    paste0(
        "Saved RDS: ",
        out_rds
    ),
    paste0(
        "Saved diagnostic Word: ",
        out_doc
    )
)

writeLines(
    audit_lines,
    out_txt
)

# ============================================================
# 13. Diagnostic Word output
# ============================================================

doc <- read_docx()

doc <- body_set_default_section(
    doc,
    prop_section(
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
)

doc <- body_add_par(
    doc,
    "B1-16. Specimen-timing comparator stress test",
    style = "heading 1"
)

doc <- body_add_par(
    doc,
    paste0(
        "Classification: post-hoc (post-outcome) diagnostic only. ",
        "The frozen primary anti-EGFR analysis is unchanged. ",
        "Only the previously defined <=730-day restriction is applied to both treatment arms."
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "A. Shared <=730-day restricted comparator population"
)

doc <- body_add_flextable(
    doc,
    three_line(
        cell730,
        7.4
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "B. Shared <=730-day treatment-by-bypass interaction"
)

doc <- body_add_flextable(
    doc,
    three_line(
        interaction_word,
        7.4
    )
)

doc <- body_add_par(
    doc,
    paste0(
        "Overlap weighting used a refitted version of the unchanged B1-14 logistic propensity model. ",
        "Maximum absolute overlap-weighted SMD after refitting = ",
        sprintf(
            "%.4f",
            ps730$max_abs_smd
        ),
        ". Exact mean balance of modeled covariates under logistic-PS overlap weighting is expected by construction; ",
        "it does not establish causal exchangeability."
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "C. Contrasts derived from the restricted M1/M2 models"
)

doc <- body_add_flextable(
    doc,
    three_line(
        contrast_word,
        7.1
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "D. Frozen manuscript-framing rule applied to restricted M1"
)

doc <- body_add_flextable(
    doc,
    three_line(
        framing_table,
        7.2
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "E. Bevacizumab bypass-positive leave-one-out influence diagnostic"
)

doc <- body_add_flextable(
    doc,
    three_line(
        loo_word,
        6.8
    )
)

doc <- body_add_par(
    doc,
    loo_interpretation
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "F. Leave-one-out summary"
)

doc <- body_add_flextable(
    doc,
    three_line(
        loo_summary,
        7.2
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "G. Existing outcome-blind high-confidence component composition by treatment arm"
)

doc <- body_add_flextable(
    doc,
    three_line(
        component_word,
        7.0
    )
)

doc <- body_add_par(
    doc,
    paste0(
        "Component counts are read directly from the outcome-blind B1-14 feasibility object. ",
        "Percentages use 23 anti-EGFR and 16 bevacizumab high-confidence-positive patients as denominators. ",
        "Patients may carry more than one qualifying component, so percentages need not sum to 100%."
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    "H. Restricted-model proportional-hazards diagnostics"
)

ph_word <- copy(
    models730$ph
)

if ("P" %in% names(ph_word)) {
    ph_word[
        ,
        P :=
            ifelse(
                is.finite(P),
                vapply(
                    P,
                    fmt_p,
                    character(1)
                ),
                NA_character_
            )
    ]
}

doc <- body_add_flextable(
    doc,
    three_line(
        ph_word,
        7.0
    )
)

doc <- body_add_par(
    doc,
    ""
)

doc <- body_add_par(
    doc,
    paste0(
        "Guardrail: B1-16 does not select a new PH remedy, does not introduce a new cutoff, ",
        "and does not redefine the primary exposure or primary endpoint. ",
        "Specimen-recency sensitivity remains a limitation regardless of the B1-16 framing tier."
    )
)

print(
    doc,
    target = out_doc
)

# ============================================================
# 14. Console summary
# ============================================================

cat("\n============================================================\n")
cat("B1-16 SPECIMEN-TIMING COMPARATOR STRESS TEST COMPLETE\n")
cat("============================================================\n\n")

cat("HARD ANCHORS\n")
cat("------------\n")
print(hard_anchor)

cat("\nSHARED <=730-DAY POPULATION\n")
cat("---------------------------\n")
print(cell730)

cat("\nSHARED <=730-DAY INTERACTION MODELS\n")
cat("-----------------------------------\n")
print(models730$interactions)

cat("\nFROZEN FRAMING DECISION\n")
cat("-----------------------\n")
cat(framing_tier, "\n")
cat("Figure 4: ", figure4_action, "\n", sep = "")
cat("Title: ", title_action, "\n", sep = "")
cat("Comparator: ", comparator_action, "\n", sep = "")

cat("\nBEV BYPASS+ LOO SUMMARY\n")
cat("-----------------------\n")
print(loo_summary)

cat("\nB1-16 RULE REMINDER\n")
cat("-------------------\n")
cat(
    "Do not add <=540/365 comparator models, non-prespecified subgroups, ",
    "three-way interactions, or outcome-driven model changes after seeing this result.\n",
    sep = ""
)

cat("\nSaved audit:\n", out_txt, "\n", sep = "")
cat("Saved RDS:\n", out_rds, "\n", sep = "")
cat("Saved diagnostic Word:\n", out_doc, "\n", sep = "")
