# ============================================================
# A7 / B-line
# 08_B1_08_robustness_validation.R
#
# POST-UNBLINDING ROBUSTNESS VALIDATION
#
# IMPORTANT:
#   These analyses are explicitly POST-UNBLINDING robustness checks.
#   They do NOT redefine the primary genomic exposure.
#
# Checks:
#   1) Patient-level bootstrap stability of the primary adjusted rwPFS HR
#   2) Leave-one-component-out high-confidence composite
#   3) Exclude carriers of each individual component
#   4) 14-day and 30-day landmark sensitivity for very early events
#   5) Descriptive RMST differences at 180 and 365 days
#
# Input:
#   03_intermediate/A7_B1_locked_analysis.rds
#
# Outputs:
#   04_results/B1_08_robustness_validation.rds
#   06_logs_and_audit/B1_08_robustness_validation.txt
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

locked_file <- file.path(
    project_root,
    "03_intermediate",
    "A7_B1_locked_analysis.rds"
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

if (!file.exists(locked_file)) {
    stop(
        "Missing locked B1 dataset:\n",
        locked_file
    )
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
cat("B1-08 ROBUSTNESS VALIDATION - STANDALONE RMST FIX\n")
cat("============================================================\n\n")

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
        "Locked analysis anchors failed."
    )
}


# ============================================================
# Helpers
# ============================================================

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

fmt_num <- function(
    x,
    digits = 2
) {

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

extract_exposure <- function(
    fit,
    term,
    analysis,
    component = NA_character_
) {

    s <- summary(fit)

    if (
        !term %in%
            rownames(
                s$coefficients
            )
    ) {

        return(
            data.table(
                ANALYSIS = analysis,
                COMPONENT = component,
                N = fit$n,
                EVENTS = fit$nevent,
                BYPASS_POS = NA_integer_,
                HR = NA_real_,
                LCL95 = NA_real_,
                UCL95 = NA_real_,
                P = NA_real_
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
        ANALYSIS = analysis,
        COMPONENT = component,
        N = fit$n,
        EVENTS = fit$nevent,
        BYPASS_POS = NA_integer_,
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

fit_adjusted <- function(
    data,
    exposure_col
) {

    form <- as.formula(
        paste0(
            "Surv(RWPFS_DAYS_LOCKED_R, RWPFS_EVENT_LOCKED_R) ~ ",
            exposure_col,
            " + GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
            "LOG_STAGE4_TO_ANTI_EGFR_R + PRIOR_HEAVY_F_R"
        )
    )

    coxph(
        form,
        data = data,
        ties = "efron"
    )
}

rmst_from_survfit <- function(
    sf,
    tau
) {

    tau <- as.numeric(tau)

    if (
        !is.finite(tau) ||
        tau <= 0
    ) {
        stop("tau must be positive.")
    }

    # One-stratum Kaplan-Meier fit.
    if (
        length(sf$time) == 0L
    ) {
        return(tau)
    }

    tm <- as.numeric(sf$time)
    sv <- as.numeric(sf$surv)

    ok <- is.finite(tm) & is.finite(sv)

    tm <- tm[ok]
    sv <- sv[ok]

    if (
        length(tm) == 0L
    ) {
        return(tau)
    }

    ord <- order(tm)
    tm <- tm[ord]
    sv <- sv[ord]

    # survfit is usually unique by time; if duplicate times occur,
    # retain the last survival estimate at that time.
    u <- unique(tm)

    s_u <- vapply(
        u,
        function(z) {
            tail(
                sv[tm == z],
                1L
            )
        },
        numeric(1)
    )

    # Exact area under KM step function from 0 to tau.
    change_times <- u[
        u > 0 &
        u < tau
    ]

    left <- c(
        0,
        change_times
    )

    right <- c(
        change_times,
        tau
    )

    surv_left <- vapply(
        left,
        function(z) {

            if (z == 0) {
                return(1)
            }

            j <- max(
                which(
                    u <= z
                )
            )

            s_u[j]
        },
        numeric(1)
    )

    area <- sum(
        (right - left) *
        surv_left
    )

    as.numeric(area)
}

rmst_by_group <- function(
    data,
    tau
) {

    d0 <- data[
        BYPASS_HIGH_CONFIDENCE_R == 0L
    ]

    d1 <- data[
        BYPASS_HIGH_CONFIDENCE_R == 1L
    ]

    sf0 <- survfit(
        Surv(
            RWPFS_DAYS_LOCKED_R,
            RWPFS_EVENT_LOCKED_R
        ) ~ 1,
        data = d0
    )

    sf1 <- survfit(
        Surv(
            RWPFS_DAYS_LOCKED_R,
            RWPFS_EVENT_LOCKED_R
        ) ~ 1,
        data = d1
    )

    r0 <- rmst_from_survfit(
        sf0,
        tau
    )

    r1 <- rmst_from_survfit(
        sf1,
        tau
    )

    data.table(
        TAU_DAYS = tau,
        RMST_BYPASS_NEG_DAYS = r0,
        RMST_BYPASS_POS_DAYS = r1,
        DIFFERENCE_POS_MINUS_NEG_DAYS =
            r1 - r0
    )
}


# ============================================================
# RMST function self-check
# ============================================================

test_sf <- survfit(
    Surv(
        c(10, 20, 30),
        c(1, 1, 0)
    ) ~ 1
)

test_rmst <- rmst_from_survfit(
    test_sf,
    25
)

if (
    !is.finite(test_rmst) ||
    test_rmst <= 0 ||
    test_rmst > 25
) {
    stop(
        "RMST helper self-check failed."
    )
}


# ============================================================
# 1. Reproduce primary adjusted rwPFS model
# ============================================================

primary_fit <- fit_adjusted(
    dt,
    "BYPASS_HIGH_CONFIDENCE_R"
)

primary_result <- extract_exposure(
    primary_fit,
    "BYPASS_HIGH_CONFIDENCE_R",
    "Primary adjusted rwPFS"
)

primary_result[
    ,
    BYPASS_POS := dt[
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        .N
    ]
]


# ============================================================
# 2. Patient-level bootstrap
# ============================================================

set.seed(
    20260831
)

B <- 2000L

boot_beta <- rep(
    NA_real_,
    B
)

boot_hr <- rep(
    NA_real_,
    B
)

for (b in seq_len(B)) {

    idx <- sample.int(
        nrow(dt),
        size = nrow(dt),
        replace = TRUE
    )

    db <- copy(
        dt[idx]
    )

    # Skip pathological resamples.
    if (
        uniqueN(
            db$BYPASS_HIGH_CONFIDENCE_R
        ) < 2L
    ) {
        next
    }

    fit_b <- try(
        fit_adjusted(
            db,
            "BYPASS_HIGH_CONFIDENCE_R"
        ),
        silent = TRUE
    )

    if (
        inherits(
            fit_b,
            "try-error"
        )
    ) {
        next
    }

    cf <- coef(
        fit_b
    )

    if (
        !"BYPASS_HIGH_CONFIDENCE_R" %in%
            names(cf)
    ) {
        next
    }

    beta <- unname(
        cf[
            "BYPASS_HIGH_CONFIDENCE_R"
        ]
    )

    if (
        !is.finite(
            beta
        )
    ) {
        next
    }

    boot_beta[b] <- beta
    boot_hr[b] <- exp(beta)
}

boot_hr_valid <- boot_hr[
    is.finite(
        boot_hr
    )
]

boot_beta_valid <- boot_beta[
    is.finite(
        boot_beta
    )
]

if (
    length(
        boot_hr_valid
    ) <
        0.90 *
            B
) {
    warning(
        "Fewer than 90% of bootstrap models converged."
    )
}

bootstrap_summary <- data.table(
    REQUESTED_ITERATIONS = B,
    SUCCESSFUL_ITERATIONS =
        length(
            boot_hr_valid
        ),

    PRIMARY_HR =
        primary_result$HR,

    BOOTSTRAP_HR_MEDIAN =
        median(
            boot_hr_valid
        ),

    BOOTSTRAP_HR_MEAN =
        mean(
            boot_hr_valid
        ),

    BOOTSTRAP_HR_P2_5 =
        unname(
            quantile(
                boot_hr_valid,
                0.025,
                na.rm = TRUE
            )
        ),

    BOOTSTRAP_HR_P97_5 =
        unname(
            quantile(
                boot_hr_valid,
                0.975,
                na.rm = TRUE
            )
        ),

    PROP_HR_GT_1 =
        mean(
            boot_hr_valid >
                1
        ),

    PROP_HR_GT_1_5 =
        mean(
            boot_hr_valid >
                1.5
        ),

    BOOTSTRAP_BETA_SD =
        sd(
            boot_beta_valid
        )
)


# ============================================================
# 3. Leave-one-component-out composite
# ============================================================

components <- c(
    "ERBB2_AMP_R",
    "ERBB2_ACTIVATING_R",
    "MET_AMP_R",
    "PIK3CA_EX20_R",
    "PTEN_LOF_R",
    "AKT1_E17K_R",
    "CANONICAL_PRESSING_FUSION_R",
    "KRAS_AMP_R",
    "MAP2K1_HOTSPOT_R",
    "NF1_LOF_R"
)

component_counts <- data.table(
    COMPONENT = components,
    N_CARRIERS = vapply(
        components,
        function(v) {
            sum(
                dt[[v]] == 1L,
                na.rm = TRUE
            )
        },
        integer(1)
    )
)

nonzero_components <- component_counts[
    N_CARRIERS > 0L,
    COMPONENT
]

loo_results <- rbindlist(
    lapply(
        nonzero_components,
        function(remove_component) {

            keep_components <- setdiff(
                components,
                remove_component
            )

            x <- rowSums(
                as.data.frame(
                    dt[
                        ,
                        ..keep_components
                    ]
                ),
                na.rm = TRUE
            )

            d <- copy(dt)

            d[
                ,
                BYPASS_LOO_R :=
                    as.integer(
                        x > 0
                    )
            ]

            if (
                d[
                    BYPASS_LOO_R == 1L,
                    .N
                ] <
                    10L
            ) {

                return(
                    data.table(
                        ANALYSIS =
                            "Leave-one-component-out",
                        COMPONENT =
                            remove_component,
                        N =
                            nrow(d),
                        EVENTS =
                            d[
                                RWPFS_EVENT_LOCKED_R == 1L,
                                .N
                            ],
                        BYPASS_POS =
                            d[
                                BYPASS_LOO_R == 1L,
                                .N
                            ],
                        HR = NA_real_,
                        LCL95 = NA_real_,
                        UCL95 = NA_real_,
                        P = NA_real_
                    )
                )
            }

            fit <- fit_adjusted(
                d,
                "BYPASS_LOO_R"
            )

            out <- extract_exposure(
                fit,
                "BYPASS_LOO_R",
                "Leave-one-component-out",
                remove_component
            )

            out[
                ,
                BYPASS_POS :=
                    d[
                        BYPASS_LOO_R == 1L,
                        .N
                    ]
            ]

            out
        }
    )
)


# ============================================================
# 4. Exclude carriers of each component
# ============================================================

carrier_exclusion_results <- rbindlist(
    lapply(
        nonzero_components,
        function(component) {

            d <- dt[
                get(component) == 0L
            ]

            bypass_pos <- d[
                BYPASS_HIGH_CONFIDENCE_R == 1L,
                .N
            ]

            if (
                bypass_pos <
                    10L
            ) {

                return(
                    data.table(
                        ANALYSIS =
                            "Exclude component carriers",
                        COMPONENT =
                            component,
                        N =
                            nrow(d),
                        EVENTS =
                            d[
                                RWPFS_EVENT_LOCKED_R == 1L,
                                .N
                            ],
                        BYPASS_POS =
                            bypass_pos,
                        HR = NA_real_,
                        LCL95 = NA_real_,
                        UCL95 = NA_real_,
                        P = NA_real_
                    )
                )
            }

            fit <- fit_adjusted(
                d,
                "BYPASS_HIGH_CONFIDENCE_R"
            )

            out <- extract_exposure(
                fit,
                "BYPASS_HIGH_CONFIDENCE_R",
                "Exclude component carriers",
                component
            )

            out[
                ,
                BYPASS_POS :=
                    bypass_pos
            ]

            out
        }
    )
)


# ============================================================
# 5. Early-event landmark sensitivity
#
# Explicitly post hoc.
# Retain only patients alive/progression-free beyond landmark,
# reset the risk clock to the landmark.
# ============================================================

landmark_fit <- function(
    landmark_days
) {

    d <- dt[
        RWPFS_DAYS_LOCKED_R >
            landmark_days
    ]

    d[
        ,
        RWPFS_LM_DAYS_R :=
            RWPFS_DAYS_LOCKED_R -
            landmark_days
    ]

    fit <- coxph(
        Surv(
            RWPFS_LM_DAYS_R,
            RWPFS_EVENT_LOCKED_R
        ) ~
            BYPASS_HIGH_CONFIDENCE_R +
            GENDER_F_R +
            SUBSITE_F_R +
            ANTI_EGFR_AGENT_F_R +
            LOG_STAGE4_TO_ANTI_EGFR_R +
            PRIOR_HEAVY_F_R,
        data = d,
        ties = "efron"
    )

    out <- extract_exposure(
        fit,
        "BYPASS_HIGH_CONFIDENCE_R",
        paste0(
            landmark_days,
            "-day landmark"
        )
    )

    out[
        ,
        BYPASS_POS :=
            d[
                BYPASS_HIGH_CONFIDENCE_R == 1L,
                .N
            ]
    ]

    out
}

landmark_results <- rbindlist(
    list(
        landmark_fit(
            14
        ),
        landmark_fit(
            30
        )
    )
)


# ============================================================
# 6. Descriptive RMST
# ============================================================

rmst_results <- rbindlist(
    list(
        rmst_by_group(
            dt,
            180
        ),

        rmst_by_group(
            dt,
            365
        )
    )
)


# ============================================================
# 7. Save outputs
# ============================================================

result_file <- file.path(
    results_dir,
    "B1_08_robustness_validation.rds"
)

audit_file <- file.path(
    audit_dir,
    "B1_08_robustness_validation.txt"
)

saveRDS(
    list(
        primary_result =
            primary_result,

        bootstrap_summary =
            bootstrap_summary,

        component_counts =
            component_counts,

        leave_one_component_out =
            loo_results,

        carrier_exclusion =
            carrier_exclusion_results,

        landmark =
            landmark_results,

        rmst =
            rmst_results,

        bootstrap_hr =
            boot_hr_valid
    ),
    result_file
)

lines <- c(
    "B1-08 POST-UNBLINDING ROBUSTNESS VALIDATION",
    "===========================================",
    "",
    "IMPORTANT",
    "---------",
    "These are explicitly post-unblinding robustness analyses.",
    "They do NOT redefine the frozen primary exposure.",
    "",
    "1. Reproduced primary adjusted rwPFS result",
    "-------------------------------------------",
    capture.output(
        print(
            primary_result
        )
    ),
    "",
    "2. Patient-level bootstrap",
    "--------------------------",
    capture.output(
        print(
            bootstrap_summary
        )
    ),
    "",
    "3. High-confidence component counts",
    "-----------------------------------",
    capture.output(
        print(
            component_counts
        )
    ),
    "",
    "4. Leave-one-component-out composite",
    "------------------------------------",
    capture.output(
        print(
            loo_results
        )
    ),
    "",
    "5. Exclude carriers of each component",
    "-------------------------------------",
    capture.output(
        print(
            carrier_exclusion_results
        )
    ),
    "",
    "6. Early-event landmark sensitivity",
    "-----------------------------------",
    capture.output(
        print(
            landmark_results
        )
    ),
    "",
    "7. Descriptive rwPFS RMST",
    "------------------------",
    capture.output(
        print(
            rmst_results
        )
    ),
    "",
    "Interpretation guardrails",
    "-------------------------",
    "- Do not present component-specific p-values as gene-level efficacy claims.",
    "- Leave-one-component-out results assess whether one component dominates the composite.",
    "- Landmark analyses are post hoc data-quality sensitivities, not replacement primary analyses.",
    "- Bootstrap results assess stability, not external validity.",
    "",
    "Next gate:",
    "If the HR remains directionally stable across component-removal, bootstrap and landmark checks, proceed to an independent endpoint validation (TTNTD) before manuscript figures."
)

writeLines(
    lines,
    audit_file
)

cat("\n============================================================\n")
cat("B1-08 COMPLETE\n")
cat("============================================================\n")

cat("\nPrimary adjusted rwPFS:\n")
print(
    primary_result
)

cat("\nBootstrap summary:\n")
print(
    bootstrap_summary
)

cat("\nLeave-one-component-out:\n")
print(
    loo_results
)

cat("\nAudit:\n", audit_file, "\n")
