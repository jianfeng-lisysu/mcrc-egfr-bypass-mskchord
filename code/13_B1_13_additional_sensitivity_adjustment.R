# ============================================================
# A7 / B-line
# 13_B1_13_additional_sensitivity_adjustment.R
#
# POST-HOC SENSITIVITY ANALYSIS
#
# PURPOSE
#   Probe residual confounding without
#   changing the locked genomic exposure, primary endpoint, or the
#   original primary analysis.
#
# IMPORTANT ANALYSIS HIERARCHY
#   - Model A is the original LOCKED primary adjusted rwPFS model.
#   - Models B-E are post-hoc sensitivity models.
#   - ALL Models A-E are reported regardless of direction/significance.
#   - No model is selected after seeing the result.
#   - Model E is the strongest pressure test and sets the ceiling for
#     manuscript wording.
#   - These analyses remain ASSOCIATIONAL, not predictive/causal.
#
# PREDECLARED B1-13 WORDING GATE
#   This gate is declared in this script BEFORE B1-13 is run.
#
#   Model E HR >= 1.90 AND lower 95% CI > 1:
#       strong post-hoc robustness;
#       stronger associative wording remains defensible.
#
#   Model E HR >= 1.70 AND lower 95% CI > 1:
#       attenuated but robust;
#       use moderate wording ("associated with shorter rwPFS"),
#       avoid "markedly" / resistance-like language.
#
#   Otherwise:
#       material attenuation / uncertainty;
#       more cautious manuscript wording is required.
#
# MODEL DEFINITIONS
#
# A. Original locked primary model
#    bypass + sex + subsite + anti-EGFR agent +
#    log(Stage IV -> anti-EGFR interval) +
#    prior both oxaliplatin + irinotecan
#
# B. Treatment-context model
#    Model A core covariates, but replace PRIOR_HEAVY with an honest
#    0/1/2 prior-major-cytotoxic-class proxy and add concurrent
#    irinotecan / oxaliplatin backbone context at the anti-EGFR index.
#
# C. Expanded clinical/genomic-context model
#    Model B + panel generation + sample type +
#    specimen-to-index interval + recent metastatic-site count.
#
# D. ECOG complete-case model
#    Model C + nearest pre-index ECOG within 90 days (0 vs >=1).
#    Always attempted and always reported. If not estimable because of
#    sparse complete cases, the non-estimability is explicitly reported.
#
# E. Strongest pressure test
#    Model C restricted to specimens acquired <=730 days before
#    anti-EGFR index.
#
# ADDITIONAL REQUIRED AUDITS
#   1) bypass+ distribution across prior major cytotoxic proxy 0/1/2.
#   2) within-stratum rwPFS estimates for each 0/1/2 stratum:
#      - unadjusted
#      - core-adjusted if estimable
#      Wide CIs are reported rather than hidden.
#   3) PH diagnostics for Models A-E when estimable.
#
# NOT DONE HERE
#   - No anti-EGFR vs bevacizumab interaction.
#   - No redefinition of the high-confidence bypass exposure.
#   - No broad-only / neither / high-confidence three-group analysis.
#     That analysis, if pursued, belongs in a separate B1-13b and must
#     be labeled hypothesis-generating and underpowered.
#
# INPUTS
#   03_intermediate/A7_B1_locked_analysis.rds
#   01_raw_data/msk_chord_2024.tar
#
# OUTPUTS (OVERWRITE)
#   04_results/B1_13_additional_sensitivity_adjustment.rds
#   06_logs_and_audit/B1_13_additional_sensitivity_adjustment.txt
#   07_tables/TableS5_additional_sensitivity_analyses.docx
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

table_dir <- file.path(
    project_root,
    "07_tables"
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
    locked_file,
    tar_file
)) {
    if (!file.exists(f0)) {
        stop(
            "Missing required input:\n",
            f0
        )
    }
}

required_packages <- c(
    "data.table",
    "survival"
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

cat("\n============================================================\n")
cat("B1-13 ADDITIONAL SENSITIVITY ADJUSTMENT\n")
cat("============================================================\n\n")


# ============================================================
# 1. Helpers
# ============================================================

safe_chr <- function(x) {
    z <- as.character(x)
    z[is.na(z)] <- ""
    z
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

    z
}

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

has_variation <- function(
    d,
    v
) {

    if (!v %in% names(d)) {
        return(FALSE)
    }

    x <- d[[v]]

    if (is.factor(x)) {
        x <- droplevels(x)
    }

    length(
        unique(
            x[
                !is.na(x)
            ]
        )
    ) >= 2L
}

extract_exposure <- function(
    fit,
    model_id,
    model_label,
    status = "OK"
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
                PH_BYPASS_P = NA_real_,
                PH_GLOBAL_P = NA_real_,
                STATUS =
                    paste0(
                        status,
                        "; exposure term not estimable"
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

    ph_obj <- tryCatch(
        cox.zph(
            fit,
            transform = "km"
        ),
        error = function(e) NULL
    )

    if (!is.null(ph_obj)) {

        ph_tab <- as.data.frame(
            ph_obj$table
        )

        if (
            term %in%
                rownames(
                    ph_tab
                )
        ) {
            ph_bypass <-
                ph_tab[
                    term,
                    "p"
                ]
        }

        if (
            "GLOBAL" %in%
                rownames(
                    ph_tab
                )
        ) {
            ph_global <-
                ph_tab[
                    "GLOBAL",
                    "p"
                ]
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
        PH_BYPASS_P = ph_bypass,
        PH_GLOBAL_P = ph_global,
        STATUS = status
    )
}

fit_safe <- function(
    formula,
    data,
    model_id,
    model_label
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
        error = function(e) {
            attr(
                NULL,
                "error_message"
            ) <- conditionMessage(e)
            NULL
        }
    )

    if (is.null(fit)) {

        out <- extract_exposure(
            NULL,
            model_id,
            model_label,
            status =
                "NOT ESTIMABLE / Cox model failed"
        )

        out[
            ,
            BYPASS_POS :=
                bypass_pos
        ]

        return(
            list(
                fit = NULL,
                result = out
            )
        )
    }

    out <- extract_exposure(
        fit,
        model_id,
        model_label
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


# ============================================================
# 2. Read locked cohort
# ============================================================

d <- as.data.table(
    readRDS(
        locked_file
    )
)

if (
    nrow(d) != 191L ||
    d[
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        .N
    ] != 23L
) {
    stop(
        "Locked cohort hard anchors failed: expected N=191 and bypass+=23."
    )
}

required_locked_cols <- c(
    "PATIENT_ID",
    "BYPASS_HIGH_CONFIDENCE_R",
    "FIRST_BIOLOGIC_DAY_R",
    "FOLLOWUP_DAY_R",
    "DEATH_AFTER_INDEX_R",
    "DEATH_DAY_R",
    "GENDER_F_R",
    "SUBSITE_F_R",
    "ANTI_EGFR_AGENT_F_R",
    "LOG_STAGE4_TO_ANTI_EGFR_R",
    "PRIOR_HEAVY_F_R",
    "PANEL_OLD_R",
    "SAMPLE_TYPE",
    "SPECIMEN_TO_ANTI_EGFR_DAYS_R",
    "MET_SITE_COUNT_90D_R"
)

missing_locked_cols <- setdiff(
    required_locked_cols,
    names(d)
)

if (
    length(
        missing_locked_cols
    ) > 0L
) {
    stop(
        "Locked RDS is missing required columns:\n",
        paste(
            missing_locked_cols,
            collapse = "\n"
        )
    )
}


# ============================================================
# 3. Extract raw timelines
# ============================================================

tmp_dir <- tempfile(
    "A7_B1_13_"
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
        ) == filename
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

prog <- fread(
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

ps <- fread(
    find_raw(
        "data_timeline_performance_status.txt"
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
        d$PATIENT_ID
]

prog <- prog[
    PATIENT_ID %chin%
        d$PATIENT_ID
]

ps <- ps[
    PATIENT_ID %chin%
        d$PATIENT_ID
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

prog[
    ,
    PROG_DAY_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

ps[
    ,
    PS_DAY_R :=
        suppressWarnings(
            as.numeric(
                START_DATE
            )
        )
]

ps[
    ,
    ECOG_NUM_R :=
        suppressWarnings(
            as.numeric(
                ECOG
            )
        )
]


# ============================================================
# 4. Reconstruct rwPFS exactly as B1-07
# ============================================================

prog_y <- prog[
    toupper(
        trimws(
            safe_chr(
                PROGRESSION
            )
        )
    ) == "Y" &
        is.finite(
            PROG_DAY_R
        )
]

first_prog <- rbindlist(
    lapply(
        seq_len(
            nrow(d)
        ),
        function(i) {

            pid <- d$PATIENT_ID[i]
            idx <- d$FIRST_BIOLOGIC_DAY_R[i]
            fup <- d$FOLLOWUP_DAY_R[i]

            cand <- prog_y[
                PATIENT_ID == pid &
                    PROG_DAY_R > idx &
                    PROG_DAY_R <= fup
            ]

            day <- if (
                nrow(cand) > 0L
            ) {
                min(
                    cand$PROG_DAY_R,
                    na.rm = TRUE
                )
            } else {
                NA_real_
            }

            data.table(
                PATIENT_ID = pid,
                FIRST_PROG_DAY_B13_R = day
            )
        }
    ),
    fill = TRUE
)

d <- merge(
    d,
    first_prog,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

d[
    ,
    DEATH_EVENT_DAY_B13_R :=
        fifelse(
            DEATH_AFTER_INDEX_R == 1L &
                is.finite(
                    DEATH_DAY_R
                ) &
                DEATH_DAY_R >
                    FIRST_BIOLOGIC_DAY_R &
                DEATH_DAY_R <=
                    FOLLOWUP_DAY_R,
            DEATH_DAY_R,
            NA_real_
        )
]

d[
    ,
    RWPFS_EVENT_DAY_B13_R :=
        pmin(
            FIRST_PROG_DAY_B13_R,
            DEATH_EVENT_DAY_B13_R,
            na.rm = TRUE
        )
]

d[
    !is.finite(
        RWPFS_EVENT_DAY_B13_R
    ),
    RWPFS_EVENT_DAY_B13_R :=
        NA_real_
]

d[
    ,
    RWPFS_EVENT_B13_R :=
        as.integer(
            is.finite(
                RWPFS_EVENT_DAY_B13_R
            )
        )
]

d[
    ,
    RWPFS_END_DAY_B13_R :=
        fifelse(
            RWPFS_EVENT_B13_R == 1L,
            RWPFS_EVENT_DAY_B13_R,
            FOLLOWUP_DAY_R
        )
]

d[
    ,
    RWPFS_DAYS_B13_R :=
        RWPFS_END_DAY_B13_R -
            FIRST_BIOLOGIC_DAY_R
]

if (
    d[
        !is.finite(
            RWPFS_DAYS_B13_R
        ) |
            RWPFS_DAYS_B13_R <= 0,
        .N
    ] > 0L
) {
    stop(
        "Invalid reconstructed rwPFS durations."
    )
}

if (
    sum(
        d$RWPFS_EVENT_B13_R
    ) != 183L
) {
    stop(
        "rwPFS reconstruction hard audit failed: expected 183 events, observed ",
        sum(
            d$RWPFS_EVENT_B13_R
        ),
        "."
    )
}


# ============================================================
# 5. Build 0/1/2 prior-major-cytotoxic proxy
# ============================================================

# Prefer an outcome-blind, previously constructed variable if present.
if (
    "PRIOR_CYTOTOXIC_PATTERN_R" %in%
        names(d)
) {

    z <- tolower(
        safe_chr(
            d$PRIOR_CYTOTOXIC_PATTERN_R
        )
    )

    d[
        ,
        PRIOR_MAJOR_N_B13_R :=
            fcase(
                grepl(
                    "oxaliplatin.*irinotecan|irinotecan.*oxaliplatin",
                    z
                ),
                2L,

                grepl(
                    "one major",
                    z
                ),
                1L,

                grepl(
                    "no prior major|fluoropyrimidine only",
                    z
                ),
                0L,

                default =
                    NA_integer_
            )
    ]

    prior_source <-
        "Existing outcome-blind PRIOR_CYTOTOXIC_PATTERN_R"

} else {

    # Fallback: reconstruct prior exposure to the two major cytotoxic
    # classes from raw treatment starts strictly before anti-EGFR index.
    prior_dt <- rbindlist(
        lapply(
            seq_len(
                nrow(d)
            ),
            function(i) {

                pid <- d$PATIENT_ID[i]
                idx <- d$FIRST_BIOLOGIC_DAY_R[i]

                p <- tx[
                    PATIENT_ID == pid &
                        is.finite(
                            START_DAY_R
                        ) &
                        START_DAY_R < idx
                ]

                prior_ox <- as.integer(
                    any(
                        p$AGENT_KEY_R ==
                            "OXALIPLATIN"
                    )
                )

                prior_iri <- as.integer(
                    any(
                        p$AGENT_KEY_R ==
                            "IRINOTECAN"
                    )
                )

                data.table(
                    PATIENT_ID = pid,
                    PRIOR_OX_B13_R = prior_ox,
                    PRIOR_IRI_B13_R = prior_iri,
                    PRIOR_MAJOR_N_B13_R =
                        prior_ox +
                            prior_iri
                )
            }
        ),
        fill = TRUE
    )

    d <- merge(
        d,
        prior_dt,
        by = "PATIENT_ID",
        all.x = TRUE,
        sort = FALSE
    )

    prior_source <-
        "Reconstructed from raw treatment timeline: prior oxaliplatin/irinotecan exposure before index"
}

if (
    any(
        !d$PRIOR_MAJOR_N_B13_R %in%
            0:2
    )
) {
    stop(
        "Prior-major-cytotoxic proxy contains missing/invalid values."
    )
}

d[
    ,
    PRIOR_MAJOR_N_F_B13_R :=
        factor(
            PRIOR_MAJOR_N_B13_R,
            levels = 0:2,
            labels = c(
                "0",
                "1",
                "2"
            )
        )
]

# Integrity cross-check against the original binary PRIOR_HEAVY variable.
prior_heavy_num <- suppressWarnings(
    as.integer(
        as.character(
            d$PRIOR_HEAVY_F_R
        )
    )
)

if (
    all(
        is.na(
            prior_heavy_num
        )
    )
) {

    prior_heavy_num <-
        as.integer(
            tolower(
                as.character(
                    d$PRIOR_HEAVY_F_R
                )
            ) %chin%
                c(
                    "1",
                    "yes",
                    "true"
                )
        )
}

prior_two_num <- as.integer(
    d$PRIOR_MAJOR_N_B13_R == 2L
)

prior_heavy_agreement <- mean(
    prior_two_num ==
        prior_heavy_num,
    na.rm = TRUE
)

if (
    is.finite(
        prior_heavy_agreement
    ) &&
    prior_heavy_agreement < 0.95
) {
    stop(
        "Prior-treatment proxy disagrees materially with PRIOR_HEAVY_F_R. Agreement = ",
        sprintf(
            "%.3f",
            prior_heavy_agreement
        ),
        ". Review before modeling."
    )
}


# ============================================================
# 6. Concurrent backbone at anti-EGFR index
# ============================================================

# Prefer the pre-existing outcome-blind index-backbone field when present.
if (
    "ACTIVE_BACKBONE_AT_INDEX_R" %in%
        names(d)
) {

    btxt <- toupper(
        safe_chr(
            d$ACTIVE_BACKBONE_AT_INDEX_R
        )
    )

    d[
        ,
        ACTIVE_IRI_B13_R :=
            as.integer(
                grepl(
                    "IRINOTECAN|FOLFIRI|CAPIRI|FOLFOXIRI",
                    btxt
                )
            )
    ]

    d[
        ,
        ACTIVE_OX_B13_R :=
            as.integer(
                grepl(
                    "OXALIPLATIN|FOLFOX|CAPOX|FOLFOXIRI",
                    btxt
                )
            )
    ]

    backbone_source <-
        "Existing outcome-blind ACTIVE_BACKBONE_AT_INDEX_R"

} else {

    # Conservative fallback using only treatment information available
    # on/before the anti-EGFR index:
    #   start within 30 days before/same day OR earlier and still active.
    active_dt <- rbindlist(
        lapply(
            seq_len(
                nrow(d)
            ),
            function(i) {

                pid <- d$PATIENT_ID[i]
                idx <- d$FIRST_BIOLOGIC_DAY_R[i]

                p <- tx[
                    PATIENT_ID == pid &
                        (
                            (
                                START_DAY_R >=
                                    idx - 30 &
                                START_DAY_R <=
                                    idx
                            ) |
                            (
                                START_DAY_R <
                                    idx - 30 &
                                is.finite(
                                    STOP_DAY_R
                                ) &
                                STOP_DAY_R >=
                                    idx
                            )
                        )
                ]

                data.table(
                    PATIENT_ID = pid,

                    ACTIVE_IRI_B13_R =
                        as.integer(
                            any(
                                p$AGENT_KEY_R ==
                                    "IRINOTECAN"
                            )
                        ),

                    ACTIVE_OX_B13_R =
                        as.integer(
                            any(
                                p$AGENT_KEY_R ==
                                    "OXALIPLATIN"
                            )
                        )
                )
            }
        ),
        fill = TRUE
    )

    d <- merge(
        d,
        active_dt,
        by = "PATIENT_ID",
        all.x = TRUE,
        sort = FALSE
    )

    backbone_source <-
        "Fallback reconstructed from on/before-index treatment timeline; no post-index drugs used"
}

d[
    ,
    ACTIVE_IRI_F_B13_R :=
        factor(
            ACTIVE_IRI_B13_R,
            levels = c(
                0,
                1
            ),
            labels = c(
                "No",
                "Yes"
            )
        )
]

d[
    ,
    ACTIVE_OX_F_B13_R :=
        factor(
            ACTIVE_OX_B13_R,
            levels = c(
                0,
                1
            ),
            labels = c(
                "No",
                "Yes"
            )
        )
]


# ============================================================
# 7. Expanded context covariates
# ============================================================

d[
    ,
    PANEL_OLD_F_B13_R :=
        factor(
            PANEL_OLD_R,
            levels = c(
                0,
                1
            ),
            labels = c(
                "Newer panel",
                "IMPACT341/410"
            )
        )
]

sample_txt <- toupper(
    trimws(
        safe_chr(
            d$SAMPLE_TYPE
        )
    )
)

d[
    ,
    SAMPLE_METASTATIC_B13_R :=
        fifelse(
            !nzchar(
                sample_txt
            ),
            NA_integer_,
            as.integer(
                grepl(
                    "METAST",
                    sample_txt
                )
            )
        )
]

d[
    ,
    SAMPLE_METASTATIC_F_B13_R :=
        factor(
            SAMPLE_METASTATIC_B13_R,
            levels = c(
                0,
                1
            ),
            labels = c(
                "Other/nonmetastatic",
                "Metastatic"
            )
        )
]

d[
    ,
    SPECIMEN_AGE_DAYS_B13_R :=
        suppressWarnings(
            as.numeric(
                SPECIMEN_TO_ANTI_EGFR_DAYS_R
            )
        )
]

d[
    ,
    LOG_SPECIMEN_AGE_B13_R :=
        fifelse(
            is.finite(
                SPECIMEN_AGE_DAYS_B13_R
            ) &
                SPECIMEN_AGE_DAYS_B13_R >= 0,
            log1p(
                SPECIMEN_AGE_DAYS_B13_R
            ),
            NA_real_
        )
]

d[
    ,
    MET_SITE_COUNT_B13_R :=
        suppressWarnings(
            as.numeric(
                MET_SITE_COUNT_90D_R
            )
        )
]

d[
    ,
    SPECIMEN_WITHIN_730_B13_R :=
        as.integer(
            is.finite(
                SPECIMEN_AGE_DAYS_B13_R
            ) &
                SPECIMEN_AGE_DAYS_B13_R >= 0 &
                SPECIMEN_AGE_DAYS_B13_R <= 730
        )
]


# ============================================================
# 8. ECOG complete-case variable
# ============================================================

ps <- ps[
    is.finite(
        PS_DAY_R
    ) &
        is.finite(
            ECOG_NUM_R
        ) &
        ECOG_NUM_R >= 0 &
        ECOG_NUM_R <= 4
]

ecog_nearest <- rbindlist(
    lapply(
        seq_len(
            nrow(d)
        ),
        function(i) {

            pid <- d$PATIENT_ID[i]
            idx <- d$FIRST_BIOLOGIC_DAY_R[i]

            p <- ps[
                PATIENT_ID == pid &
                    PS_DAY_R <= idx &
                    PS_DAY_R >=
                        idx - 90
            ]

            if (
                nrow(p) == 0L
            ) {
                return(
                    data.table(
                        PATIENT_ID = pid,
                        ECOG_B13_R = NA_real_,
                        ECOG_DAY_B13_R = NA_real_,
                        ECOG_DAYS_BEFORE_INDEX_B13_R =
                            NA_real_
                    )
                )
            }

            # Nearest pre-index record. If several records occur on the
            # same nearest day, conservatively use the highest ECOG.
            max_day <- max(
                p$PS_DAY_R,
                na.rm = TRUE
            )

            pp <- p[
                PS_DAY_R == max_day
            ]

            ecog_val <- max(
                pp$ECOG_NUM_R,
                na.rm = TRUE
            )

            data.table(
                PATIENT_ID = pid,
                ECOG_B13_R = ecog_val,
                ECOG_DAY_B13_R = max_day,
                ECOG_DAYS_BEFORE_INDEX_B13_R =
                    idx - max_day
            )
        }
    ),
    fill = TRUE
)

d <- merge(
    d,
    ecog_nearest,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

d[
    ,
    ECOG_0_VS_GE1_B13_R :=
        fcase(
            ECOG_B13_R == 0,
            "0",

            is.finite(
                ECOG_B13_R
            ) &
                ECOG_B13_R >= 1,
            ">=1",

            default =
                NA_character_
        )
]

d[
    ,
    ECOG_F_B13_R :=
        factor(
            ECOG_0_VS_GE1_B13_R,
            levels = c(
                "0",
                ">=1"
            )
        )
]


# ============================================================
# 9. Model formulas A-E
# ============================================================

surv_rhs <-
    "Surv(RWPFS_DAYS_B13_R, RWPFS_EVENT_B13_R)"

formula_A <- as.formula(
    paste0(
        surv_rhs,
        " ~ BYPASS_HIGH_CONFIDENCE_R + ",
        "GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
        "LOG_STAGE4_TO_ANTI_EGFR_R + PRIOR_HEAVY_F_R"
    )
)

formula_B <- as.formula(
    paste0(
        surv_rhs,
        " ~ BYPASS_HIGH_CONFIDENCE_R + ",
        "GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
        "LOG_STAGE4_TO_ANTI_EGFR_R + ",
        "PRIOR_MAJOR_N_F_B13_R + ",
        "ACTIVE_IRI_F_B13_R + ACTIVE_OX_F_B13_R"
    )
)

formula_C <- as.formula(
    paste0(
        surv_rhs,
        " ~ BYPASS_HIGH_CONFIDENCE_R + ",
        "GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
        "LOG_STAGE4_TO_ANTI_EGFR_R + ",
        "PRIOR_MAJOR_N_F_B13_R + ",
        "ACTIVE_IRI_F_B13_R + ACTIVE_OX_F_B13_R + ",
        "PANEL_OLD_F_B13_R + SAMPLE_METASTATIC_F_B13_R + ",
        "LOG_SPECIMEN_AGE_B13_R + MET_SITE_COUNT_B13_R"
    )
)

formula_D <- as.formula(
    paste0(
        surv_rhs,
        " ~ BYPASS_HIGH_CONFIDENCE_R + ",
        "GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + ",
        "LOG_STAGE4_TO_ANTI_EGFR_R + ",
        "PRIOR_MAJOR_N_F_B13_R + ",
        "ACTIVE_IRI_F_B13_R + ACTIVE_OX_F_B13_R + ",
        "PANEL_OLD_F_B13_R + SAMPLE_METASTATIC_F_B13_R + ",
        "LOG_SPECIMEN_AGE_B13_R + MET_SITE_COUNT_B13_R + ",
        "ECOG_F_B13_R"
    )
)

formula_E <- formula_C


# ============================================================
# 10. Fit A-E — ALL ARE REPORTED
# ============================================================

mA <- fit_safe(
    formula_A,
    d,
    "A",
    "Original locked primary adjusted model"
)

# HARD REPRODUCTION GATE:
# Model A must reproduce the original B1-07 HR = 2.334372.
if (
    is.na(
        mA$result$HR
    ) ||
    abs(
        mA$result$HR -
            2.334372
    ) > 0.01
) {
    stop(
        "Model A failed to reproduce locked B1-07 HR=2.334372. Observed HR=",
        paste(
            mA$result$HR,
            collapse = ","
        ),
        ". Stop before interpreting B1-13."
    )
}

mB <- fit_safe(
    formula_B,
    d,
    "B",
    "Treatment-context expanded adjustment"
)

mC <- fit_safe(
    formula_C,
    d,
    "C",
    "Expanded treatment + clinical/genomic-context adjustment"
)

mD <- fit_safe(
    formula_D,
    d,
    "D",
    "Model C + nearest pre-index ECOG within 90 days"
)

d_E <- d[
    SPECIMEN_WITHIN_730_B13_R == 1L
]

mE <- fit_safe(
    formula_E,
    d_E,
    "E",
    "Model C restricted to specimen <=730 days before anti-EGFR"
)

model_results <- rbindlist(
    list(
        mA$result,
        mB$result,
        mC$result,
        mD$result,
        mE$result
    ),
    fill = TRUE
)


# ============================================================
# 11. Required treatment-depth distribution audit
# ============================================================

prior_depth_distribution <- d[
    ,
    .(
        N = .N,
        EVENTS =
            sum(
                RWPFS_EVENT_B13_R
            ),
        BYPASS_POS =
            sum(
                BYPASS_HIGH_CONFIDENCE_R == 1L
            ),
        BYPASS_NEG =
            sum(
                BYPASS_HIGH_CONFIDENCE_R == 0L
            )
    ),
    by =
        PRIOR_MAJOR_N_F_B13_R
][
    ,
    BYPASS_POS_PCT :=
        100 *
            BYPASS_POS /
            N
]

setnames(
    prior_depth_distribution,
    "PRIOR_MAJOR_N_F_B13_R",
    "PRIOR_MAJOR_N"
)


# ============================================================
# 12. Within-stratum estimates for prior-major 0/1/2
# ============================================================

fit_stratum_one <- function(
    k
) {

    ds <- d[
        PRIOR_MAJOR_N_B13_R == k
    ]

    n_pos <- ds[
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        .N
    ]

    n_neg <- ds[
        BYPASS_HIGH_CONFIDENCE_R == 0L,
        .N
    ]

    ev_pos <- ds[
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        sum(
            RWPFS_EVENT_B13_R
        )
    ]

    ev_neg <- ds[
        BYPASS_HIGH_CONFIDENCE_R == 0L,
        sum(
            RWPFS_EVENT_B13_R
        )
    ]

    base_meta <- data.table(
        PRIOR_MAJOR_N = as.character(k),
        N = nrow(ds),
        EVENTS =
            sum(
                ds$RWPFS_EVENT_B13_R
            ),
        BYPASS_POS = n_pos,
        BYPASS_NEG = n_neg,
        BYPASS_POS_EVENTS = ev_pos,
        BYPASS_NEG_EVENTS = ev_neg
    )

    if (
        n_pos == 0L ||
        n_neg == 0L
    ) {

        return(
            rbindlist(
                list(
                    cbind(
                        base_meta,
                        data.table(
                            ANALYSIS =
                                "Within-stratum unadjusted",
                            HR = NA_real_,
                            LCL95 = NA_real_,
                            UCL95 = NA_real_,
                            P = NA_real_,
                            STATUS =
                                "NOT ESTIMABLE: one exposure group absent"
                        )
                    ),

                    cbind(
                        base_meta,
                        data.table(
                            ANALYSIS =
                                "Within-stratum core-adjusted",
                            HR = NA_real_,
                            LCL95 = NA_real_,
                            UCL95 = NA_real_,
                            P = NA_real_,
                            STATUS =
                                "NOT ESTIMABLE: one exposure group absent"
                        )
                    )
                ),
                fill = TRUE
            )
        )
    }

    f_unadj <- as.formula(
        paste0(
            surv_rhs,
            " ~ BYPASS_HIGH_CONFIDENCE_R"
        )
    )

    u <- fit_safe(
        f_unadj,
        ds,
        paste0(
            "STRATUM_",
            k,
            "_UNADJ"
        ),
        paste0(
            "Prior-major ",
            k,
            ": unadjusted"
        )
    )$result

    # Core-adjusted stratum model.
    # Covariates are removed ONLY if they have no variation in that
    # stratum. This is a structural/estimability rule, not an
    # outcome-based selection rule.
    candidate_covars <- c(
        "GENDER_F_R",
        "SUBSITE_F_R",
        "ANTI_EGFR_AGENT_F_R",
        "LOG_STAGE4_TO_ANTI_EGFR_R",
        "ACTIVE_IRI_F_B13_R",
        "ACTIVE_OX_F_B13_R"
    )

    varying_covars <- candidate_covars[
        vapply(
            candidate_covars,
            function(v) {
                has_variation(
                    ds,
                    v
                )
            },
            logical(1)
        )
    ]

    f_adj <- as.formula(
        paste0(
            surv_rhs,
            " ~ BYPASS_HIGH_CONFIDENCE_R",
            if (
                length(
                    varying_covars
                ) > 0L
            ) {
                paste0(
                    " + ",
                    paste(
                        varying_covars,
                        collapse = " + "
                    )
                )
            } else {
                ""
            }
        )
    )

    a <- fit_safe(
        f_adj,
        ds,
        paste0(
            "STRATUM_",
            k,
            "_ADJ"
        ),
        paste0(
            "Prior-major ",
            k,
            ": core-adjusted"
        )
    )$result

    rbindlist(
        list(
            cbind(
                base_meta,
                data.table(
                    ANALYSIS =
                        "Within-stratum unadjusted",
                    HR = u$HR,
                    LCL95 = u$LCL95,
                    UCL95 = u$UCL95,
                    P = u$P,
                    STATUS = u$STATUS
                )
            ),

            cbind(
                base_meta,
                data.table(
                    ANALYSIS =
                        "Within-stratum core-adjusted",
                    HR = a$HR,
                    LCL95 = a$LCL95,
                    UCL95 = a$UCL95,
                    P = a$P,
                    STATUS = a$STATUS
                )
            )
        ),
        fill = TRUE
    )
}

stratum_results <- rbindlist(
    lapply(
        0:2,
        fit_stratum_one
    ),
    fill = TRUE
)


# ============================================================
# 13. ECOG feasibility audit
# ============================================================

ecog_feasibility <- d[
    ,
    .(
        N = .N,
        ECOG_AVAILABLE =
            sum(
                !is.na(
                    ECOG_F_B13_R
                )
            ),
        ECOG_0 =
            sum(
                ECOG_0_VS_GE1_B13_R == "0",
                na.rm = TRUE
            ),
        ECOG_GE1 =
            sum(
                ECOG_0_VS_GE1_B13_R == ">=1",
                na.rm = TRUE
            )
    ),
    by =
        BYPASS_HIGH_CONFIDENCE_R
]

ecog_feasibility[
    ,
    ECOG_AVAILABLE_PCT :=
        100 *
            ECOG_AVAILABLE /
            N
]


# ============================================================
# 14. Context-variable distributions
# ============================================================

backbone_distribution <- d[
    ,
    .(
        N = .N,
        BYPASS_POS =
            sum(
                BYPASS_HIGH_CONFIDENCE_R == 1L
            )
    ),
    by = .(
        ACTIVE_IRI =
            ACTIVE_IRI_F_B13_R,
        ACTIVE_OX =
            ACTIVE_OX_F_B13_R
    )
]

model_e_counts <- d_E[
    ,
    .(
        N = .N,
        EVENTS =
            sum(
                RWPFS_EVENT_B13_R
            ),
        BYPASS_POS =
            sum(
                BYPASS_HIGH_CONFIDENCE_R == 1L
            ),
        BYPASS_NEG =
            sum(
                BYPASS_HIGH_CONFIDENCE_R == 0L
            )
    )
]


# ============================================================
# 15. Model-E wording gate
# ============================================================

e_row <- model_results[
    MODEL_ID == "E"
]

if (
    nrow(e_row) != 1L ||
    !is.finite(
        e_row$HR
    ) ||
    !is.finite(
        e_row$LCL95
    )
) {

    wording_gate <-
        paste0(
            "MODEL E NOT RELIABLY ESTIMABLE. ",
            "More cautious manuscript wording is required until the reason is resolved; ",
            "do not use Model B/C as a substitute."
        )

    wording_category <-
        "NOT ESTIMABLE / DOWNGRADE"

} else if (
    e_row$HR >= 1.90 &&
    e_row$LCL95 > 1
) {

    wording_gate <-
        paste0(
            "Model E met the strongest predeclared B1-13 robustness gate ",
            "(HR >=1.90 and lower 95% CI >1). ",
            "Strong ASSOCIATIVE wording remains defensible, but predictive, ",
            "causal, treatment-benefit, or gene-specific claims remain prohibited."
        )

    wording_category <-
        "STRONG POST-HOC ROBUSTNESS"

} else if (
    e_row$HR >= 1.70 &&
    e_row$LCL95 > 1
) {

    wording_gate <-
        paste0(
            "Model E shows attenuation but retains a robust association ",
            "(HR >=1.70 and lower 95% CI >1). ",
            "Use moderate wording such as 'associated with shorter rwPFS'; ",
            "avoid 'markedly poorer', 'resistance biomarker', or equivalent language."
        )

    wording_category <-
        "ATTENUATED BUT ROBUST"

} else {

    wording_gate <-
        paste0(
            "Model E shows material attenuation and/or substantial uncertainty. ",
            "The manuscript wording ceiling must be lowered regardless of stronger ",
            "Model A-C estimates."
        )

    wording_category <-
        "MATERIAL ATTENUATION / UNCERTAINTY"
}


# ============================================================
# 16. Save RDS
# ============================================================

result_file <- file.path(
    results_dir,
    "B1_13_additional_sensitivity_adjustment.rds"
)

saveRDS(
    list(
        patient_data = d,
        model_results = model_results,
        models = list(
            A = mA$fit,
            B = mB$fit,
            C = mC$fit,
            D = mD$fit,
            E = mE$fit
        ),
        formulas = list(
            A = formula_A,
            B = formula_B,
            C = formula_C,
            D = formula_D,
            E = formula_E
        ),
        prior_depth_distribution =
            prior_depth_distribution,
        within_prior_depth_strata =
            stratum_results,
        ecog_feasibility =
            ecog_feasibility,
        backbone_distribution =
            backbone_distribution,
        model_e_counts =
            model_e_counts,
        prior_proxy_source =
            prior_source,
        backbone_source =
            backbone_source,
        prior_heavy_agreement =
            prior_heavy_agreement,
        wording_category =
            wording_category,
        wording_gate =
            wording_gate
    ),
    result_file
)


# ============================================================
# 17. Audit TXT
# ============================================================

audit_file <- file.path(
    audit_dir,
    "B1_13_additional_sensitivity_adjustment.txt"
)

audit_lines <- c(
    "B1-13 ADDITIONAL SENSITIVITY ADJUSTMENT",
    "===================================",
    "",
    "ANALYSIS STATUS",
    "---------------",
    "POST-HOC sensitivity analysis.",
    "Frozen high-confidence genomic exposure is unchanged.",
    "Primary rwPFS endpoint is unchanged.",
    "Original Model A remains the primary adjusted analysis.",
    "Models B-E do NOT replace Model A.",
    "ALL Models A-E are reported regardless of result.",
    "Model E sets the ceiling for manuscript wording.",
    "",
    "B1-13b is intentionally NOT run here.",
    "Any broad-only / neither / high-confidence three-group analysis must be",
    "a separate, explicitly underpowered, hypothesis-generating analysis.",
    "",
    "LOCKED ANCHORS",
    "--------------",
    paste0(
        "N = ",
        nrow(d)
    ),
    paste0(
        "Bypass+ = ",
        d[
            BYPASS_HIGH_CONFIDENCE_R == 1L,
            .N
        ]
    ),
    paste0(
        "Bypass- = ",
        d[
            BYPASS_HIGH_CONFIDENCE_R == 0L,
            .N
        ]
    ),
    paste0(
        "rwPFS events = ",
        sum(
            d$RWPFS_EVENT_B13_R
        )
    ),
    paste0(
        "Model A reproduced HR = ",
        sprintf(
            "%.6f",
            mA$result$HR
        ),
        " [target 2.334372]"
    ),
    "",
    "VARIABLE CONSTRUCTION",
    "---------------------",
    paste0(
        "Prior-major 0/1/2 proxy source: ",
        prior_source
    ),
    paste0(
        "Agreement of prior-major=2 with locked PRIOR_HEAVY: ",
        sprintf(
            "%.1f%%",
            100 *
                prior_heavy_agreement
        )
    ),
    paste0(
        "Concurrent backbone source: ",
        backbone_source
    ),
    "ECOG: nearest pre-index ECOG within 90 days; same-day ties use highest ECOG.",
    "",
    "MODELS A-E",
    "----------",
    capture.output(
        print(
            model_results
        )
    ),
    "",
    "MODEL DEFINITIONS",
    "-----------------",
    "A = original locked primary adjusted model.",
    "B = 0/1/2 prior-major proxy + concurrent irinotecan/oxaliplatin context.",
    "C = B + panel + sample type + specimen age + recent metastatic-site count.",
    "D = C + nearest pre-index ECOG within 90 days, complete-case.",
    "E = C restricted to specimen <=730 days before anti-EGFR.",
    "",
    "PRIOR MAJOR CYTOTOXIC 0/1/2 DISTRIBUTION",
    "-----------------------------------------",
    capture.output(
        print(
            prior_depth_distribution
        )
    ),
    "",
    "WITHIN-STRATUM ESTIMATES",
    "------------------------",
    "Both unadjusted and core-adjusted estimates are attempted for every 0/1/2 stratum.",
    "Sparse/unstable strata remain visible; they are not deleted after seeing results.",
    capture.output(
        print(
            stratum_results
        )
    ),
    "",
    "CONCURRENT BACKBONE DISTRIBUTION",
    "--------------------------------",
    capture.output(
        print(
            backbone_distribution
        )
    ),
    "",
    "ECOG FEASIBILITY",
    "----------------",
    capture.output(
        print(
            ecog_feasibility
        )
    ),
    "",
    "MODEL E COUNTS",
    "--------------",
    capture.output(
        print(
            model_e_counts
        )
    ),
    "",
    "PREDECLARED B1-13 WORDING GATE",
    "------------------------------",
    "E HR >=1.90 and lower95>1: strong post-hoc robustness.",
    "E HR >=1.70 and lower95>1: attenuated but robust; moderate wording only.",
    "Otherwise: material attenuation/uncertainty; use more cautious wording.",
    "",
    paste0(
        "RESULT CATEGORY: ",
        wording_category
    ),
    wording_gate,
    "",
    "INTERPRETATION GUARDRAILS",
    "-------------------------",
    "- This remains an association study within anti-EGFR-treated patients.",
    "- Do NOT call the result predictive without a valid comparator interaction.",
    "- Do NOT redefine the genomic exposure after B1-13.",
    "- Do NOT hide Model E if it is weaker than Models A-C.",
    "- Within-stratum estimates may be imprecise and are treatment-depth diagnostics,",
    "  not independent replication.",
    "- Bootstrap / multiple sensitivity analyses do not substitute for external validation.",
    "",
    paste0(
        "Saved RDS: ",
        result_file
    )
)

writeLines(
    audit_lines,
    audit_file
)


# ============================================================
# 18. Optional manuscript-facing Word table
# ============================================================

word_file <- file.path(
    table_dir,
    "TableS5_additional_sensitivity_analyses.docx"
)

word_status <- "Not attempted"

if (
    requireNamespace(
        "officer",
        quietly = TRUE
    ) ||
    tryCatch(
        {
            install.packages(
                "officer"
            )
            TRUE
        },
        error = function(e) FALSE
    )
) {

    if (
        requireNamespace(
            "flextable",
            quietly = TRUE
        ) ||
        tryCatch(
            {
                install.packages(
                    "flextable"
                )
                TRUE
            },
            error = function(e) FALSE
        )
    ) {

        word_status <- tryCatch(
            {

                library(officer)
                library(flextable)

                model_table_word <- model_results[
                    ,
                    .(
                        Model = MODEL_ID,
                        Description = MODEL,
                        N,
                        Events = EVENTS,
                        `Bypass+, n` = BYPASS_POS,
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

                stratum_table_word <- stratum_results[
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
                            mapply(
                                fmt_hr,
                                HR,
                                LCL95,
                                UCL95
                            ),
                        Status = STATUS
                    )
                ]

                ft_three_line <- function(
                    df,
                    font_size = 8
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
                            top = 0.55,
                            bottom = 0.55,
                            left = 0.55,
                            right = 0.55
                        )
                )

                doc <- body_set_default_section(
                    doc,
                    sec
                )

                doc <- body_add_par(
                    doc,
                    "Supplementary Table S5. Additional treatment-context sensitivity analyses",
                    style = "Normal"
                )

                doc <- body_add_par(
                    doc,
                    "A-E models",
                    style = "Normal"
                )

                doc <- body_add_flextable(
                    doc,
                    ft_three_line(
                        model_table_word,
                        7.5
                    )
                )

                doc <- body_add_par(
                    doc,
                    ""
                )

                doc <- body_add_par(
                    doc,
                    "Within-stratum estimates by prior major cytotoxic exposure count",
                    style = "Normal"
                )

                doc <- body_add_flextable(
                    doc,
                    ft_three_line(
                        stratum_table_word,
                        7.5
                    )
                )

                doc <- body_add_par(
                    doc,
                    ""
                )

                doc <- body_add_par(
                    doc,
                    paste0(
                        "Notes: Model A is the original locked primary adjusted model. ",
                        "Models B-E are post-hoc sensitivity analyses and ",
                        "do not replace Model A. Model E is the strongest pressure test. ",
                        "The 0/1/2 variable is an honest proxy for prior exposure to the two ",
                        "major cytotoxic classes (oxaliplatin and irinotecan), not a formal ",
                        "line-of-therapy number. ECOG Model D uses the nearest ECOG recorded ",
                        "within 90 days before the anti-EGFR index. All A-E models are shown ",
                        "regardless of result."
                    ),
                    style = "Normal"
                )

                print(
                    doc,
                    target = word_file
                )

                "Created successfully"
            },
            error = function(e) {
                paste0(
                    "Word table failed, but RDS/TXT remain valid: ",
                    conditionMessage(e)
                )
            }
        )
    }
}


# ============================================================
# 19. Console summary
# ============================================================

cat("\n============================================================\n")
cat("B1-13 COMPLETE\n")
cat("============================================================\n\n")

print(
    model_results
)

cat("\nPrior-major 0/1/2 distribution:\n")
print(
    prior_depth_distribution
)

cat("\nWithin-stratum estimates:\n")
print(
    stratum_results
)

cat("\nECOG feasibility:\n")
print(
    ecog_feasibility
)

cat("\nWORDING GATE:\n")
cat(
    wording_category,
    "\n"
)
cat(
    wording_gate,
    "\n"
)

cat(
    "\nRDS:\n",
    result_file,
    "\n"
)

cat(
    "\nAudit:\n",
    audit_file,
    "\n"
)

cat(
    "\nWord table:\n",
    word_file,
    "\nStatus: ",
    word_status,
    "\n"
)
