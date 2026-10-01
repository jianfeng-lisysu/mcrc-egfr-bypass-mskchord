# ============================================================
# A7 / B-line
# 14C_B1_14_full_pipeline_bootstrap.R
#
# Post-outcome robustness only.
# Frozen pipeline:
# patient resample -> refit B1-14 PS -> recalculate overlap weights
# -> refit unchanged B1-14B M1 treatment*bypass weighted Cox.
#
# B = 2000; patient-level bootstrap; fixed seed.
# M1 remains the primary exploratory comparator model.
# This bootstrap does not establish causality or a predictive biomarker.
# ============================================================

options(stringsAsFactors = FALSE)

root <- "/path/to/A7"  # set to your local project directory

f14  <- file.path(root, "04_results", "B1_14_bevacizumab_comparator_feasibility.rds")
f14b <- file.path(root, "04_results", "B1_14B_exploratory_comparator_rwPFS.rds")

out_rds <- file.path(root, "04_results", "B1_14C_full_pipeline_bootstrap.rds")
out_txt <- file.path(root, "06_logs_and_audit", "B1_14C_full_pipeline_bootstrap.txt")
out_doc <- file.path(root, "07_tables", "TableS9_full_pipeline_bootstrap.docx")

for (p in c("data.table", "survival", "officer", "flextable")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(data.table)
library(survival)
library(officer)
library(flextable)

B <- 2000L
SEED <- 20260904L

term_t <- "ANTI_EGFR_TREAT_R"
term_b <- "BYPASS_HIGH_CONFIDENCE_R"
term_i <- "ANTI_EGFR_TREAT_R:BYPASS_HIGH_CONFIDENCE_R"

fmt_hr <- function(hr, lo, hi) {
    if (any(!is.finite(c(hr, lo, hi)))) return("NA")
    sprintf("%.2f (%.2f–%.2f)", hr, lo, hi)
}

q3 <- function(x) {
    z <- x[is.finite(x)]
    if (!length(z)) return(c(NA_real_, NA_real_, NA_real_))
    as.numeric(quantile(z, c(.025, .5, .975), names = FALSE))
}

# ------------------------------------------------------------
# 1. Load frozen objects and anchors
# ------------------------------------------------------------

b14 <- readRDS(f14)
b14b <- readRDS(f14b)

if (!identical(as.character(b14$gate_decision), "GO")) {
    stop("B1-14 gate is not GO.")
}

d <- as.data.table(copy(b14b$analysis_data))
m1_saved <- as.data.table(b14b$model_contrasts)[MODEL_ID == "M1"]

if (
    nrow(d) != 319L ||
    d[FIRST_BIOLOGIC_R == "anti-EGFR", .N] != 188L ||
    d[FIRST_BIOLOGIC_R == "bevacizumab", .N] != 131L ||
    d[FIRST_BIOLOGIC_R == "anti-EGFR" & BYPASS_HIGH_CONFIDENCE_R == 1L, .N] != 23L ||
    d[FIRST_BIOLOGIC_R == "bevacizumab" & BYPASS_HIGH_CONFIDENCE_R == 1L, .N] != 14L ||
    sum(d$RWPFS_EVENT_B14B_R) != 305L ||
    nrow(m1_saved) != 1L
) {
    stop("B1-14B frozen analysis anchors failed.")
}

obs_hr <- m1_saved$INTERACTION_HR
obs_lo <- m1_saved$INTERACTION_LCL95
obs_hi <- m1_saved$INTERACTION_UCL95
obs_p  <- m1_saved$INTERACTION_P

if (abs(obs_hr - 3.083038) > .001) {
    stop("Unexpected saved M1 interaction HR.")
}

ps_formula <- b14$ps_formula

# ------------------------------------------------------------
# 2. Reproduce original PS + M1 point estimate before bootstrap
# ------------------------------------------------------------

ps0 <- glm(ps_formula, data = d, family = binomial())

if (!isTRUE(ps0$converged) || any(!is.finite(coef(ps0)))) {
    stop("Frozen PS reproduction failed.")
}

d[, PS_B14C := as.numeric(predict(ps0, type = "response"))]
d[, OW_B14C := fifelse(ANTI_EGFR_TREAT_R == 1L, 1 - PS_B14C, PS_B14C)]

fit0 <- coxph(
    Surv(RWPFS_DAYS_B14B_R, RWPFS_EVENT_B14B_R) ~
        ANTI_EGFR_TREAT_R * BYPASS_HIGH_CONFIDENCE_R,
    data = d,
    weights = OW_B14C,
    ties = "efron"
)

refit_hr <- exp(coef(fit0)[term_i])

if (!is.finite(refit_hr) || abs(refit_hr - obs_hr) > 1e-5) {
    stop(
        sprintf(
            "Frozen M1 reproduction failed: saved %.8f vs refit %.8f",
            obs_hr,
            refit_hr
        )
    )
}

cat("Frozen M1 reproduction: PASS\n")
cat(sprintf("Observed interaction HR = %.4f\n\n", obs_hr))

# ------------------------------------------------------------
# 3. One bootstrap replicate
# ------------------------------------------------------------

boot_one <- function(idx, id) {

    x <- copy(d[idx])

    # Require all four treatment x bypass cells.
    cells <- x[, .N, by = .(ANTI_EGFR_TREAT_R, BYPASS_HIGH_CONFIDENCE_R)]
    if (nrow(cells) < 4L || min(cells$N) < 1L) {
        return(data.table(
            REP = id, STATUS = "FAIL_ZERO_CELL",
            INTERACTION_HR = NA_real_,
            TREAT_NEG = NA_real_,
            TREAT_POS = NA_real_,
            BYPASS_BEV = NA_real_,
            BYPASS_ANTI = NA_real_
        ))
    }

    ps <- tryCatch(
        suppressWarnings(glm(ps_formula, data = x, family = binomial())),
        error = function(e) NULL
    )

    if (is.null(ps) || !isTRUE(ps$converged) || any(!is.finite(coef(ps)))) {
        return(data.table(
            REP = id, STATUS = "FAIL_PS",
            INTERACTION_HR = NA_real_,
            TREAT_NEG = NA_real_,
            TREAT_POS = NA_real_,
            BYPASS_BEV = NA_real_,
            BYPASS_ANTI = NA_real_
        ))
    }

    pr <- tryCatch(
        as.numeric(predict(ps, type = "response")),
        error = function(e) NULL
    )

    if (is.null(pr) || length(pr) != nrow(x) || any(!is.finite(pr))) {
        return(data.table(
            REP = id, STATUS = "FAIL_PS_PRED",
            INTERACTION_HR = NA_real_,
            TREAT_NEG = NA_real_,
            TREAT_POS = NA_real_,
            BYPASS_BEV = NA_real_,
            BYPASS_ANTI = NA_real_
        ))
    }

    x[, W := fifelse(ANTI_EGFR_TREAT_R == 1L, 1 - pr, pr)]

    fit <- tryCatch(
        suppressWarnings(
            coxph(
                Surv(RWPFS_DAYS_B14B_R, RWPFS_EVENT_B14B_R) ~
                    ANTI_EGFR_TREAT_R * BYPASS_HIGH_CONFIDENCE_R,
                data = x,
                weights = W,
                ties = "efron"
            )
        ),
        error = function(e) NULL
    )

    if (is.null(fit)) {
        return(data.table(
            REP = id, STATUS = "FAIL_COX",
            INTERACTION_HR = NA_real_,
            TREAT_NEG = NA_real_,
            TREAT_POS = NA_real_,
            BYPASS_BEV = NA_real_,
            BYPASS_ANTI = NA_real_
        ))
    }

    cf <- coef(fit)

    if (
        !all(c(term_t, term_b, term_i) %in% names(cf)) ||
        any(!is.finite(cf[c(term_t, term_b, term_i)]))
    ) {
        return(data.table(
            REP = id, STATUS = "FAIL_COEF",
            INTERACTION_HR = NA_real_,
            TREAT_NEG = NA_real_,
            TREAT_POS = NA_real_,
            BYPASS_BEV = NA_real_,
            BYPASS_ANTI = NA_real_
        ))
    }

    bt <- unname(cf[term_t])
    bb <- unname(cf[term_b])
    bi <- unname(cf[term_i])

    data.table(
        REP = id,
        STATUS = "SUCCESS",
        INTERACTION_HR = exp(bi),
        TREAT_NEG = exp(bt),
        TREAT_POS = exp(bt + bi),
        BYPASS_BEV = exp(bb),
        BYPASS_ANTI = exp(bb + bi)
    )
}

# ------------------------------------------------------------
# 4. Run 2,000 deterministic patient-level bootstrap replicates
# ------------------------------------------------------------

set.seed(SEED)

res <- vector("list", B)

for (b in seq_len(B)) {

    idx <- sample.int(nrow(d), nrow(d), replace = TRUE)

    res[[b]] <- boot_one(idx, b)

    if (b %% 100L == 0L) {
        cat("Bootstrap ", b, "/", B, "\n", sep = "")
    }
}

boot <- rbindlist(res)
ok <- boot[STATUS == "SUCCESS" & is.finite(INTERACTION_HR)]

n_ok <- nrow(ok)
rate_ok <- n_ok / B

qi <- q3(ok$INTERACTION_HR)
qn <- q3(ok$TREAT_NEG)
qp <- q3(ok$TREAT_POS)
qb <- q3(ok$BYPASS_BEV)
qa <- q3(ok$BYPASS_ANTI)

p_gt1  <- mean(ok$INTERACTION_HR > 1)
p_gt15 <- mean(ok$INTERACTION_HR > 1.5)
p_gt2  <- mean(ok$INTERACTION_HR > 2)

if (rate_ok < .95) {
    gate <- "NUMERICAL BOOTSTRAP INSTABILITY"
    wording <- paste0(
        "Fewer than 95% of bootstrap replicates were successful. ",
        "Do not rely on the percentile interval until failures are diagnosed."
    )
} else if (is.finite(qi[1]) && qi[1] > 1) {
    gate <- "FULL-PIPELINE STABILITY SUPPORT"
    wording <- paste0(
        "The full estimation-pipeline bootstrap supports stability of the ",
        "exploratory treatment-specific heterogeneity signal; the percentile ",
        "95% interval remained above 1."
    )
} else {
    gate <- "FULL-PIPELINE CI INCLUDES NULL"
    wording <- paste0(
        "The point estimate remained directionally stable, but the full ",
        "estimation-pipeline percentile 95% interval included 1."
    )
}

primary <- data.table(
    B = B,
    Successful = n_ok,
    Success_rate = rate_ok,
    Observed_HR = obs_hr,
    Observed_robust_LCL95 = obs_lo,
    Observed_robust_UCL95 = obs_hi,
    Bootstrap_median_HR = qi[2],
    Bootstrap_LCL95 = qi[1],
    Bootstrap_UCL95 = qi[3],
    Prop_HR_gt_1 = p_gt1,
    Prop_HR_gt_1_5 = p_gt15,
    Prop_HR_gt_2 = p_gt2
)

contrasts <- data.table(
    Contrast = c(
        "Interaction HR",
        "Anti-EGFR vs bevacizumab HR in bypass-",
        "Anti-EGFR vs bevacizumab HR in bypass+",
        "Bypass HR within bevacizumab",
        "Bypass HR within anti-EGFR"
    ),
    Observed = c(
        m1_saved$INTERACTION_HR,
        m1_saved$ANTI_VS_BEV_BYPASS_NEG_HR,
        m1_saved$ANTI_VS_BEV_BYPASS_POS_HR,
        m1_saved$BYPASS_HR_IN_BEV,
        m1_saved$BYPASS_HR_IN_ANTI
    ),
    Bootstrap_median = c(qi[2], qn[2], qp[2], qb[2], qa[2]),
    Bootstrap_LCL95  = c(qi[1], qn[1], qp[1], qb[1], qa[1]),
    Bootstrap_UCL95  = c(qi[3], qn[3], qp[3], qb[3], qa[3])
)

failures <- boot[, .N, by = STATUS][order(-N)]

# ------------------------------------------------------------
# 5. Save RDS
# ------------------------------------------------------------

saveRDS(
    list(
        status = "POST-OUTCOME FULL ESTIMATION-PIPELINE BOOTSTRAP",
        B = B,
        seed = SEED,
        original_M1 = m1_saved,
        refit_interaction_HR = refit_hr,
        primary = primary,
        contrasts = contrasts,
        failure_summary = failures,
        bootstrap_results = boot,
        bootstrap_gate = gate,
        bootstrap_wording = wording
    ),
    out_rds
)

# ------------------------------------------------------------
# 6. Audit TXT
# ------------------------------------------------------------

txt <- c(
    "B1-14C FULL ESTIMATION-PIPELINE BOOTSTRAP",
    "=========================================",
    "",
    "STATUS",
    "------",
    "Post-outcome robustness analysis only.",
    "M1 remains the primary exploratory comparator model.",
    "No exposure, endpoint, PS covariate, weight formula, or interaction model was changed.",
    "",
    "DESIGN",
    "------",
    paste0("Replicates = ", B),
    paste0("Seed = ", SEED),
    paste0("Bootstrap population N = ", nrow(d)),
    "Each replicate: patient resample -> refit frozen PS -> recalculate overlap weights -> refit frozen M1 weighted Cox.",
    "",
    "ORIGINAL M1 REPRODUCTION",
    "------------------------",
    paste0("Saved interaction HR = ", sprintf("%.8f", obs_hr)),
    paste0("Refit interaction HR = ", sprintf("%.8f", refit_hr)),
    "Reproduction: PASS",
    "",
    "BOOTSTRAP STATUS",
    "----------------",
    capture.output(print(failures)),
    paste0(
        "Successful = ", n_ok, "/", B,
        " (", sprintf("%.1f%%", 100 * rate_ok), ")"
    ),
    "",
    "PRIMARY INTERACTION",
    "-------------------",
    capture.output(print(primary)),
    "",
    "DERIVED M1 CONTRASTS",
    "--------------------",
    capture.output(print(contrasts)),
    "",
    "INTERPRETATION",
    "--------------",
    paste0("BOOTSTRAP GATE: ", gate),
    wording,
    "",
    "GUARDRAILS",
    "----------",
    "- This bootstrap propagates PS-estimation and overlap-weighting uncertainty through M1.",
    "- It does not remove residual or unmeasured confounding.",
    "- It does not establish causal treatment selection.",
    "- It does not establish a predictive biomarker.",
    "- Independent non-MSK replication remains necessary.",
    "",
    paste0("Saved RDS: ", out_rds),
    paste0("Saved Table S9: ", out_doc)
)

writeLines(txt, out_txt)

# ------------------------------------------------------------
# 7. Word Table S9
# ------------------------------------------------------------

pword <- data.table(
    `Observed M1 interaction HR (robust 95% CI)` =
        fmt_hr(obs_hr, obs_lo, obs_hi),
    `Successful bootstrap replicates` =
        paste0(n_ok, "/", B, " (", sprintf("%.1f%%", 100 * rate_ok), ")"),
    `Bootstrap median HR` =
        sprintf("%.2f", qi[2]),
    `Bootstrap percentile 95% CI` =
        sprintf("%.2f–%.2f", qi[1], qi[3]),
    `HR>1` =
        sprintf("%.1f%%", 100 * p_gt1),
    `HR>1.5` =
        sprintf("%.1f%%", 100 * p_gt15),
    `HR>2` =
        sprintf("%.1f%%", 100 * p_gt2)
)

cword <- contrasts[, .(
    Contrast,
    `Observed estimate` = sprintf("%.2f", Observed),
    `Bootstrap median` = sprintf("%.2f", Bootstrap_median),
    `Bootstrap percentile 95% CI` =
        sprintf("%.2f–%.2f", Bootstrap_LCL95, Bootstrap_UCL95)
)]

fword <- copy(failures)
setnames(fword, c("STATUS", "N"), c("Bootstrap status", "Replicates"))

three_line <- function(x, size = 7.4) {
    ft <- flextable(x)
    ft <- border_remove(ft)
    top <- fp_border(color = "black", width = 1)
    mid <- fp_border(color = "black", width = .6)
    ft <- hline_top(ft, border = top, part = "header")
    ft <- hline_bottom(ft, border = mid, part = "header")
    ft <- hline_bottom(ft, border = top, part = "body")
    ft <- font(ft, fontname = "Arial", part = "all")
    ft <- fontsize(ft, size = size, part = "all")
    ft <- bold(ft, part = "header")
    ft <- align(ft, align = "center", part = "all")
    ft <- align(ft, j = 1, align = "left", part = "all")
    ft <- autofit(ft)
    ft
}

doc <- read_docx()

doc <- body_set_default_section(
    doc,
    prop_section(
        page_size = page_size(orient = "landscape"),
        page_margins = page_mar(
            top = .45,
            bottom = .45,
            left = .45,
            right = .45
        )
    )
)

doc <- body_add_par(
    doc,
    "Supplementary Table S9. Full estimation-pipeline bootstrap of the exploratory overlap-weighted treatment-by-bypass interaction",
    style = "heading 1"
)

doc <- body_add_par(doc, "A. Primary interaction bootstrap")
doc <- body_add_flextable(doc, three_line(pword, 7.0))

doc <- body_add_par(doc, "")
doc <- body_add_par(doc, "B. Contrasts derived from the same M1 bootstrap")
doc <- body_add_flextable(doc, three_line(cword, 7.4))

doc <- body_add_par(doc, "")
doc <- body_add_par(doc, "C. Bootstrap estimation status")
doc <- body_add_flextable(doc, three_line(fword, 7.4))

doc <- body_add_par(
    doc,
    paste0(
        "Notes: Two thousand patient-level bootstrap samples were drawn with replacement from the frozen 319-patient ",
        "propensity-score-complete B1-14B analysis population. In each replicate, the fixed propensity model was refit, ",
        "overlap weights were recalculated, and the unchanged weighted Cox treatment-by-bypass model was refit. ",
        "Percentile 95% confidence intervals summarize the empirical bootstrap distribution. This analysis does not ",
        "address unmeasured confounding and does not establish causal or predictive-biomarker inference. Classification: ",
        gate, "."
    )
)

print(doc, target = out_doc)

# ------------------------------------------------------------
# 8. Console summary
# ------------------------------------------------------------

cat("\nBOOTSTRAP STATUS\n")
print(failures)

cat(
    "\nSuccessful: ", n_ok, "/", B,
    " (", sprintf("%.1f%%", 100 * rate_ok), ")\n",
    sep = ""
)

cat(
    "Observed M1 interaction HR: ",
    sprintf("%.3f", obs_hr),
    " (robust 95% CI ",
    sprintf("%.3f", obs_lo),
    "–",
    sprintf("%.3f", obs_hi),
    ")\n",
    sep = ""
)

cat(
    "Bootstrap median HR: ",
    sprintf("%.3f", qi[2]),
    "\nBootstrap percentile 95% CI: ",
    sprintf("%.3f", qi[1]),
    "–",
    sprintf("%.3f", qi[3]),
    "\n",
    sep = ""
)

cat(
    "Replicates HR>1: ",
    sprintf("%.1f%%", 100 * p_gt1),
    "\n",
    sep = ""
)

cat("\nBOOTSTRAP GATE: ", gate, "\n", sep = "")
cat(wording, "\n")

cat("\nB1-14C COMPLETE\n")
cat("Audit: ", out_txt, "\n", sep = "")
cat("Table S9: ", out_doc, "\n", sep = "")
cat("RDS: ", out_rds, "\n", sep = "")
