# ============================================================
# A7 / B-line
# 16A_B1_Figure4_exploratory_treatment_heterogeneity.R
#
# REVISED MAIN MANUSCRIPT FIGURE 4
# Exploratory treatment-by-bypass heterogeneity in rwPFS
#
# This revision fixes only figure layout/readability.
# Statistical definitions and frozen estimates are unchanged.
#
# PANELS
#   A. Overlap-weighted descriptive rwPFS curves
#   B. Anti-EGFR vs bevacizumab HR within bypass strata
#   C. Treatment-by-bypass interaction robustness
#
# INPUTS
#   04_results/B1_14_bevacizumab_comparator_feasibility.rds
#   04_results/B1_14B_exploratory_comparator_rwPFS.rds
#   04_results/B1_14C_full_pipeline_bootstrap.rds
#
# OUTPUTS -- overwrite existing files
#   05_figures/Figure4_exploratory_treatment_heterogeneity.png
#   05_figures/Figure4_exploratory_treatment_heterogeneity.tiff
#   06_logs_and_audit/Figure4_exploratory_treatment_heterogeneity_audit.txt
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

results_dir <- file.path(project_root, "04_results")
figures_dir <- file.path(project_root, "05_figures")
audit_dir   <- file.path(project_root, "06_logs_and_audit")

for (d0 in c(results_dir, figures_dir, audit_dir)) {
    if (!dir.exists(d0)) {
        dir.create(d0, recursive = TRUE, showWarnings = FALSE)
    }
}

b14_file <- file.path(results_dir, "B1_14_bevacizumab_comparator_feasibility.rds")
b14b_file <- file.path(results_dir, "B1_14B_exploratory_comparator_rwPFS.rds")
b14c_file <- file.path(results_dir, "B1_14C_full_pipeline_bootstrap.rds")

png_file <- file.path(figures_dir, "Figure4_exploratory_treatment_heterogeneity.png")
tiff_file <- file.path(figures_dir, "Figure4_exploratory_treatment_heterogeneity.tiff")
audit_file <- file.path(audit_dir, "Figure4_exploratory_treatment_heterogeneity_audit.txt")

for (f0 in c(b14_file, b14b_file, b14c_file)) {
    if (!file.exists(f0)) stop("Missing required input:\n", f0)
}

pkgs <- c("data.table", "survival", "ggplot2", "patchwork", "scales")
for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(data.table)
library(survival)
library(ggplot2)
library(patchwork)
library(scales)

cat("\n============================================================\n")
cat("A7 B1-16A: REBUILD MAIN FIGURE 4\n")
cat("============================================================\n\n")

# 1. Load frozen objects
b14  <- readRDS(b14_file)
b14b <- readRDS(b14b_file)
b14c <- readRDS(b14c_file)

if (!identical(as.character(b14$gate_decision), "GO")) {
    stop("B1-14 feasibility gate is not GO.")
}

if (!all(c("analysis_data", "model_contrasts") %in% names(b14b))) {
    stop("B1-14B RDS is missing analysis_data/model_contrasts.")
}
if (!"ps_formula" %in% names(b14)) {
    stop("B1-14 RDS is missing the frozen PS formula.")
}
if (!all(c("primary", "contrasts") %in% names(b14c))) {
    stop("B1-14C RDS is missing primary/contrasts.")
}

d  <- as.data.table(copy(b14b$analysis_data))
mc <- as.data.table(copy(b14b$model_contrasts))

m0 <- mc[MODEL_ID == "M0"]
m1 <- mc[MODEL_ID == "M1"]
m2 <- mc[MODEL_ID == "M2"]

if (nrow(m0) != 1L || nrow(m1) != 1L || nrow(m2) != 1L) {
    stop("Could not uniquely identify M0/M1/M2.")
}

boot_primary <- as.data.table(copy(b14c$primary))
if (nrow(boot_primary) != 1L) {
    stop("B1-14C primary bootstrap summary is not a single row.")
}

# 2. Hard anchors
anchors <- d[
    ,
    .(N = .N, EVENTS = sum(RWPFS_EVENT_B14B_R)),
    by = .(FIRST_BIOLOGIC_R, BYPASS_HIGH_CONFIDENCE_R)
][order(FIRST_BIOLOGIC_R, BYPASS_HIGH_CONFIDENCE_R)]

expected <- data.table(
    FIRST_BIOLOGIC_R = c("anti-EGFR", "anti-EGFR", "bevacizumab", "bevacizumab"),
    BYPASS_HIGH_CONFIDENCE_R = c(0L, 1L, 0L, 1L),
    N = c(165L, 23L, 117L, 14L),
    EVENTS = c(160L, 21L, 111L, 13L)
)

anchors_check <- merge(
    expected,
    anchors,
    by = c("FIRST_BIOLOGIC_R", "BYPASS_HIGH_CONFIDENCE_R"),
    suffixes = c("_EXPECTED", "_OBSERVED"),
    all = TRUE
)

if (
    nrow(d) != 319L ||
    sum(d$RWPFS_EVENT_B14B_R) != 305L ||
    any(anchors_check$N_EXPECTED != anchors_check$N_OBSERVED) ||
    any(anchors_check$EVENTS_EXPECTED != anchors_check$EVENTS_OBSERVED)
) {
    stop("Four-cell B1-14B anchors failed. Figure not generated.")
}

if (abs(m1$INTERACTION_HR - 3.083038) > 0.001) {
    stop("M1 interaction HR does not match frozen B1-14B.")
}

if (
    boot_primary$Successful != 2000L ||
    abs(boot_primary$Bootstrap_median_HR - 3.12982) > 0.01 ||
    abs(boot_primary$Bootstrap_LCL95 - 1.407512) > 0.01 ||
    abs(boot_primary$Bootstrap_UCL95 - 7.294777) > 0.02
) {
    stop("B1-14C bootstrap anchors failed.")
}

# 3. Reproduce frozen PS + overlap weights
ps_fit <- glm(b14$ps_formula, data = d, family = binomial())

if (!isTRUE(ps_fit$converged) || any(!is.finite(coef(ps_fit)))) {
    stop("Frozen propensity model could not be reproduced.")
}

d[, PS_FIG4_R := as.numeric(predict(ps_fit, type = "response"))]
d[, OW_FIG4_R := fifelse(ANTI_EGFR_TREAT_R == 1L, 1 - PS_FIG4_R, PS_FIG4_R)]

m1_check <- coxph(
    Surv(RWPFS_DAYS_B14B_R, RWPFS_EVENT_B14B_R) ~
        ANTI_EGFR_TREAT_R * BYPASS_HIGH_CONFIDENCE_R,
    data = d,
    weights = OW_FIG4_R,
    ties = "efron"
)

interaction_term <- "ANTI_EGFR_TREAT_R:BYPASS_HIGH_CONFIDENCE_R"
m1_reproduced_hr <- exp(coef(m1_check)[interaction_term])

if (
    !is.finite(m1_reproduced_hr) ||
    abs(m1_reproduced_hr - m1$INTERACTION_HR) > 1e-5
) {
    stop(
        sprintf(
            "M1 reproduction failed: saved %.8f vs reproduced %.8f",
            m1$INTERACTION_HR,
            m1_reproduced_hr
        )
    )
}

# 4. Four-group labels
d[, GROUP_FIG4_R := factor(
    paste(FIRST_BIOLOGIC_R, BYPASS_HIGH_CONFIDENCE_R, sep = "__"),
    levels = c(
        "bevacizumab__0",
        "bevacizumab__1",
        "anti-EGFR__0",
        "anti-EGFR__1"
    ),
    labels = c(
        "Bevacizumab / bypass−",
        "Bevacizumab / bypass+",
        "Anti-EGFR / bypass−",
        "Anti-EGFR / bypass+"
    )
)]

# 5. Panel A
sf_w <- survfit(
    Surv(RWPFS_DAYS_B14B_R, RWPFS_EVENT_B14B_R) ~ GROUP_FIG4_R,
    data = d,
    weights = OW_FIG4_R,
    conf.type = "log"
)

s <- summary(sf_w)

curve_dt <- data.table(
    time = s$time,
    surv = s$surv,
    strata = as.character(s$strata)
)

curve_dt[, group := sub("^GROUP_FIG4_R=", "", strata)]
curve_dt[, group := factor(group, levels = levels(d$GROUP_FIG4_R))]

curve_dt <- rbind(
    data.table(
        time = 0,
        surv = 1,
        strata = paste0("GROUP_FIG4_R=", levels(d$GROUP_FIG4_R)),
        group = factor(levels(d$GROUP_FIG4_R), levels = levels(d$GROUP_FIG4_R))
    ),
    curve_dt,
    fill = TRUE
)

curve_dt[, treatment := fifelse(
    grepl("^Bevacizumab", as.character(group)),
    "Bevacizumab",
    "Anti-EGFR"
)]

curve_dt[, bypass := fifelse(
    grepl("bypass\\+$", as.character(group)),
    "Bypass+",
    "Bypass−"
)]

curve_dt[, treatment := factor(treatment, levels = c("Bevacizumab", "Anti-EGFR"))]
curve_dt[, bypass := factor(bypass, levels = c("Bypass−", "Bypass+"))]
curve_dt <- curve_dt[order(group, time)]

# [19] p-value style aligned with the manuscript: lowercase p, leading zero,
#      spaces around the operator (e.g. "p = 0.003", "p < 0.001).
interaction_subtitle <- sprintf(
    "Primary overlap-weighted interaction HR %.2f (95%% CI %.2f–%.2f); %s",
    m1$INTERACTION_HR,
    m1$INTERACTION_LCL95,
    m1$INTERACTION_UCL95,
    ifelse(
        m1$INTERACTION_P < 0.001,
        "p < 0.001",
        sprintf("p = %.3f", m1$INTERACTION_P)
    )
)

pA <- ggplot(
    curve_dt,
    aes(
        x = time,
        y = surv,
        group = group,
        color = treatment,
        linetype = bypass
    )
) +
    geom_step(linewidth = 0.9, direction = "hv") +
    scale_color_manual(
        values = c(
            "Bevacizumab" = "#0072B2",
            "Anti-EGFR" = "#D55E00"
        ),
        name = "Treatment"
    ) +
    scale_linetype_manual(
        values = c(
            "Bypass−" = "solid",
            "Bypass+" = "22"
        ),
        name = "Bypass status"
    ) +
    scale_x_continuous(
        limits = c(0, 730),
        breaks = c(0, 180, 365, 540, 730),
        expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = percent_format(accuracy = 1),
        expand = expansion(mult = c(0, 0.025))
    ) +
    labs(
        tag = "A",
        title = "Overlap-weighted rwPFS by treatment and bypass status",
        subtitle = interaction_subtitle,
        x = "Days from first relevant biologic exposure",
        y = "rwPFS probability"
    ) +
    guides(
        color = guide_legend(
            order = 1,
            override.aes = list(linewidth = 1.1)
        ),
        linetype = guide_legend(
            order = 2,
            override.aes = list(linewidth = 1.1)
        )
    ) +
    theme_classic(base_size = 10.5) +
    theme(
        plot.tag = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 11.2, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 9.0, margin = margin(b = 7)),
        axis.title = element_text(size = 9.6),
        axis.text = element_text(size = 9.0),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold", size = 8.8),
        legend.text = element_text(size = 8.6),
        legend.key.width = unit(19, "pt"),
        plot.margin = margin(5, 8, 5, 8)
    )

# 6. Panel B
forest_b <- data.table(
    label = c("Bypass−", "Bypass+"),
    hr = c(
        m1$ANTI_VS_BEV_BYPASS_NEG_HR,
        m1$ANTI_VS_BEV_BYPASS_POS_HR
    ),
    lo = c(
        m1$ANTI_VS_BEV_BYPASS_NEG_LCL95,
        m1$ANTI_VS_BEV_BYPASS_POS_LCL95
    ),
    hi = c(
        m1$ANTI_VS_BEV_BYPASS_NEG_UCL95,
        m1$ANTI_VS_BEV_BYPASS_POS_UCL95
    )
)

forest_b[, label := factor(label, levels = rev(c("Bypass−", "Bypass+")))]
forest_b[, display := sprintf("%.2f (%.2f–%.2f)", hr, lo, hi)]

pB <- ggplot(forest_b, aes(x = hr, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45) +
    geom_errorbarh(
        aes(xmin = lo, xmax = hi),
        height = 0.12,
        linewidth = 0.7
    ) +
    geom_point(size = 2.8) +
    geom_text(
        aes(label = display),
        x = 10.8,
        hjust = 1,
        size = 3.1
    ) +
    scale_x_log10(
        limits = c(0.6, 12),
        breaks = c(0.5, 1, 2, 4, 8),
        labels = c("0.5", "1", "2", "4", "8")
    ) +
    labs(
        tag = "B",
        title = "Treatment contrast by bypass status",
        subtitle = "HR >1 indicates a higher rwPFS event rate with anti-EGFR",
        x = "Anti-EGFR vs bevacizumab HR (log scale)",
        y = NULL
    ) +
    theme_classic(base_size = 10.0) +
    theme(
        plot.tag = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 10.2, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 8.2, margin = margin(b = 6)),
        axis.title.x = element_text(size = 8.8),
        axis.text = element_text(size = 8.6),
        plot.margin = margin(5, 8, 5, 8)
    )

# 7. Panel C
forest_c <- data.table(
    label = c(
        "Crude",
        "Overlap-weighted M1",
        "Doubly adjusted M2",
        "Full-pipeline bootstrap"
    ),
    hr = c(
        m0$INTERACTION_HR,
        m1$INTERACTION_HR,
        m2$INTERACTION_HR,
        boot_primary$Bootstrap_median_HR
    ),
    lo = c(
        m0$INTERACTION_LCL95,
        m1$INTERACTION_LCL95,
        m2$INTERACTION_LCL95,
        boot_primary$Bootstrap_LCL95
    ),
    hi = c(
        m0$INTERACTION_UCL95,
        m1$INTERACTION_UCL95,
        m2$INTERACTION_UCL95,
        boot_primary$Bootstrap_UCL95
    ),
    type = c("Model", "Primary", "Sensitivity", "Bootstrap")
)

forest_c[, label := factor(
    label,
    levels = rev(c(
        "Crude",
        "Overlap-weighted M1",
        "Doubly adjusted M2",
        "Full-pipeline bootstrap"
    ))
)]

forest_c[, display := sprintf("%.2f (%.2f–%.2f)", hr, lo, hi)]

pC <- ggplot(forest_c, aes(x = hr, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45) +
    geom_errorbarh(
        aes(xmin = lo, xmax = hi),
        height = 0.12,
        linewidth = 0.7
    ) +
    geom_point(aes(shape = type), size = 2.8) +
    geom_text(
        aes(label = display),
        x = 10.8,
        hjust = 1,
        size = 3.0
    ) +
    scale_x_log10(
        limits = c(0.9, 12),
        breaks = c(1, 2, 4, 8),
        labels = c("1", "2", "4", "8")
    ) +
    scale_shape_manual(
        values = c(
            "Model" = 16,
            "Primary" = 17,
            "Sensitivity" = 15,
            "Bootstrap" = 18
        ),
        guide = "none"
    ) +
    labs(
        tag = "C",
        title = "Interaction robustness",
        subtitle = "Bootstrap: median estimate with percentile 95% CI",
        x = "Treatment-by-bypass interaction HR (log scale)",
        y = NULL
    ) +
    theme_classic(base_size = 10.0) +
    theme(
        plot.tag = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 10.2, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 8.2, margin = margin(b = 6)),
        axis.title.x = element_text(size = 8.8),
        axis.text = element_text(size = 8.3),
        plot.margin = margin(5, 8, 5, 8)
    )

# 8. Assemble
bottom_row <- pB + pC + plot_layout(widths = c(1, 1.12))

fig4 <- pA / bottom_row +
    plot_layout(heights = c(1.28, 1)) &
    theme(
        plot.background = element_rect(fill = "white", color = NA)
    )

# 9. Save -- overwrite
ggsave(
    filename = png_file,
    plot = fig4,
    width = 8.5,
    height = 7.8,
    units = "in",
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = tiff_file,
    plot = fig4,
    width = 8.5,
    height = 7.8,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

# 10. Audit
audit_lines <- c(
    "A7 B1-16A MAIN FIGURE 4 AUDIT",
    "==============================",
    "",
    "FIGURE",
    "------",
    "Figure 4. Exploratory treatment-by-bypass heterogeneity in rwPFS.",
    "Panel A: overlap-weighted descriptive rwPFS curves.",
    "Panel B: anti-EGFR vs bevacizumab HRs within bypass strata.",
    "Panel C: interaction robustness across M0, M1, M2, and full-pipeline bootstrap.",
    "",
    "HARD COHORT ANCHORS",
    "-------------------",
    paste0("PS-complete N = ", nrow(d)),
    paste0("rwPFS events = ", sum(d$RWPFS_EVENT_B14B_R)),
    capture.output(print(anchors)),
    "",
    "M1 REPRODUCTION",
    "---------------",
    paste0("Saved M1 interaction HR = ", sprintf("%.8f", m1$INTERACTION_HR)),
    paste0("Reproduced M1 interaction HR = ", sprintf("%.8f", m1_reproduced_hr)),
    "Reproduction: PASS",
    "",
    "PANEL B VALUES",
    "--------------",
    capture.output(print(forest_b[, .(label, hr, lo, hi)])),
    "",
    "PANEL C VALUES",
    "--------------",
    capture.output(print(forest_c[, .(label, hr, lo, hi)])),
    "",
    "INTERPRETATION GUARDRAILS",
    "-------------------------",
    "- M1 is the primary exploratory comparator interaction.",
    "- Panel A is overlap-weighted descriptive visualization, not a causal standardized treatment effect.",
    "- Interaction HR >1 means the anti-EGFR-versus-bevacizumab HR is relatively less favorable in bypass-positive than bypass-negative patients.",
    "- The comparator analysis remains observational and exploratory.",
    "- No treatment-selection or established predictive-biomarker claim is made.",
    "",
    paste0("PNG: ", png_file),
    paste0("TIFF: ", tiff_file)
)

writeLines(audit_lines, audit_file)

cat("\n============================================================\n")
cat("FIGURE 4 COMPLETE\n")
cat("============================================================\n")

cat("\nHard anchors:\n")
print(anchors)

cat(
    "\nM1 interaction HR: ",
    sprintf("%.3f", m1$INTERACTION_HR),
    " (95% CI ",
    sprintf("%.3f", m1$INTERACTION_LCL95),
    "–",
    sprintf("%.3f", m1$INTERACTION_UCL95),
    "), P=",
    sprintf("%.4f", m1$INTERACTION_P),
    "\n",
    sep = ""
)

cat(
    "Full-pipeline bootstrap: median HR ",
    sprintf("%.3f", boot_primary$Bootstrap_median_HR),
    " (percentile 95% CI ",
    sprintf("%.3f", boot_primary$Bootstrap_LCL95),
    "–",
    sprintf("%.3f", boot_primary$Bootstrap_UCL95),
    ")\n",
    sep = ""
)

cat("\nSaved PNG:\n", png_file, "\n")
cat("\nSaved TIFF:\n", tiff_file, "\n")
cat("\nSaved audit:\n", audit_file, "\n")
