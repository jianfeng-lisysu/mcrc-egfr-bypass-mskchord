# ============================================================
# A7 / B-line
# 11C_B1_11_refine_main_figures_for_submission.R
# FINAL FIGURE REFINEMENT ONLY -- no new primary analysis
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory
intermediate_dir <- file.path(project_root, "03_intermediate")
results_dir <- file.path(project_root, "04_results")
fig_dir <- file.path(project_root, "05_figures")
audit_dir <- file.path(project_root, "06_logs_and_audit")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(audit_dir)) dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

locked_file <- file.path(intermediate_dir, "A7_B1_locked_analysis.rds")
b107_file <- file.path(results_dir, "B1_07_primary_models.rds")
b108_file <- file.path(results_dir, "B1_08_robustness_validation.rds")
b110_file <- file.path(results_dir, "B1_10_refined_TTNTD_PH_validation.rds")

for (f in c(locked_file, b107_file, b108_file, b110_file)) {
  if (!file.exists(f)) stop("Missing required input:\n", f)
}

packages <- c("data.table", "survival", "ggplot2", "patchwork", "scales", "grid")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(data.table)
library(survival)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)

cat("\n============================================================\n")
cat("B1-11C FINAL FIGURE REFINEMENT\n")
cat("============================================================\n\n")

dt <- as.data.table(readRDS(locked_file))
b107 <- readRDS(b107_file)
b108 <- readRDS(b108_file)
b110 <- readRDS(b110_file)
ttdt <- as.data.table(b110$patient_data)

if (nrow(dt) != 191L || dt[BYPASS_HIGH_CONFIDENCE_R == 1L, .N] != 23L) {
  stop("Locked cohort anchor failed.")
}

COL_NEG        <- "#0072B2"
COL_POS        <- "#D55E00"
COL_TEAL       <- "#009E73"
COL_GOLD       <- "#E69F00"
COL_DARK       <- "#222222"
COL_LIGHT      <- "#EDF1F4"
COL_LIGHT_BLUE <- "#E7F2F8"

fmt_num_vec <- function(x, digits = 2) {
  vapply(x, function(z) {
    if (length(z) == 0L || is.na(z) || !is.finite(z)) return("NA")
    formatC(z, format = "f", digits = digits)
  }, character(1))
}

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = COL_DARK),
      axis.title = element_text(face = "bold", color = COL_DARK),
      plot.title = element_text(face = "bold", color = COL_DARK, size = base_size + 1),
      legend.title = element_blank(),
      legend.position = "top",
      plot.margin = margin(6, 8, 6, 8)
    )
}

save_pub <- function(plot, basename, width, height) {
  ggsave(file.path(fig_dir, paste0(basename, ".png")), plot = plot,
         width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(basename, ".tiff")), plot = plot,
         width = width, height = height, units = "in", dpi = 600,
         compression = "lzw", bg = "white")
}

survfit_to_df <- function(sf) {
  s <- summary(sf)
  out <- data.table(time = s$time, surv = s$surv, strata = as.character(s$strata))
  out[, group := fifelse(grepl("=1$", strata),
                         "High-confidence bypass+",
                         "High-confidence bypass−")]
  out[, group := factor(group,
                        levels = c("High-confidence bypass−", "High-confidence bypass+"))]
  out
}

make_km_with_risk <- function(data, time_col, event_col, title, x_breaks, x_limit) {
  sf <- survfit(
    Surv(get(time_col), get(event_col)) ~ BYPASS_HIGH_CONFIDENCE_R,
    data = data,
    conf.type = "log-log"
  )
  kd <- survfit_to_df(sf)

  p_km <- ggplot(kd, aes(x = time, y = surv, color = group, linetype = group)) +
    geom_step(linewidth = 0.95) +
    scale_color_manual(
      values = c("High-confidence bypass−" = COL_NEG,
                 "High-confidence bypass+" = COL_POS),
      breaks = c("High-confidence bypass−", "High-confidence bypass+")
    ) +
    scale_linetype_manual(
      values = c("High-confidence bypass−" = "solid",
                 "High-confidence bypass+" = "dashed"),
      breaks = c("High-confidence bypass−", "High-confidence bypass+")
    ) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1),
                       breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    scale_x_continuous(limits = c(0, x_limit), breaks = x_breaks,
                       expand = expansion(mult = c(0, 0.01))) +
    labs(title = title,
         x = "Days from first anti-EGFR exposure",
         y = "Survival probability") +
    theme_pub(10) +
    theme(legend.position = "top")

  risk_dt <- rbindlist(lapply(c(0L, 1L), function(g) {
    d0 <- data[BYPASS_HIGH_CONFIDENCE_R == g]
    data.table(
      time = x_breaks,
      n_risk = vapply(
        x_breaks,
        function(t0) sum(d0[[time_col]] >= t0),
        integer(1)
      ),
      group = if (g == 0L) "Bypass−" else "Bypass+"
    )
  }))

  risk_dt[, group := factor(
    group,
    # ggplot2 draws the last discrete y level on top:
    # this yields Bypass− on top and Bypass+ below.
    levels = c("Bypass+", "Bypass−")
  )]

  # Prevent time-zero counts from being clipped at the left plotting boundary.
  risk_dt[, plot_time := fcase(
    time == min(x_breaks),
    time + 0.025 * x_limit,

    time == max(x_breaks),
    time - 0.025 * x_limit,

    default = time
  )]

  p_risk <- ggplot(
    risk_dt,
    aes(
      x = plot_time,
      y = group,
      label = n_risk,
      color = group
    )
  ) +
    geom_text(size = 3.0, show.legend = FALSE) +
    scale_color_manual(values = c(
      "Bypass−" = COL_NEG,
      "Bypass+" = COL_POS
    )) +
    scale_x_continuous(
      limits = c(0, x_limit),
      breaks = x_breaks,
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(x = NULL, y = "No. at risk") +
    theme_void(base_size = 9, base_family = "Arial") +
    theme(
      axis.text.y = element_text(size = 8, color = COL_DARK),
      axis.title.y = element_text(face = "bold", angle = 0, vjust = 0.5,
                                  margin = margin(r = 8)),
      plot.margin = margin(0, 8, 2, 8)
    )

  p_km / p_risk + plot_layout(heights = c(4.7, 0.8))
}

# ============================================================
# FIGURE 1 -- A BLACK/WHITE ONLY
# ============================================================

flow_nodes <- data.table(
  x = 1,
  y = c(5, 4, 3, 2, 1),
  label = c(
    "MSK-CHORD colorectal cancer\nn = 5,543",
    "Clinical-genomic eligibility screen\nleft-sided stage IV, MSS,\nno detected KRAS/NRAS mutation,\nBRAF V600E-negative\nn = 803",
    "First relevant biologic: anti-EGFR\nn = 230",
    "Sequencing completed on/before\nfirst anti-EGFR exposure\nn = 191",
    "Locked analysis cohort\nhigh-confidence bypass+ n = 23\nbypass− n = 168"
  )
)

arrow_dt <- data.table(
  x = 1, xend = 1,
  y = c(4.67, 3.67, 2.67, 1.67),
  yend = c(4.33, 3.33, 2.33, 1.33)
)

p1a <- ggplot() +
  geom_segment(
    data = arrow_dt,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.6, color = "black",
    arrow = arrow(length = unit(0.08, "inches"))
  ) +
  geom_label(
    data = flow_nodes,
    aes(x = x, y = y, label = label),
    size = 3.1, color = "black", fill = "white",
    label.size = 0.45, label.padding = unit(0.20, "lines")
  ) +
  coord_cartesian(xlim = c(0.25, 1.75), ylim = c(0.55, 5.45)) +
  labs(title = "A  Cohort derivation") +
  theme_void(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 12))

# Figure 1B can retain restrained colors.
timeline <- data.table(
  x = c(0, 1.5, 3, 4.5, 6),
  label = c(
    "Stage IV\nregistry diagnosis",
    "Pre-treatment\ntumor specimen",
    "MSK-IMPACT\nsequencing",
    "First anti-EGFR\n(index)",
    "rwPFS / OS /\nTTNTD follow-up"
  ),
  point_color = c(COL_DARK, COL_NEG, COL_NEG, COL_POS, COL_TEAL)
)

p1b <- ggplot() +
  geom_segment(aes(x = 0, xend = 6, y = 1, yend = 1),
               linewidth = 0.75, color = COL_DARK) +
  geom_point(data = timeline, aes(x = x, y = 1, color = point_color), size = 3) +
  scale_color_identity() +
  geom_text(data = timeline, aes(x = x, y = 0.83, label = label),
            size = 3.0, vjust = 1) +
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 1.11, ymax = 1.31,
           fill = COL_LIGHT_BLUE, color = COL_NEG, linewidth = 0.35) +
  annotate("text", x = 2.25, y = 1.21,
           label = "Genomic ascertainment must precede anti-EGFR index", size = 3.0) +
  annotate("text", x = 3.0, y = 1.50,
           label = "Frozen exposure: high-confidence EGFR-bypass alterations",
           size = 3.2, fontface = "bold") +
  coord_cartesian(xlim = c(-0.35, 6.35), ylim = c(0.15, 1.70)) +
  labs(title = "B  Study time zero and endpoint framework") +
  theme_void(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 12))

fig1 <- p1a / p1b + plot_layout(heights = c(1.45, 0.68))
save_pub(fig1, "Figure1_study_design_and_cohort_flow", 8.2, 8.4)

# ============================================================
# FIGURE 2
# ============================================================

component_vars <- c(
  "ERBB2_AMP_R", "ERBB2_ACTIVATING_R", "MET_AMP_R", "PIK3CA_EX20_R",
  "PTEN_LOF_R", "KRAS_AMP_R", "MAP2K1_HOTSPOT_R", "NF1_LOF_R"
)

component_labels <- c(
  "ERBB2 amplification", "ERBB2 activating mutation", "MET amplification",
  "PIK3CA exon 20", "PTEN loss-of-function", "KRAS amplification",
  "MAP2K1 hotspot", "NF1 loss-of-function"
)
names(component_labels) <- component_vars

comp_counts <- rbindlist(lapply(component_vars, function(v) {
  data.table(component = v,
             label = component_labels[[v]],
             n = sum(dt[[v]] == 1L, na.rm = TRUE))
}))
comp_counts[, prevalence := n / nrow(dt)]
comp_counts[, label := factor(label, levels = rev(component_labels))]

p2a <- ggplot(comp_counts, aes(x = prevalence, y = label)) +
  geom_col(fill = COL_NEG, width = 0.62) +
  geom_text(aes(label = paste0(n, " (", percent(prevalence, accuracy = 0.1), ")")),
            hjust = -0.08, size = 3.0) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.22))) +
  labs(title = "A  Prevalence of high-confidence bypass components",
       x = "Prevalence in the locked cohort", y = NULL) +
  theme_pub(9)

carriers <- copy(dt[BYPASS_HIGH_CONFIDENCE_R == 1L])
carriers[, burden := rowSums(.SD, na.rm = TRUE), .SDcols = component_vars]
setorderv(carriers,
          c("burden", "ERBB2_AMP_R", "ERBB2_ACTIVATING_R", "PIK3CA_EX20_R", "KRAS_AMP_R"),
          c(-1, -1, -1, -1, -1))
carriers[, patient_order := seq_len(.N)]

heat <- rbindlist(lapply(component_vars, function(v) {
  data.table(
    patient_order = carriers$patient_order,
    patient = sprintf("P%02d", carriers$patient_order),
    component = component_labels[[v]],
    present = carriers[[v]]
  )
}))
heat[, patient := factor(patient, levels = sprintf("P%02d", seq_len(nrow(carriers))))]
heat[, component := factor(component, levels = rev(component_labels))]

p2b <- ggplot(heat, aes(x = patient, y = component, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.30) +
  scale_fill_manual(values = c("0" = COL_LIGHT, "1" = COL_POS),
                    labels = c("Absent", "Present")) +
  labs(title = "B  Component landscape among 23 high-confidence bypass+ patients",
       x = "High-confidence bypass+ patients", y = NULL, fill = NULL) +
  theme_pub(8.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.8),
        legend.position = "top")

fig2 <- p2a / p2b + plot_layout(heights = c(0.85, 1.25))
save_pub(fig2, "Figure2_high_confidence_bypass_landscape", 8.6, 8.0)

# ============================================================
# FIGURE 3
# ============================================================

p3a <- make_km_with_risk(
  dt, "RWPFS_DAYS_LOCKED_R", "RWPFS_EVENT_LOCKED_R",
  "A  Real-world progression-free survival",
  c(0, 200, 400, 600, 800), 800
)

p3b <- make_km_with_risk(
  dt, "OS_DAYS_LOCKED_R", "OS_EVENT_LOCKED_R",
  "B  Overall survival",
  c(0, 500, 1000, 1400), 1400
)

forest_main <- as.data.table(b107$primary_results)
forest_main <- forest_main[ANALYSIS %chin% c("Crude", "Primary adjusted")]
forest_main[, LABEL := paste0(OUTCOME, " — ", ANALYSIS)]
forest_main[, TEXT := paste0(fmt_num_vec(HR, 2), " (",
                            fmt_num_vec(LCL95, 2), "–",
                            fmt_num_vec(UCL95, 2), ")")]
forest_main[, LABEL := factor(
  LABEL,
  levels = rev(c("rwPFS — Crude", "rwPFS — Primary adjusted",
                 "OS — Crude", "OS — Primary adjusted"))
)]

p3c <- ggplot(forest_main, aes(x = HR, y = LABEL)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = COL_DARK) +
  geom_errorbarh(aes(xmin = LCL95, xmax = UCL95),
                 height = 0.13, linewidth = 0.7, color = COL_POS) +
  geom_point(size = 2.5, color = COL_POS) +
  geom_text(aes(x = UCL95 * 1.16, label = TEXT),
            hjust = 0, size = 3.0, color = COL_DARK) +
  scale_x_log10(breaks = c(1, 2, 4),
                limits = c(0.75, max(forest_main$UCL95) * 1.75)) +
  labs(title = "C  Hazard-ratio estimates",
       x = "Hazard ratio (log scale)", y = NULL) +
  theme_pub(9) +
  theme(legend.position = "none")

rmst <- as.data.table(b108$rmst)
rmst_long <- rbindlist(list(
  rmst[, .(TAU_DAYS, group = "High-confidence bypass−", value = RMST_BYPASS_NEG_DAYS)],
  rmst[, .(TAU_DAYS, group = "High-confidence bypass+", value = RMST_BYPASS_POS_DAYS)]
))
rmst_long[, group := factor(group,
                            levels = c("High-confidence bypass−", "High-confidence bypass+"))]
rmst_long[, horizon := factor(paste0(TAU_DAYS, "-day horizon"),
                              levels = c("180-day horizon", "365-day horizon"))]
rmst_long[, label := sprintf("%.1f", value)]

p3d <- ggplot(rmst_long, aes(x = horizon, y = value, fill = group)) +
  geom_col(width = 0.68, position = position_dodge(width = 0.74)) +
  geom_text(aes(label = label), position = position_dodge(width = 0.74),
            vjust = -0.30, size = 3.0) +
  scale_fill_manual(
    values = c("High-confidence bypass−" = COL_NEG,
               "High-confidence bypass+" = COL_POS),
    breaks = c("High-confidence bypass−", "High-confidence bypass+")
  ) +
  labs(title = "D  Restricted mean rwPFS",
       x = NULL, y = "RMST, days", fill = NULL) +
  theme_pub(9) +
  theme(legend.position = "top")

fig3 <- (p3a | p3b) / (p3c | p3d) + plot_layout(heights = c(1.08, 0.92))
save_pub(fig3, "Figure3_primary_rwPFS_OS", 11.4, 9.0)

# ============================================================
# FIGURE 4 MAIN -- 2 PANELS ONLY
# ============================================================

primary_robust <- as.data.table(b108$primary_result)
primary_robust[, label := "Primary adjusted"]

loo <- as.data.table(b108$leave_one_component_out)
pretty_component <- c(
  ERBB2_AMP_R = "ERBB2 amplification",
  ERBB2_ACTIVATING_R = "ERBB2 activating mutation",
  MET_AMP_R = "MET amplification",
  PIK3CA_EX20_R = "PIK3CA exon 20",
  PTEN_LOF_R = "PTEN loss-of-function",
  AKT1_E17K_R = "AKT1 E17K",
  CANONICAL_PRESSING_FUSION_R = "Canonical ALK/ROS1/NTRK/RET fusion",
  KRAS_AMP_R = "KRAS amplification",
  MAP2K1_HOTSPOT_R = "MAP2K1 hotspot",
  NF1_LOF_R = "NF1 loss-of-function"
)
loo[, label := paste0("Leave out ", pretty_component[COMPONENT])]

landmark <- as.data.table(b108$landmark)
landmark[, label := ANALYSIS]

robust <- rbindlist(list(
  primary_robust[, .(label, HR, LCL95, UCL95)],
  loo[, .(label, HR, LCL95, UCL95)],
  landmark[, .(label, HR, LCL95, UCL95)]
), fill = TRUE)
robust <- robust[is.finite(HR) & is.finite(LCL95) & is.finite(UCL95)]
robust[, label := factor(label, levels = rev(unique(label)))]

p4a <- ggplot(robust, aes(x = HR, y = label)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = COL_DARK) +
  geom_errorbarh(aes(xmin = LCL95, xmax = UCL95),
                 height = 0.12, color = COL_POS, linewidth = 0.65) +
  geom_point(color = COL_POS, size = 2.4) +
  scale_x_log10(breaks = c(1, 2, 4)) +
  labs(title = "A  Robustness of adjusted rwPFS association",
       x = "Adjusted rwPFS hazard ratio (log scale)", y = NULL) +
  theme_pub(9) +
  theme(legend.position = "none")

boot_hr <- as.numeric(b108$bootstrap_hr)
boot_hr <- boot_hr[is.finite(boot_hr)]
primary_hr <- primary_robust$HR[1]

p4b <- ggplot(data.table(HR = boot_hr), aes(x = HR)) +
  geom_histogram(bins = 35, boundary = 1,
                 fill = COL_NEG, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "dotted", color = COL_DARK, linewidth = 0.8) +
  geom_vline(xintercept = primary_hr, linetype = "dashed", color = COL_POS, linewidth = 0.95) +
  labs(
    title = "B  Patient-level bootstrap distribution",
    subtitle = paste0(length(boot_hr), "/2,000 successful; median HR ",
                      sprintf("%.2f", median(boot_hr))),
    x = "Adjusted rwPFS hazard ratio",
    y = "Bootstrap iterations"
  ) +
  theme_pub(9) +
  theme(legend.position = "none")

fig4 <- p4a | p4b + plot_layout(widths = c(1.05, 0.95))
save_pub(fig4, "Figure4_robustness_validation", 11.0, 5.8)

# ============================================================
# FIGURE S1 -- former Figure 4C / 4D
# ============================================================

component_counts <- as.data.table(b108$component_counts)
component_counts[, label := fifelse(
  COMPONENT %chin% names(pretty_component),
  pretty_component[COMPONENT],
  gsub("_R$", "", COMPONENT)
)]
component_counts[, label := factor(label, levels = rev(label))]

pS1a <- ggplot(component_counts, aes(x = N_CARRIERS, y = label)) +
  geom_col(width = 0.62, fill = COL_NEG) +
  geom_text(aes(label = N_CARRIERS), hjust = -0.18, size = 3.0) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(title = "A  Carrier counts by bypass component", x = "Patients", y = NULL) +
  theme_pub(9)

rmst_diff <- as.data.table(b108$rmst)
rmst_diff[, horizon := factor(paste0(TAU_DAYS, " days"),
                              levels = c("180 days", "365 days"))]
rmst_diff[, label := sprintf("%.1f", DIFFERENCE_POS_MINUS_NEG_DAYS)]

pS1b <- ggplot(rmst_diff, aes(x = horizon, y = DIFFERENCE_POS_MINUS_NEG_DAYS)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_col(width = 0.56, fill = COL_POS) +
  geom_text(aes(label = label), vjust = 1.45, color = "white", size = 3.2) +
  labs(title = "B  Absolute rwPFS deficit",
       x = "RMST horizon", y = "Bypass+ minus bypass−, days") +
  theme_pub(9)

figS1 <- pS1a | pS1b
save_pub(figS1, "FigureS1_robustness_context", 10.5, 4.8)

# ============================================================
# FIGURE 5 MAIN -- 2 PANELS ONLY
# ============================================================

tt_models <- as.data.table(b110$model_results)
tt_models[, label := paste0(ENDPOINT, " — ", MODEL)]
tt_models[, label := factor(label, levels = rev(unique(label)))]

p5a <- ggplot(tt_models, aes(x = HR, y = label)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = COL_DARK) +
  geom_errorbarh(aes(xmin = LCL95, xmax = UCL95),
                 height = 0.12, color = COL_POS, linewidth = 0.65) +
  geom_point(color = COL_POS, size = 2.4) +
  scale_x_log10(breaks = c(1, 2, 4)) +
  labs(title = "A  Refined treatment-based validation",
       x = "Hazard ratio (log scale)", y = NULL) +
  theme_pub(8.8) +
  theme(legend.position = "none")

if (!all(c("CRC30_DAYS_R", "CRC30_EVENT_R") %in% names(ttdt))) {
  stop("B1-10 patient_data does not contain CRC30 endpoint columns.")
}

p5b <- make_km_with_risk(
  ttdt, "CRC30_DAYS_R", "CRC30_EVENT_R",
  "B  Conservative CRC TTNTD30",
  c(0, 250, 500, 750, 1000, 1250), 1250
)

fig5 <- p5a | p5b + plot_layout(widths = c(1.05, 0.95))
save_pub(fig5, "Figure5_TTNTD_validation", 11.0, 5.8)

# ============================================================
# FIGURE S2 -- former Figure 5C / 5D
# ============================================================

tt_event <- as.data.table(b110$event_summary)
tt_event[, event_rate := EVENTS / N]
tt_event[, definition := fifelse(
  grepl("^Conservative", ENDPOINT),
  "Conservative CRC-compatible",
  "Broad cancer-directed"
)]

pS2a <- ggplot(tt_event, aes(x = ENDPOINT, y = event_rate, fill = definition)) +
  geom_col(width = 0.60) +
  geom_text(aes(label = paste0(EVENTS, "/", N)), vjust = -0.32, size = 3.0) +
  scale_fill_manual(values = c(
    "Conservative CRC-compatible" = COL_TEAL,
    "Broad cancer-directed" = COL_GOLD
  )) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "A  Event completeness across TTNTD definitions",
       x = NULL, y = "Event proportion", fill = NULL) +
  theme_pub(8.5) +
  theme(axis.text.x = element_text(angle = 23, hjust = 1))

tt_median <- ttdt[, .(
  median_days = median(CRC30_DAYS_R, na.rm = TRUE)
), by = BYPASS_HIGH_CONFIDENCE_R]

tt_median[, group := fifelse(
  BYPASS_HIGH_CONFIDENCE_R == 0L,
  "High-confidence bypass−",
  "High-confidence bypass+"
)]
tt_median[, group := factor(group,
                            levels = c("High-confidence bypass−", "High-confidence bypass+"))]
tt_median[, label := sprintf("%.1f", median_days)]

pS2b <- ggplot(tt_median, aes(x = group, y = median_days, fill = group)) +
  geom_col(width = 0.56) +
  geom_text(aes(label = label), vjust = -0.30, size = 3.1) +
  scale_fill_manual(values = c(
    "High-confidence bypass−" = COL_NEG,
    "High-confidence bypass+" = COL_POS
  )) +
  labs(title = "B  Raw median conservative TTNTD30",
       x = NULL, y = "Median days", fill = NULL) +
  theme_pub(9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 12, hjust = 1))

figS2 <- pS2a | pS2b
save_pub(figS2, "FigureS2_TTNTD_endpoint_audit", 10.5, 4.8)

# ============================================================
# AUDIT
# ============================================================

expected <- file.path(fig_dir, c(
  "Figure1_study_design_and_cohort_flow.png",
  "Figure1_study_design_and_cohort_flow.tiff",
  "Figure2_high_confidence_bypass_landscape.png",
  "Figure2_high_confidence_bypass_landscape.tiff",
  "Figure3_primary_rwPFS_OS.png",
  "Figure3_primary_rwPFS_OS.tiff",
  "Figure4_robustness_validation.png",
  "Figure4_robustness_validation.tiff",
  "Figure5_TTNTD_validation.png",
  "Figure5_TTNTD_validation.tiff",
  "FigureS1_robustness_context.png",
  "FigureS1_robustness_context.tiff",
  "FigureS2_TTNTD_endpoint_audit.png",
  "FigureS2_TTNTD_endpoint_audit.tiff"
))

missing <- expected[!file.exists(expected)]
audit_file <- file.path(audit_dir, "B1_11C_final_figure_refinement_audit.txt")

writeLines(c(
  "B1-11C FINAL FIGURE REFINEMENT AUDIT",
  "====================================",
  "",
  paste0("Locked cohort N: ", nrow(dt)),
  paste0("High-confidence bypass+: ", dt[BYPASS_HIGH_CONFIDENCE_R == 1L, .N]),
  paste0("High-confidence bypass-: ", dt[BYPASS_HIGH_CONFIDENCE_R == 0L, .N]),
  "",
  "Figure 1A: black-and-white only.",
  "Figure 3: group order standardized; risk tables and HR text added; time-zero risk-count clipping corrected.",
  "Figure 4: reduced to 2 main panels.",
  "Figure 5: reduced to 2 main panels.",
  "Figure S1/S2 contain panels removed from main figures.",
  "",
  paste0("Missing outputs: ", length(missing)),
  if (length(missing) > 0L) missing else "All expected outputs exist.",
  "",
  "No new primary analysis introduced."
), audit_file)

cat("\n============================================================\n")
cat("B1-11C COMPLETE\n")
cat("============================================================\n")
cat("Missing outputs: ", length(missing), "\n", sep = "")
cat("\nAudit:\n", audit_file, "\n")
