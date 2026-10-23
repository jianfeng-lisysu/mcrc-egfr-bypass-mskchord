################################################################################
# A7 — Literature-based external consistency analysis
# (Online Resource: Supplementary Methods/Results, Table S10, Figure S4)
#
# What this script does
#   1. Holds the published hazard ratios (HRs) and 95% CIs transcribed from the
#      source publications (no patient-level data are used).
#   2. Pools them on the log scale with DerSimonian-Laird random effects
#      (metafor), with modified (truncated) Hartung-Knapp intervals as a
#      small-k sensitivity analysis.
#   3. Compares the MSK-CHORD estimate with the pooled external estimate.
#   4. Writes result tables (CSV + Word three-line tables) and draws Figure S4
#      (TIFF and PNG at 600 dpi, plus PDF).
#
# Requirements: R >= 4.1; packages metafor, ggplot2, patchwork, flextable, officer
#   install.packages(c("metafor", "ggplot2", "patchwork", "flextable", "officer"))
# Usage: outputs are written to OUT_DIR below (the A7 project folder); run
#   Rscript 21_literature_consistency_analysis.R
################################################################################

suppressPackageStartupMessages({
  library(metafor)
  library(ggplot2)
  library(patchwork)
})

project_root <- "/path/to/A7"   # set to your own directory
OUT_DIR <- file.path(project_root, "11_literature_consistency")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIG_DIR <- file.path(project_root, "05_figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
z <- qnorm(0.975)
se_from_ci <- function(lo, hi) (log(hi) - log(lo)) / (2 * z)

# ------------------------------------------------------------------------------
# 1. Source estimates (verify against the original publications before use)
# ------------------------------------------------------------------------------
# A. Gene-altered vs non-altered, within anti-EGFR-treated patients
#    MSK-CHORD values are this study's locked primary adjusted estimates.
prog <- data.frame(
  endpoint = rep(c("PFS", "OS"), each = 4),
  study = rep(c("MSK-CHORD (this study)", "Emmett 2024 (Flatiron-FMI)",
                "PANDA (Pietrantonio 2023)", "PRESSING2 (Randon 2022)"), 2),
  source = rep(c("This study: adjusted Cox, rwPFS / OS",
                 "ESMO Real World Data Digit Oncol 2024;4:100036 - first-line EGFR mAb, adjusted rwPFS / rwOS",
                 "Eur J Cancer 2023;195:113396 - PFS / OS, gene-altered vs hyperselected",
                 "JCO Precis Oncol 2022;6:e2200037 - multivariable PFS / OS, whole cohort"), 2),
  hr = c(2.33, 1.52, 2.08, 4.19,   1.47, 1.64, 1.82, 2.98),
  lo = c(1.45, 1.14, 1.43, 2.58,   0.88, 1.14, 1.23, 1.49),
  hi = c(3.76, 2.02, 3.03, 6.79,   2.48, 2.36, 2.69, 5.96),
  stringsAsFactors = FALSE)
prog$yi  <- log(prog$hr)
prog$sei <- se_from_ci(prog$lo, prog$hi)

# B. Treatment-by-alteration interaction
#    External studies report anti-EGFR vs bevacizumab HRs separately in gene-altered
#    (a_*) and non-altered (n_*) strata; the interaction is approximated as the
#    ratio a/n with var(log ratio) = var(log a) + var(log n) (disjoint strata).
strata <- data.frame(
  endpoint = c("PFS", "PFS", "PFS", "OS", "OS", "OS"),
  study = c("Emmett 2024 (Flatiron-FMI)", "PARADIGM, left-sided (Shitara 2024)", "FIRE-3, left-sided (Weiss 2025)",
            "Emmett 2024 (Flatiron-FMI)", "PARADIGM, left-sided (Shitara 2024)", "FIRE-3, left-sided (Weiss 2025)"),
  source = c("Emmett 2024: EGFR mAb vs bev, adjusted", "Nat Med 2024;30:730-739, left-sided subset",
             "Eur J Cancer 2025;221:115399, left-sided, PRESSING-3 panel",
             "Emmett 2024: EGFR mAb vs bev, adjusted", "Nat Med 2024;30:730-739, left-sided subset",
             "Eur J Cancer 2025;221:115399, left-sided, PRESSING-3 panel"),
  a_hr = c(1.37, 1.45, 1.22, 1.32, 1.08, 1.10),
  a_lo = c(1.05, 0.94, 0.66, 0.96, 0.71, 0.59),
  a_hi = c(1.78, 2.23, 2.26, 1.83, 1.64, 2.07),
  n_hr = c(0.93, 0.91, 1.02, 0.82, 0.76, 0.70),
  n_lo = c(0.77, 0.73, 0.67, 0.64, 0.61, 0.43),
  n_hi = c(1.12, 1.13, 1.57, 1.05, 0.95, 1.12),
  stringsAsFactors = FALSE)
strata$yi  <- log(strata$a_hr) - log(strata$n_hr)
strata$sei <- sqrt(se_from_ci(strata$a_lo, strata$a_hi)^2 + se_from_ci(strata$n_lo, strata$n_hi)^2)

# MSK-CHORD: directly estimated overlap-weighted interaction (model M1), rwPFS only
msk_int <- data.frame(endpoint = "PFS", study = "MSK-CHORD (this study)",
                      source = "This study: overlap-weighted M1 interaction",
                      yi = log(3.08), sei = se_from_ci(1.45, 6.55), stringsAsFactors = FALSE)
inter <- rbind(msk_int, strata[, c("endpoint", "study", "source", "yi", "sei")])

# Sensitivity source: 7-trial IPD pooled analysis (ASCO 2026 abstract 3550), OS.
# Point estimates per JCO abstract (aHR 0.83 / 1.16); 95% CIs available only from the
# meeting summary, which reported 0.84 (0.73-0.97) and 1.17 (0.76-1.76); the latter are used.
ipd <- data.frame(endpoint = "OS", study = "Germani IPD, 7 trials (sensitivity)",
                  source = "J Clin Oncol 2026;44(16 suppl):3550",
                  yi = log(1.17) - log(0.84),
                  sei = sqrt(se_from_ci(0.76, 1.76)^2 + se_from_ci(0.73, 0.97)^2))

# ------------------------------------------------------------------------------
# 2. Pooling helpers
# ------------------------------------------------------------------------------
pool <- function(d, label) {
  dl <- rma(yi = d$yi, sei = d$sei, method = "DL")
  hk <- rma(yi = d$yi, sei = d$sei, method = "DL", test = "adhoc")   # modified (truncated) Knapp-Hartung
  data.frame(analysis = label, k = nrow(d),
             est = exp(dl$b[1]), lo = exp(dl$ci.lb), hi = exp(dl$ci.ub), p = dl$pval,
             I2 = dl$I2, tau2 = dl$tau2,
             hksj_lo = exp(hk$ci.lb), hksj_hi = exp(hk$ci.ub),
             log_est = dl$b[1], log_se = dl$se)
}
msk_vs_ext <- function(msk, ext, label) {
  d <- msk$yi - ext$log_est; s <- sqrt(msk$sei^2 + ext$log_se^2)
  data.frame(comparison = label, ratio = exp(d), lo = exp(d - z * s), hi = exp(d + z * s), p = 2 * pnorm(-abs(d / s)))
}
is_msk <- function(d) grepl("MSK-CHORD", d$study)

# ------------------------------------------------------------------------------
# 3. Analyses
# ------------------------------------------------------------------------------
res <- list(); cmp <- list()
for (ep in c("PFS", "OS")) {
  d <- prog[prog$endpoint == ep, ]
  e <- pool(d[!is_msk(d), ], paste("Altered vs non-altered,", ep, "- external only"))
  res[[length(res) + 1]] <- e
  res[[length(res) + 1]] <- pool(d, paste("Altered vs non-altered,", ep, "- all incl. MSK-CHORD"))
  cmp[[length(cmp) + 1]] <- msk_vs_ext(d[is_msk(d), ], e, paste("MSK-CHORD vs external, altered vs non-altered,", ep))
}
res[[length(res) + 1]] <- pool(prog[prog$endpoint == "PFS" & !grepl("PRESSING2", prog$study), ],
                               "Altered vs non-altered, PFS - all excl. PRESSING2")
for (ep in c("PFS", "OS")) {
  d <- inter[inter$endpoint == ep, ]
  e <- pool(d[!is_msk(d), ], paste("Interaction,", ep, "- external only"))
  res[[length(res) + 1]] <- e
  if (any(is_msk(d))) {
    res[[length(res) + 1]] <- pool(d, paste("Interaction,", ep, "- all incl. MSK-CHORD"))
    cmp[[length(cmp) + 1]] <- msk_vs_ext(d[is_msk(d), ], e, paste("MSK-CHORD vs external, interaction,", ep))
  }
}
res[[length(res) + 1]] <- pool(rbind(inter[inter$endpoint == "OS" & grepl("Emmett", inter$study), c("yi", "sei")], ipd[, c("yi", "sei")]),
                               "Interaction, OS - Emmett + IPD (sensitivity; IPD may overlap FIRE-3/PARADIGM)")
res <- do.call(rbind, res); cmp <- do.call(rbind, cmp)

per_study <- rbind(
  transform(prog[, c("endpoint", "study", "source")], analysis = "Altered vs non-altered within anti-EGFR",
            est = prog$hr, lo = prog$lo, hi = prog$hi),   # as reported
  transform(rbind(inter, ipd)[, c("endpoint", "study", "source")], analysis = "Interaction (anti-EGFR vs bev HR ratio, altered/non-altered)",
            est = exp(c(inter$yi, ipd$yi)), lo = exp(c(inter$yi, ipd$yi) - z * c(inter$sei, ipd$sei)),
            hi = exp(c(inter$yi, ipd$yi) + z * c(inter$sei, ipd$sei))))

fmt <- function(x) sprintf("%.2f", x)
cat("\n== Per-study estimates ==\n"); print(transform(per_study[, c("analysis", "endpoint", "study")], est_95CI = paste0(fmt(per_study$est), " (", fmt(per_study$lo), "-", fmt(per_study$hi), ")")), row.names = FALSE)
cat("\n== Pooled (DerSimonian-Laird) with modified Hartung-Knapp CI ==\n")
print(data.frame(analysis = res$analysis, k = res$k,
                 DL = paste0(fmt(res$est), " (", fmt(res$lo), "-", fmt(res$hi), ")"),
                 I2 = sprintf("%.0f%%", res$I2), HKSJ_CI = paste0(fmt(res$hksj_lo), "-", fmt(res$hksj_hi))), row.names = FALSE)
cat("\n== MSK-CHORD vs pooled external ==\n")
print(transform(cmp, ratio = fmt(ratio), lo = fmt(lo), hi = fmt(hi), p = sprintf("%.3f", p)), row.names = FALSE)

write.csv(per_study, file.path(OUT_DIR, "21_per_study_estimates.csv"), row.names = FALSE)
write.csv(res[, setdiff(names(res), c("log_est", "log_se"))], file.path(OUT_DIR, "21_pooled_estimates.csv"), row.names = FALSE)
write.csv(cmp, file.path(OUT_DIR, "21_msk_vs_external.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Figure S4
# ------------------------------------------------------------------------------
ORANGE <- "#D55E00"; INK <- "#1a1a1a"; MUTED <- "#6b6b6b"
build_panel <- function(blocks, xlab, hdr) {
  rows <- list(); y <- 0
  for (b in blocks) {
    rows[[length(rows) + 1]] <- data.frame(y = y, label = b$title, est = NA, lo = NA, hi = NA, kind = "header"); y <- y - 1
    for (i in seq_len(nrow(b$d))) {
      rep_ok <- all(c("hr", "lo", "hi") %in% names(b$d))   # reported CI available (panel a)
      rows[[length(rows) + 1]] <- data.frame(y = y, label = b$d$study[i], est = exp(b$d$yi[i]),
        lo = if (rep_ok) b$d$lo[i] else exp(b$d$yi[i] - z * b$d$sei[i]),
        hi = if (rep_ok) b$d$hi[i] else exp(b$d$yi[i] + z * b$d$sei[i]),
        kind = ifelse(grepl("MSK", b$d$study[i]), "msk", "ext")); y <- y - 1
    }
    for (pp in b$pooled) {
      rows[[length(rows) + 1]] <- data.frame(y = y, label = sprintf("%s  (I2 = %.0f%%)", pp$lab, pp$r$I2),
        est = pp$r$est, lo = pp$r$lo, hi = pp$r$hi, kind = pp$kind); y <- y - 1
    }
    y <- y - 0.6
  }
  df <- do.call(rbind, rows)
  df$txt <- ifelse(is.na(df$est), "", sprintf("%.2f (%.2f-%.2f)", df$est, df$lo, df$hi))
  col <- c(header = INK, msk = ORANGE, ext = INK, pext = INK, pall = MUTED)
  dia <- do.call(rbind, lapply(which(df$kind %in% c("pext", "pall")), function(i)
    data.frame(id = i, x = c(df$lo[i], df$est[i], df$hi[i], df$est[i]), yy = df$y[i] + c(0, .32, 0, -.32), kind = df$kind[i])))
  pt <- df[df$kind %in% c("msk", "ext"), ]
  ggplot() +
    geom_vline(xintercept = 1, linetype = "dotted", colour = MUTED) +
    geom_segment(data = pt, aes(x = lo, xend = hi, y = y, yend = y, colour = kind), linewidth = 0.6) +
    geom_point(data = pt, aes(x = est, y = y, colour = kind, shape = kind), size = 2.2) +
    geom_polygon(data = dia, aes(x = x, y = yy, group = id, fill = kind)) +
    geom_text(data = df, aes(x = 11, y = y, label = txt, colour = kind,
              fontface = ifelse(kind %in% c("pext", "pall"), "bold", "plain")), hjust = 0, size = 2.7) +
    annotate("text", x = 11, y = 1, label = hdr, hjust = 0, fontface = "bold", size = 2.8) +
    scale_colour_manual(values = col, guide = "none") + scale_fill_manual(values = col, guide = "none") +
    scale_shape_manual(values = c(msk = 18, ext = 15), guide = "none") +
    scale_x_log10(breaks = c(0.5, 1, 2, 4, 8), labels = c("0.5", "1", "2", "4", "8")) +
    scale_y_continuous(breaks = df$y, labels = df$label, expand = expansion(add = 0.6)) +
    coord_cartesian(xlim = c(0.5, 8), clip = "off") +
    labs(x = xlab, y = NULL) +
    theme_classic(base_size = 8.5) +
    theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.y = element_text(colour = INK),
          plot.margin = margin(14, 95, 4, 4), plot.title = element_text(face = "bold", size = 9.5))
}
P <- function(ep, kind) { d <- if (kind == "prog") prog[prog$endpoint == ep, ] else inter[inter$endpoint == ep, ]; d }
R <- function(lbl) res[res$analysis == lbl, ]
pA <- build_panel(list(
  list(title = "Progression-free survival", d = P("PFS", "prog"),
       pooled = list(list(lab = "Pooled, external studies only", r = R("Altered vs non-altered, PFS - external only"), kind = "pext"),
                     list(lab = "Pooled, all studies", r = R("Altered vs non-altered, PFS - all incl. MSK-CHORD"), kind = "pall"))),
  list(title = "Overall survival", d = P("OS", "prog"),
       pooled = list(list(lab = "Pooled, external studies only", r = R("Altered vs non-altered, OS - external only"), kind = "pext"),
                     list(lab = "Pooled, all studies", r = R("Altered vs non-altered, OS - all incl. MSK-CHORD"), kind = "pall")))),
  "Hazard ratio, altered vs non-altered (log scale)", "HR (95% CI)") +
  ggtitle("a  Altered vs non-altered, within anti-EGFR-treated patients")
pB <- build_panel(list(
  list(title = "Progression-free survival", d = P("PFS", "int"),
       pooled = list(list(lab = "Pooled, external studies only", r = R("Interaction, PFS - external only"), kind = "pext"),
                     list(lab = "Pooled, all studies", r = R("Interaction, PFS - all incl. MSK-CHORD"), kind = "pall"))),
  list(title = "Overall survival", d = P("OS", "int"),
       pooled = list(list(lab = "Pooled, external studies only", r = R("Interaction, OS - external only"), kind = "pext")))),
  "Ratio of anti-EGFR vs bevacizumab HRs, altered / non-altered (log scale)", "Ratio (95% CI)") +
  ggtitle("b  Treatment-by-alteration interaction, anti-EGFR vs bevacizumab")
fig <- pA / pB + plot_annotation(caption = paste(
  "Random-effects (DerSimonian-Laird) pooling. MSK-CHORD (this study, orange) shown for comparison;",
  "no MSK-CHORD OS interaction was estimated. Alteration panels differ across studies; see Table S10.", sep = "\n"),
  theme = theme(plot.caption = element_text(hjust = 0, size = 7, colour = MUTED)))
suppressWarnings(ggsave(file.path(FIG_DIR, "FigureS4_literature_consistency.tiff"), fig, width = 174, height = 210, units = "mm", dpi = 600, compression = "lzw"))
suppressWarnings(ggsave(file.path(FIG_DIR, "FigureS4_literature_consistency.png"), fig, width = 174, height = 210, units = "mm", dpi = 600, bg = "white"))
suppressWarnings(ggsave(file.path(FIG_DIR, "FigureS4_literature_consistency.pdf"), fig, width = 174, height = 210, units = "mm"))

# ------------------------------------------------------------------------------
# 5. Word three-line tables (top and bottom rules 1.5 pt, header rule 0.75 pt)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(flextable); library(officer) })
ci <- function(e, l, h) sprintf("%.2f (%.2f-%.2f)", e, l, h)
pv <- function(p) ifelse(p < 0.001, "<.001", sub("^0", "", sprintf("%.3f", p)))
three_line <- function(df) {
  ft <- flextable(df)
  ft <- border_remove(ft)
  ft <- hline_top(ft, part = "header", border = fp_border(width = 1.5))
  ft <- hline_bottom(ft, part = "header", border = fp_border(width = 0.75))
  ft <- hline_bottom(ft, part = "body", border = fp_border(width = 1.5))
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "left", part = "all")
  ft <- padding(ft, padding.top = 2, padding.bottom = 2, part = "all")
  set_table_properties(ft, layout = "autofit", width = 1)
}
tA <- data.frame(Analysis = per_study$analysis, Endpoint = per_study$endpoint, Study = per_study$study,
                 `Estimate (95% CI)` = ci(per_study$est, per_study$lo, per_study$hi), Source = per_study$source, check.names = FALSE)
tB <- data.frame(Analysis = res$analysis, k = res$k, `Pooled estimate (95% CI)` = ci(res$est, res$lo, res$hi),
                 P = pv(res$p), `I2, %` = sprintf("%.0f", res$I2), `tau2` = sprintf("%.3f", res$tau2),
                 `Modified Hartung-Knapp 95% CI` = sprintf("%.2f-%.2f", res$hksj_lo, res$hksj_hi), check.names = FALSE)
tC <- data.frame(Comparison = cmp$comparison, `Ratio (95% CI)` = ci(cmp$ratio, cmp$lo, cmp$hi), P = pv(cmp$p), check.names = FALSE)
note_txt <- paste("Estimates are hazard ratios (gene-altered vs non-altered within anti-EGFR-treated patients) or ratios of",
  "anti-EGFR-versus-bevacizumab hazard ratios (gene-altered / non-altered). Pooling used DerSimonian-Laird random effects.",
  "CI, confidence interval; I2, heterogeneity; tau2, between-study variance.")
doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "landscape")))
add_tab <- function(doc, cap, df) {
  doc <- body_add_par(doc, cap, style = "Normal")
  doc <- body_add_flextable(doc, three_line(df))
  doc <- body_add_par(doc, note_txt, style = "Normal")
  body_add_par(doc, "", style = "Normal")
}
doc <- add_tab(doc, "Table A. Per-study estimates used in the literature-based external consistency analysis", tA)
doc <- add_tab(doc, "Table B. Pooled estimates", tB)
doc <- add_tab(doc, "Table C. MSK-CHORD estimate versus pooled external estimate", tC)
print(doc, target = file.path(OUT_DIR, "21_tables_three_line.docx"))
cat("\nTables written to:", normalizePath(OUT_DIR), "\nFigure S4 written to:", normalizePath(FIG_DIR), "\n")
sessionInfo()
