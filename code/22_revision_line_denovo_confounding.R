# ==============================================================================
# A7 | Script 22 | BMC Gastroenterology revision (round 1) analyses
# ------------------------------------------------------------------------------
# Version 2 (2026-11-24): imaging window fixed to +/-120 days; treatment-history
# proxies renamed and restricted to post-stage-IV exposure; ECOG 180-day model
# matched to model D; M4 weights re-estimated; M3 full-pipeline bootstrap;
# within-stratum comparator diagnostics; interaction E-value removed.
# Version 4 (2026-11-24): bootstrap fit-quality diagnostics; zero-variance SMD edge case;
# assertion that all 188 overlapping patients share the index date.
# Version 3 (2026-11-24): proportional-hazards term diagnostics and panel-stratified
# refits (Table S12b); per-covariate balance in the bypass-positive comparator
# stratum (Table S13c); audit of propensity-score-incomplete comparator patients.
#
# Purpose: respond to reviewer/editor comments on
#   (1) treatment-line heterogeneity (late-line anti-EGFR use),
#   (2) de novo (synchronous) versus metachronous stage IV disease,
#   (3) residual confounding in the primary and comparator analyses.
#
# All analyses are post hoc revision analyses. The locked exposure definition
# (BYPASS_HIGH_CONFIDENCE_R), endpoints and primary covariate set are NOT changed.
#
# Inputs (already in the A7 folder):
#   03_intermediate/A7_B1_locked_analysis.rds          locked anti-EGFR cohort (n=191)
#   04_results/B1_14B_exploratory_comparator_rwPFS.rds  comparator cohort (n=319) + frozen weights
#   01_raw_data/msk_chord_2024.tar                      public cBioPortal release (timelines)
# Outputs (all written inside the A7 folder):
#   12_revision_analyses/22_*.csv (Tables S11-S14 source data)
#   12_revision_analyses/22_revision_tables_three_line.docx
#   12_revision_analyses/figures/FigureS6_revision_sensitivity_forest.tiff / .png (600 dpi) / .pdf
# Original locked inputs are only READ; nothing in 02_code, 03_intermediate,
# 04_results, 05_figures or 07_tables is modified.
#
# Requirements: R >= 4.1; packages data.table, survival, ggplot2, flextable, officer
# ==============================================================================

A7 <- "/path/to/A7"   # set to your own directory
OUT_DIR <- file.path(A7, "12_revision_analyses")
FIG_DIR <- file.path(OUT_DIR, "figures")   # revision figures kept separate from 05_figures
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0. Read locked data and the public raw timelines
# ------------------------------------------------------------------------------
x  <- as.data.table(readRDS(file.path(A7, "03_intermediate", "A7_B1_locked_analysis.rds")))
cmp_obj <- readRDS(file.path(A7, "04_results", "B1_14B_exploratory_comparator_rwPFS.rds"))
dw <- as.data.table(cmp_obj$analysis_data)
stopifnot(nrow(x) == 191, uniqueN(x$PATIENT_ID) == 191, sum(x$BYPASS_HIGH_CONFIDENCE_R) == 23,
          sum(x$RWPFS_EVENT_LOCKED_R) == 183, nrow(dw) == 319, uniqueN(dw$PATIENT_ID) == 319)

# index dates must agree for patients present in both cohorts
chk <- merge(x[, .(PATIENT_ID, a = FIRST_BIOLOGIC_DAY_R)], dw[, .(PATIENT_ID, b = as.numeric(FIRST_BIOLOGIC_DAY_R))])
stopifnot(nrow(chk) == 188, all(chk$a == chk$b))   # all 188 overlapping anti-EGFR patients share the same index date
# anti-EGFR patients absent from the comparator: excluded from the propensity-score-complete population
fe_audit <- as.data.table(readRDS(file.path(A7, "04_results", "B1_14_bevacizumab_comparator_feasibility.rds"))$patient_audit)
not_in_cmp <- fe_audit[FIRST_BIOLOGIC_R == "anti-EGFR" & !PATIENT_ID %in% dw$PATIENT_ID]
stopifnot(nrow(not_in_cmp) == 3, all(not_in_cmp$PS_COMPLETE_R == FALSE), all(not_in_cmp$PATIENT_ID %in% x$PATIENT_ID))
ps_incomplete <- fe_audit[PS_COMPLETE_R == FALSE, .N, by = .(FIRST_BIOLOGIC_R, BYPASS_HIGH_CONFIDENCE_R)]

raw_dir <- file.path(tempdir(), "a7_raw")
dir.create(raw_dir, showWarnings = FALSE)
need <- paste0("msk_chord_2024/", c("data_timeline_diagnosis.txt",
                                    "data_timeline_tumor_sites.txt",
                                    "data_timeline_performance_status.txt",
                                    "data_timeline_treatment.txt",
                                    "data_clinical_patient.txt"))
untar(file.path(A7, "01_raw_data", "msk_chord_2024.tar"), files = need, exdir = raw_dir)
rd <- function(f) fread(file.path(raw_dir, "msk_chord_2024", f), showProgress = FALSE)
dx <- rd("data_timeline_diagnosis.txt")
tx <- rd("data_timeline_treatment.txt")
pt <- fread(file.path(raw_dir, "msk_chord_2024", "data_clinical_patient.txt"), skip = 4, showProgress = FALSE)

ids_all <- unique(c(x$PATIENT_ID, dw$PATIENT_ID))
index_tab <- unique(rbind(
  x[, .(PATIENT_ID, INDEX_DAY = FIRST_BIOLOGIC_DAY_R, S4_DAY = EXPLICIT_STAGE4_CRC_DAY_R)],
  dw[!PATIENT_ID %in% x$PATIENT_ID, .(PATIENT_ID, INDEX_DAY = as.numeric(FIRST_BIOLOGIC_DAY_R),
                                      S4_DAY = as.numeric(EXPLICIT_STAGE4_CRC_DAY_R))]), by = "PATIENT_ID")

# ------------------------------------------------------------------------------
# 1. Stage at initial diagnosis (tumour registry) and radiology corroboration
#    Registry: first primary colorectal record (ICD-O topography C18-C20, and
#    C21.8 "rectum, other parts"); stage IV at diagnosis = SEER summary stage
#    "Distant ..." or AJCC clinical stage group 4. NA if neither is recorded.
#    Radiology: radiology-report NLP tumour-site timeline (rows exist only for
#    reports in which a tumour site was extracted). A patient is corroborated
#    when a distant organ site (liver, lung, adrenal, bone, CNS, pleura,
#    reproductive organs) is recorded within +/-120 days of the registry
#    diagnosis date and before the index biologic. Intra-abdominal, lymph-node
#    and "Other" sites are not counted (conservative).
# ------------------------------------------------------------------------------
crc <- dx[PATIENT_ID %in% ids_all & SUBTYPE == "Primary" & grepl("C18|C19|C20|C218", DX_DESCRIPTION)]
crc <- crc[order(PATIENT_ID, START_DATE)][, .SD[1], by = PATIENT_ID]
crc[, SUMMARY := trimws(SUMMARY)]
crc[, CG := trimws(as.character(CLINICAL_GROUP))]
crc[, STAGE4_AT_DX_R := fifelse(grepl("^Distant", SUMMARY) | grepl("^4", CG), 1L,
                        fifelse((SUMMARY == "" | is.na(SUMMARY) | grepl("Unknown|N/A", SUMMARY)) &
                                (CG == "" | is.na(CG) | CG %in% c("99", "NA")), NA_integer_, 0L))]
crc[, STAGE4_SOURCE_R := fifelse(grepl("^Distant", SUMMARY), "SEER summary stage distant",
                         fifelse(grepl("^4", CG), "AJCC clinical stage group 4 only", "Not stage IV / not recorded"))]
denovo <- crc[, .(PATIENT_ID, DX_DAY = START_DATE, STAGE4_AT_DX_R, STAGE4_SOURCE_R)]

ts <- rd("data_timeline_tumor_sites.txt")[PATIENT_ID %in% ids_all]
dist_sites <- c("Liver", "Lung", "Adrenal Glands", "Bone", "CNS/Brain", "Pleura", "Reproductive Organs")
tsw <- merge(ts, merge(denovo[, .(PATIENT_ID, DX_DAY)], index_tab[, .(PATIENT_ID, INDEX_DAY)]), by = "PATIENT_ID")
tsw <- tsw[!is.na(START_DATE) & START_DATE >= DX_DAY - 120 & START_DATE <= DX_DAY + 120 & START_DATE < INDEX_DAY]
img <- tsw[, .(ANY_SITE_WIN = TRUE, DIST_WIN = any(TUMOR_SITE %in% dist_sites)), by = PATIENT_ID]
denovo <- merge(denovo, img, by = "PATIENT_ID", all.x = TRUE)
lev_img <- c("Distant organ site within 120 days of diagnosis",
             "Tumour-site records within 120 days, no distant organ site",
             "No tumour-site record within 120 days (not assessable)")
denovo[, IMAGING_CLASS_R := factor(fifelse(DIST_WIN %in% TRUE, lev_img[1],
                                   fifelse(ANY_SITE_WIN %in% TRUE, lev_img[2], lev_img[3])), levels = lev_img)]
denovo[, STAGE4_IMG_CORROB_R := as.integer(STAGE4_AT_DX_R %in% 1L & IMAGING_CLASS_R == lev_img[1])]
denovo[, c("DX_DAY", "ANY_SITE_WIN", "DIST_WIN") := NULL]

# ------------------------------------------------------------------------------
# 2. Treatment-history proxies at biologic initiation (NOT chart-verified lines)
#    (a) Prior cytotoxic backbone exposure (from the locked data): none;
#        fluoropyrimidine only or one of oxaliplatin/irinotecan; both.
#    (b) Timeline-derived prior drug-class exposure count: number of distinct
#        line-defining drug classes (oxaliplatin, irinotecan, regorafenib,
#        trifluridine-tipiracil, investigational, immunotherapy, other agents;
#        fluoropyrimidines and anti-angiogenic antibodies excluded) started
#        after (stage IV diagnosis - 30 days) and >28 days before the index
#        biologic, plus one if a fluoropyrimidine-only period preceded the first
#        such class by >28 days. This is an exposure count, not a regimen
#        sequence (e.g., first-line FOLFOXIRI counts as two classes).
#        No recorded prior exposure is coded 0 (absence of records).
# ------------------------------------------------------------------------------
fp_like  <- c("FLUOROURACIL", "CAPECITABINE", "FLOXURIDINE", "LEUCOVORIN", "LEVOLEUCOVORIN")
not_line <- c(fp_like, "BEVACIZUMAB", "ZIV-AFLIBERCEPT", "RAMUCIRUMAB")
excl_sub <- c("Bone Treatment", "Hormone")
line_class <- function(agent, subtype) {
  fifelse(agent == "OXALIPLATIN", "OX",
  fifelse(agent == "IRINOTECAN", "IRI",
  fifelse(agent == "REGORAFENIB", "REG",
  fifelse(agent == "TIPIRACIL-TRIFLURIDINE", "TAS",
  fifelse(subtype == "Investigational" | agent == "INVESTIGATIONAL", "INV",
  fifelse(subtype == "Immuno", "IO",
  fifelse(agent %in% c("CETUXIMAB", "PANITUMUMAB"), "EGFR",
  fifelse(grepl("Y-90|YTTRIUM", agent), "NONE",
          paste0("OTHER_", agent)))))))))
}
tt <- merge(tx[PATIENT_ID %in% ids_all & !SUBTYPE %in% excl_sub & !is.na(START_DATE)], index_tab, by = "PATIENT_ID")
tt <- tt[START_DATE < INDEX_DAY - 28 & (is.na(S4_DAY) | START_DATE >= S4_DAY - 30)]
tt[, CLASS := fifelse(AGENT %in% not_line, ifelse(AGENT %in% fp_like, "FP", "ANTIVEGF"), line_class(AGENT, SUBTYPE))]
tt <- tt[CLASS != "NONE"]
cls <- tt[, {
  lc <- unique(CLASS[!CLASS %in% c("FP", "ANTIVEGF")])
  first_start <- if (length(lc)) min(START_DATE[!CLASS %in% c("FP", "ANTIVEGF")]) else Inf
  .(PRIOR_CLASS_N_R = length(lc) + as.integer(any(CLASS == "FP" & START_DATE < first_start - 28)))
}, by = PATIENT_ID]
cls <- merge(index_tab[, .(PATIENT_ID)], cls, all.x = TRUE)
cls[is.na(PRIOR_CLASS_N_R), PRIOR_CLASS_N_R := 0L]
cls[, PRIOR_CLASS_CAT_R := factor(pmin(PRIOR_CLASS_N_R, 2L), levels = 0:2,
                                  labels = c("No recorded prior class", "One prior class", "Two or more prior classes"))]

prior_out <- pt[PATIENT_ID %in% ids_all, .(PATIENT_ID, PRIOR_MED_TO_MSK)]

add_vars <- function(d) {
  d <- merge(d, denovo, by = "PATIENT_ID", all.x = TRUE)
  d <- merge(d, cls, by = "PATIENT_ID", all.x = TRUE)
  d <- merge(d, prior_out, by = "PATIENT_ID", all.x = TRUE)
  d[is.na(STAGE4_IMG_CORROB_R), STAGE4_IMG_CORROB_R := 0L]
  d[, STAGE4_AT_DX_F_R := factor(fifelse(STAGE4_AT_DX_R %in% 1L, "Yes", "No or unknown"), levels = c("Yes", "No or unknown"))]
  d[, OUTSIDE_PRIOR_R := factor(fifelse(PRIOR_MED_TO_MSK == "No prior medications", "No",
                                fifelse(PRIOR_MED_TO_MSK == "Prior medications to MSK", "Yes", "Unknown")),
                                levels = c("No", "Yes", "Unknown"))]
  stopifnot(all(is.na(d$ECOG_LAG_DAYS_R) | d$ECOG_LAG_DAYS_R >= 0))
  d[, ECOG_CAT_R := factor(fifelse(is.na(ECOG_LATEST_PRE_R) | is.na(ECOG_LAG_DAYS_R) | ECOG_LAG_DAYS_R > 180,
                                   "Not recorded within 180 days",
                           fifelse(ECOG_LATEST_PRE_R == 0, "0",
                           fifelse(ECOG_LATEST_PRE_R == 1, "1", "2 or higher"))),
                           levels = c("0", "1", "2 or higher", "Not recorded within 180 days"))]
  d
}
x <- add_vars(x)
dw <- add_vars(dw)
x[, BACKBONE_EXPO_R := factor(fifelse(PRIOR_CYTOTOXIC_PATTERN_R == "no prior major cytotoxic", "None",
                              fifelse(PRIOR_CYTOTOXIC_PATTERN_R == "prior oxaliplatin + irinotecan",
                                      "Both oxaliplatin and irinotecan",
                                      "Fluoropyrimidine only or one of oxaliplatin/irinotecan")),
                              levels = c("None", "Fluoropyrimidine only or one of oxaliplatin/irinotecan",
                                         "Both oxaliplatin and irinotecan"))]
x[, BOTH_BACKBONES_R := as.integer(BACKBONE_EXPO_R == "Both oxaliplatin and irinotecan")]
x[, COMBO_AT_INDEX_R := factor(fifelse(ACTIVE_BACKBONE_AT_INDEX_R == "no major cytotoxic active",
                                       "Without concurrent cytotoxic therapy", "With concurrent cytotoxic therapy"),
                               levels = c("With concurrent cytotoxic therapy", "Without concurrent cytotoxic therapy"))]

# ECOG from the raw performance-status timeline (script 13 definition):
# nearest pre-index record within the window, worst value on that day, 0 vs >=1.
ps_raw <- rd("data_timeline_performance_status.txt")[PATIENT_ID %in% x$PATIENT_ID]
ps_raw[, ECOG := suppressWarnings(as.numeric(ECOG))]
ps_raw <- ps_raw[is.finite(START_DATE) & is.finite(ECOG) & ECOG >= 0 & ECOG <= 4]
ecog_win <- function(win) {
  e <- merge(ps_raw, x[, .(PATIENT_ID, IDX = FIRST_BIOLOGIC_DAY_R)], by = "PATIENT_ID")[START_DATE <= IDX & START_DATE >= IDX - win]
  e[, .SD[START_DATE == max(START_DATE)], by = PATIENT_ID][, .(E = max(ECOG)), by = PATIENT_ID]
}
x <- merge(x, setnames(ecog_win(90), "E", "ECOG90_R"), by = "PATIENT_ID", all.x = TRUE)
x <- merge(x, setnames(ecog_win(180), "E", "ECOG180_R"), by = "PATIENT_ID", all.x = TRUE)
zc <- tolower(x$PRIOR_CYTOTOXIC_PATTERN_R); bb <- toupper(x$ACTIVE_BACKBONE_AT_INDEX_R)
x[, PRIOR_MAJOR_N_F_R := factor(fcase(grepl("oxaliplatin.*irinotecan", zc), 2L, grepl("one major", zc), 1L, default = 0L))]
x[, ACTIVE_IRI_F_R := factor(as.integer(grepl("IRINOTECAN|FOLFIRI|CAPIRI|FOLFOXIRI", bb)))]
x[, ACTIVE_OX_F_R := factor(as.integer(grepl("OXALIPLATIN|FOLFOX|CAPOX|FOLFOXIRI", bb)))]
x[, PANEL_OLD_F_R := factor(PANEL_OLD_R)]
x[, SAMPLE_METASTATIC_F_R := factor(as.integer(grepl("METAST", toupper(SAMPLE_TYPE))))]
x[, LOG_SPECIMEN_AGE_R := fifelse(SPECIMEN_TO_ANTI_EGFR_DAYS_R >= 0, log1p(SPECIMEN_TO_ANTI_EGFR_DAYS_R), NA_real_)]
x[, MET_SITE_COUNT_NUM_R := as.numeric(MET_SITE_COUNT_90D_R)]
x[, ECOG90_F_R := factor(as.integer(ECOG90_R >= 1))]
x[, ECOG180_F_R := factor(as.integer(ECOG180_R >= 1))]
cc_c <- !is.na(x$LOG_SPECIMEN_AGE_R) & !is.na(x$MET_SITE_COUNT_NUM_R)
x[, CC90_R := as.integer(cc_c & !is.na(ECOG90_R))]
x[, CC180_R := as.integer(cc_c & !is.na(ECOG180_R))]
modelC_cov <- paste("GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + LOG_STAGE4_TO_ANTI_EGFR_R +",
                    "PRIOR_MAJOR_N_F_R + ACTIVE_IRI_F_R + ACTIVE_OX_F_R + PANEL_OLD_F_R +",
                    "SAMPLE_METASTATIC_F_R + LOG_SPECIMEN_AGE_R + MET_SITE_COUNT_NUM_R")

# ------------------------------------------------------------------------------
# 3. Descriptive table (Table S11)
# ------------------------------------------------------------------------------
npct <- function(v, lev) sprintf("%d (%.1f)", sum(v == lev, na.rm = TRUE), 100 * sum(v == lev, na.rm = TRUE) / length(v))
medi <- function(v) sprintf("%.0f (%.0f-%.0f)", median(v, na.rm = TRUE),
                            quantile(v, .25, na.rm = TRUE), quantile(v, .75, na.rm = TRUE))
desc_rows <- function(d) {
  out <- list()
  add <- function(lab, val) out[[length(out) + 1]] <<- c(lab, val)
  add("Stage IV at initial diagnosis (registry)", npct(d$STAGE4_AT_DX_F_R, "Yes"))
  for (l in c("SEER summary stage distant", "AJCC clinical stage group 4 only")) add(paste0("  Source: ", l), npct(d$STAGE4_SOURCE_R, l))
  add("Stage I-III or not recorded at initial diagnosis", npct(d$STAGE4_AT_DX_F_R, "No or unknown"))
  for (l in levels(d$IMAGING_CLASS_R)) add(paste0("Radiology: ", l), npct(d$IMAGING_CLASS_R, l))
  for (l in levels(d$BACKBONE_EXPO_R)) add(paste0("Prior cytotoxic backbone exposure: ", l), npct(d$BACKBONE_EXPO_R, l))
  for (l in levels(d$PRIOR_CLASS_CAT_R)) add(paste0("Timeline-derived prior drug classes: ", tolower(l)), npct(d$PRIOR_CLASS_CAT_R, l))
  add("Timeline-derived prior drug classes, median (IQR)", medi(d$PRIOR_CLASS_N_R))
  for (l in levels(d$OUTSIDE_PRIOR_R)) add(paste0("Systemic therapy before MSK referral: ", tolower(l)), npct(d$OUTSIDE_PRIOR_R, l))
  for (l in levels(d$COMBO_AT_INDEX_R)) add(paste0("Anti-EGFR ", tolower(l)), npct(d$COMBO_AT_INDEX_R, l))
  for (l in levels(d$ECOG_CAT_R)) add(paste0("ECOG performance status: ", l), npct(d$ECOG_CAT_R, l))
  add("Stage IV diagnosis to anti-EGFR, days, median (IQR)", medi(d$STAGE4_TO_ANTI_EGFR_DAYS_R))
  add("Specimen acquisition to anti-EGFR, days, median (IQR)", medi(d$SPECIMEN_TO_ANTI_EGFR_DAYS_R))
  do.call(rbind, out)
}
r_all <- desc_rows(x); r_pos <- desc_rows(x[BYPASS_HIGH_CONFIDENCE_R == 1]); r_neg <- desc_rows(x[BYPASS_HIGH_CONFIDENCE_R == 0])
tabR1 <- data.frame(Characteristic = r_all[, 1], `All (n=191)` = r_all[, 2],
                    `Bypass-positive (n=23)` = r_pos[, 2], `Bypass-negative (n=168)` = r_neg[, 2], check.names = FALSE)
fwrite(tabR1, file.path(OUT_DIR, "22_TableS11_treatment_history_stage_characteristics.csv"))

# ------------------------------------------------------------------------------
# 4. rwPFS and OS sensitivity models in the anti-EGFR cohort (Tables S12a, S12b)
# ------------------------------------------------------------------------------
base_cov <- "GENDER_F_R + SUBSITE_F_R + ANTI_EGFR_AGENT_F_R + LOG_STAGE4_TO_ANTI_EGFR_R"
fit_hr <- function(d, endpoint, rhs, label) {
  f <- as.formula(sprintf("Surv(%s_DAYS_LOCKED_R, %s_EVENT_LOCKED_R) ~ BYPASS_HIGH_CONFIDENCE_R + %s",
                          endpoint, endpoint, rhs))
  m <- coxph(f, data = d, ties = "efron")
  s <- summary(m)$coefficients["BYPASS_HIGH_CONFIDENCE_R", ]
  ci <- exp(confint(m)["BYPASS_HIGH_CONFIDENCE_R", ])
  mm <- model.matrix(m)
  zp <- tryCatch(cox.zph(m)$table["GLOBAL", "p"], error = function(e) NA_real_)
  data.table(analysis = label, endpoint = ifelse(endpoint == "RWPFS", "rwPFS", "OS"),
             n = m$n, events = m$nevent, n_bypass = as.integer(sum(mm[, "BYPASS_HIGH_CONFIDENCE_R"])),
             n_param = ncol(mm), hr = unname(exp(s["coef"])), lcl = unname(ci[1]), ucl = unname(ci[2]),
             p = unname(s["Pr(>|z|)"]), ph_global_p = zp, converged = m$iter < 30)
}
specs <- list(
  list("Primary locked model (reference)", x, paste(base_cov, "+ PRIOR_HEAVY_F_R")),
  list("Locked model + ECOG category (missing as category) + concurrent cytotoxic therapy", x, paste(base_cov, "+ PRIOR_HEAVY_F_R + ECOG_CAT_R + COMBO_AT_INDEX_R")),
  list("Prior cytotoxic backbone exposure (3 categories) instead of the prior-therapy covariate", x, paste(base_cov, "+ BACKBONE_EXPO_R")),
  list("Timeline-derived prior drug-class count (0/1/>=2) instead of the prior-therapy covariate", x, paste(base_cov, "+ PRIOR_CLASS_CAT_R")),
  list("Stratified by prior cytotoxic backbone exposure", x, paste(base_cov, "+ strata(BACKBONE_EXPO_R)")),
  list("Fully extended model", x, paste(base_cov, "+ PRIOR_CLASS_CAT_R + ECOG_CAT_R + COMBO_AT_INDEX_R + OUTSIDE_PRIOR_R")),
  list("Registry stage IV at diagnosis with radiology corroboration", x[STAGE4_IMG_CORROB_R == 1], paste(base_cov, "+ PRIOR_HEAVY_F_R")),
  list("Model C + ECOG 0 vs >=1 within 90 days, complete case (Table S5 model D)", x[CC90_R == 1], paste(modelC_cov, "+ ECOG90_F_R")),
  list("Same 90-day subset, model C without ECOG", x[CC90_R == 1], modelC_cov),
  list("Model C + ECOG 0 vs >=1 within 180 days, complete case", x[CC180_R == 1], paste(modelC_cov, "+ ECOG180_F_R")),
  list("Same 180-day subset, model C without ECOG", x[CC180_R == 1], modelC_cov),
  list("No systemic therapy before MSK referral", x[OUTSIDE_PRIOR_R == "No"], paste(base_cov, "+ PRIOR_HEAVY_F_R")),
  list("Not previously exposed to both oxaliplatin and irinotecan", x[BOTH_BACKBONES_R == 0], base_cov),
  list("Previously exposed to both oxaliplatin and irinotecan", x[BOTH_BACKBONES_R == 1], base_cov)
)
tabR2 <- rbindlist(lapply(c("RWPFS", "OS"), function(ep)
  rbindlist(lapply(specs, function(s) fit_hr(s[[2]], ep, s[[3]], s[[1]])))))
# reproduction checks against the locked manuscript results
stopifnot(abs(tabR2[analysis == "Primary locked model (reference)" & endpoint == "rwPFS"]$hr - 2.334) < 0.005,
          abs(tabR2[analysis == "Primary locked model (reference)" & endpoint == "OS"]$hr - 1.474) < 0.005,
          tabR2[analysis == "Primary locked model (reference)" & endpoint == "rwPFS"]$events == 183,
          tabR2[grepl("Table S5 model D", analysis) & endpoint == "rwPFS"]$n == 87)

int_p <- sapply(c("RWPFS", "OS"), function(ep) {
  f <- as.formula(sprintf("Surv(%s_DAYS_LOCKED_R, %s_EVENT_LOCKED_R) ~ BYPASS_HIGH_CONFIDENCE_R * BOTH_BACKBONES_R + %s",
                          ep, ep, base_cov))
  m <- coxph(f, data = x, ties = "efron")
  summary(m)$coefficients["BYPASS_HIGH_CONFIDENCE_R:BOTH_BACKBONES_R", "Pr(>|z|)"]
})
fwrite(tabR2, file.path(OUT_DIR, "22_TableS12_sensitivity_models.csv"))

# Proportional-hazards diagnostics (Table S12b): term-level cox.zph for every
# model whose global test gave P < 0.05, plus PH-addressed refits in which the
# offending adjustment covariate is moved into strata()
ph_terms <- rbindlist(lapply(c("RWPFS", "OS"), function(ep) rbindlist(lapply(specs, function(s) {
  f <- as.formula(sprintf("Surv(%s_DAYS_LOCKED_R, %s_EVENT_LOCKED_R) ~ BYPASS_HIGH_CONFIDENCE_R + %s", ep, ep, s[[3]]))
  z <- cox.zph(coxph(f, data = s[[2]], ties = "efron"))$table
  if (z["GLOBAL", "p"] >= 0.05) return(NULL)
  data.table(analysis = s[[1]], endpoint = ifelse(ep == "RWPFS", "rwPFS", "OS"), term = rownames(z),
             chisq = z[, "chisq"], df = z[, "df"], p = z[, "p"])
}))))
fwrite(ph_terms, file.path(OUT_DIR, "22_PH_term_diagnostics.csv"))
ph_fix_specs <- list(
  list("Model C + ECOG 0 vs >=1 within 180 days, complete case; stratified by panel generation", x[CC180_R == 1],
       paste(sub("PANEL_OLD_F_R", "strata(PANEL_OLD_F_R)", modelC_cov, fixed = TRUE), "+ ECOG180_F_R")),
  list("Same 180-day subset, model C without ECOG; stratified by panel generation", x[CC180_R == 1],
       sub("PANEL_OLD_F_R", "strata(PANEL_OLD_F_R)", modelC_cov, fixed = TRUE)))
ph_fix <- rbindlist(lapply(ph_fix_specs, function(s) fit_hr(s[[2]], "RWPFS", s[[3]], s[[1]])))
fwrite(ph_fix, file.path(OUT_DIR, "22_PH_addressed_refits.csv"))
ph_fix_terms <- rbindlist(lapply(ph_fix_specs, function(s) {
  f <- as.formula(sprintf("Surv(RWPFS_DAYS_LOCKED_R, RWPFS_EVENT_LOCKED_R) ~ BYPASS_HIGH_CONFIDENCE_R + %s", s[[3]]))
  z <- cox.zph(coxph(f, data = s[[2]], ties = "efron"))$table
  data.table(analysis = s[[1]], endpoint = "rwPFS", term = rownames(z), chisq = z[, "chisq"], df = z[, "df"], p = z[, "p"])
}))
fwrite(ph_fix_terms, file.path(OUT_DIR, "22_PH_term_diagnostics_stratified_refits.csv"))
# which propensity covariates were missing in the excluded comparator patients
ps_vars <- all.vars(readRDS(file.path(A7, "04_results", "B1_14_bevacizumab_comparator_feasibility.rds"))$ps_formula)[-1]
ps_missing <- fe_audit[PS_COMPLETE_R == FALSE, .(missing_covariates = paste(ps_vars[is.na(unlist(.SD))], collapse = "; ")),
                       by = .(PATIENT_ID, FIRST_BIOLOGIC_R, BYPASS_HIGH_CONFIDENCE_R), .SDcols = ps_vars]

# E-values for the primary anti-EGFR-cohort associations only (common-outcome conversion)
evalue_hr <- function(hr) {
  rr <- (1 - 0.5^sqrt(hr)) / (1 - 0.5^sqrt(1 / hr))
  rr <- ifelse(rr < 1, 1 / rr, rr)
  rr + sqrt(rr * (rr - 1))
}
prim <- tabR2[analysis == "Primary locked model (reference)"]
tabEV <- data.table(analysis = paste0("Primary adjusted ", prim$endpoint, " HR"), hr = prim$hr, lcl = prim$lcl,
                    evalue_point = evalue_hr(prim$hr), evalue_ci = ifelse(prim$lcl > 1, evalue_hr(prim$lcl), 1))
fwrite(tabEV, file.path(OUT_DIR, "22_TableS14_Evalues.csv"))

# ------------------------------------------------------------------------------
# 5. Comparator: residual confounding (Table S13)
#    M1 = frozen overlap weights (reference; reproduces Table S8)
#    M3 = overlap weights from an extended propensity model (frozen 10
#         covariates + ECOG category + stage IV at diagnosis + therapy before
#         referral + timeline-derived prior drug-class count), with a full
#         estimation-pipeline bootstrap (propensity model, weights and Cox
#         model refit in each resample)
#    M4 = patients not previously exposed to both oxaliplatin and irinotecan,
#         with propensity model and overlap weights re-estimated within the
#         subgroup
# ------------------------------------------------------------------------------
fit_int <- function(d, w) {
  m <- coxph(Surv(RWPFS_DAYS_B14B_R, RWPFS_EVENT_B14B_R) ~ ANTI_EGFR_TREAT_R * BYPASS_HIGH_CONFIDENCE_R,
             data = d, weights = w, cluster = PATIENT_ID, robust = TRUE, ties = "efron")
  term <- grep(":", names(coef(m)), value = TRUE)
  list(m = m, b = unname(coef(m)[term]), se = unname(sqrt(diag(vcov(m)))[term]))
}
drop_const <- function(fml, d) {
  # drop terms that are constant in d (e.g., a factor level absent in a subgroup or resample)
  tl <- attr(terms(fml), "term.labels")
  keep <- tl[vapply(tl, function(v) length(unique(na.omit(d[[v]]))) > 1, logical(1))]
  as.formula(paste(paste(deparse(fml[[2]]), collapse = ""), "~", paste(keep, collapse = " + ")), env = environment(fml))
}
ow <- function(fml, d) {
  fml <- drop_const(fml, d)
  f <- suppressWarnings(glm(fml, data = d, family = binomial()))
  p <- predict(f, type = "response")
  list(w = ifelse(d$ANTI_EGFR_TREAT_R == 1, 1 - p, p), fit = f)
}
ps_frozen <- readRDS(file.path(A7, "04_results", "B1_14_bevacizumab_comparator_feasibility.rds"))$ps_formula
ps_ext <- update(ps_frozen, . ~ . + ECOG_CAT_R + STAGE4_AT_DX_F_R + OUTSIDE_PRIOR_R + PRIOR_CLASS_CAT_R)
ps_ext <- drop_const(ps_ext, dw)                            # stage IV at diagnosis is constant if all are stage IV
cat("Extended propensity model:", deparse(ps_ext, width.cutoff = 500), "\n")
stopifnot(nrow(model.frame(ps_ext, dw)) == nrow(dw))     # no rows dropped by the propensity model
o3 <- ow(ps_ext, dw); dw[, OW_EXT_R := o3$w]
early <- dw[PRIOR_MAJOR_N_R < 2]
o4 <- ow(ps_frozen, early); early[, OW_SUB_R := o4$w]

r1 <- fit_int(dw, dw$OVERLAP_WEIGHT_R); r3 <- fit_int(dw, dw$OW_EXT_R); r4 <- fit_int(early, early$OW_SUB_R)
stopifnot(abs(exp(r1$b) - 3.083) < 0.005)

set.seed(20261124)
B <- 2000
boot_mat <- t(vapply(seq_len(B), function(i) {
  bi <- dw[sample.int(nrow(dw), replace = TRUE)]
  bi[, PATIENT_ID := paste0(PATIENT_ID, "_", seq_len(.N))]
  tryCatch({
    o <- ow(ps_ext, bi); r <- fit_int(bi, o$w); pr <- fitted(o$fit)
    c(b = r$b, se = r$se, glm_conv = as.numeric(o$fit$converged), p_min = min(pr), p_max = max(pr),
      cox_iter = r$m$iter[1])
  }, error = function(e) c(b = NA_real_, se = NA_real_, glm_conv = NA_real_, p_min = NA_real_, p_max = NA_real_, cox_iter = NA_real_))
}, numeric(6)))
boot <- boot_mat[, "b"]
boot_ok <- boot[is.finite(boot)]                         # prespecified: all resamples with a finite estimate
boot_ci <- exp(quantile(boot_ok, c(0.025, 0.975)))
# fit-quality diagnostics: propensity model converged without near-separation (fitted PS within 1e-6 of 0 or 1),
# Cox model converged within the default 20 iterations with a finite robust SE
boot_clean <- is.finite(boot_mat[, "b"]) & is.finite(boot_mat[, "se"]) & boot_mat[, "glm_conv"] == 1 &
  boot_mat[, "p_min"] > 1e-6 & boot_mat[, "p_max"] < 1 - 1e-6 & boot_mat[, "cox_iter"] < 20
boot_clean[is.na(boot_clean)] <- FALSE
boot_ci_clean <- exp(quantile(boot[boot_clean], c(0.025, 0.975)))
boot_diag <- data.table(resamples = B, finite_estimate = length(boot_ok),
  ps_not_converged = sum(boot_mat[, "glm_conv"] == 0, na.rm = TRUE),
  ps_near_separation = sum(boot_mat[, "p_min"] <= 1e-6 | boot_mat[, "p_max"] >= 1 - 1e-6, na.rm = TRUE),
  cox_not_converged = sum(boot_mat[, "cox_iter"] >= 20 | !is.finite(boot_mat[, "se"]), na.rm = TRUE),
  all_checks_passed = sum(boot_clean),
  ci_all_finite = sprintf("%.2f-%.2f", boot_ci[1], boot_ci[2]),
  ci_all_checks_passed = sprintf("%.2f-%.2f", boot_ci_clean[1], boot_ci_clean[2]))
fwrite(boot_diag, file.path(OUT_DIR, "22_bootstrap_fit_diagnostics.csv"))

smd_w <- function(v, t, w) {
  v <- as.numeric(v); m1 <- weighted.mean(v[t == 1], w[t == 1]); m0 <- weighted.mean(v[t == 0], w[t == 0])
  s <- sqrt((var(v[t == 1]) + var(v[t == 0])) / 2)
  if (!is.finite(s) || s == 0) {                      # zero within-group variance in both groups
    smd_zero_sd_n <<- smd_zero_sd_n + 1
    return(if (isTRUE(all.equal(m1, m0))) 0 else sign(m1 - m0) * Inf)  # differing constant means = complete imbalance
  }
  (m1 - m0) / s
}
smd_zero_sd_n <- 0
ess <- function(w) sum(w)^2 / sum(w^2)
mm_ext <- model.matrix(ps_ext, data = dw)[, -1, drop = FALSE]
bal <- data.table(variable = colnames(mm_ext),
                  smd_unweighted = apply(mm_ext, 2, smd_w, t = dw$ANTI_EGFR_TREAT_R, w = rep(1, nrow(dw))),
                  smd_M1 = apply(mm_ext, 2, smd_w, t = dw$ANTI_EGFR_TREAT_R, w = dw$OVERLAP_WEIGHT_R),
                  smd_M3 = apply(mm_ext, 2, smd_w, t = dw$ANTI_EGFR_TREAT_R, w = dw$OW_EXT_R))
mm_sub <- model.matrix(drop_const(ps_frozen, early), data = early)[, -1, drop = FALSE]
mm_sub <- mm_sub[, apply(mm_sub, 2, var) > 0, drop = FALSE]
max_smd_M4 <- max(abs(apply(mm_sub, 2, smd_w, t = early$ANTI_EGFR_TREAT_R, w = early$OW_SUB_R)))
fwrite(bal, file.path(OUT_DIR, "22_comparator_balance_extended_ps.csv"))

# four treatment x bypass cells: counts, events, ESS, within-stratum balance
cells <- rbindlist(lapply(list(list("M1", dw, dw$OVERLAP_WEIGHT_R, mm_ext), list("M3", dw, dw$OW_EXT_R, mm_ext)), function(s) {
  d <- s[[2]]; w <- s[[3]]; mm <- s[[4]]
  rbindlist(lapply(0:1, function(b) {
    i <- d$BYPASS_HIGH_CONFIDENCE_R == b
    keep <- apply(mm[i, , drop = FALSE], 2, var) > 0
    data.table(weights = s[[1]], bypass = ifelse(b == 1, "Bypass-positive", "Bypass-negative"),
               n_antiEGFR = sum(i & d$ANTI_EGFR_TREAT_R == 1), n_bev = sum(i & d$ANTI_EGFR_TREAT_R == 0),
               events_antiEGFR = sum(d$RWPFS_EVENT_B14B_R[i & d$ANTI_EGFR_TREAT_R == 1]),
               events_bev = sum(d$RWPFS_EVENT_B14B_R[i & d$ANTI_EGFR_TREAT_R == 0]),
               ess_antiEGFR = ess(w[i & d$ANTI_EGFR_TREAT_R == 1]), ess_bev = ess(w[i & d$ANTI_EGFR_TREAT_R == 0]),
               max_abs_smd_within_stratum = max(abs(apply(mm[i, keep, drop = FALSE], 2, smd_w,
                                                          t = d$ANTI_EGFR_TREAT_R[i], w = w[i]))))
  }))
}))
fwrite(cells, file.path(OUT_DIR, "22_TableS13b_comparator_cells.csv"))
# per-covariate balance inside the bypass-positive stratum (Table S13c)
ib <- dw$BYPASS_HIGH_CONFIDENCE_R == 1
smd_pos <- data.table(variable = colnames(mm_ext),
  smd_unweighted = apply(mm_ext[ib, ], 2, smd_w, t = dw$ANTI_EGFR_TREAT_R[ib], w = rep(1, sum(ib))),
  smd_M1 = apply(mm_ext[ib, ], 2, smd_w, t = dw$ANTI_EGFR_TREAT_R[ib], w = dw$OVERLAP_WEIGHT_R[ib]),
  smd_M3 = apply(mm_ext[ib, ], 2, smd_w, t = dw$ANTI_EGFR_TREAT_R[ib], w = dw$OW_EXT_R[ib]))
smd_pos <- smd_pos[order(-abs(smd_M3))]
stopifnot(abs(max(abs(smd_pos$smd_M3)) - cells[weights == "M3" & bypass == "Bypass-positive"]$max_abs_smd_within_stratum) < 1e-8)
fwrite(smd_pos, file.path(OUT_DIR, "22_TableS13c_bypass_positive_balance.csv"))

row3 <- function(lab, r, d, w, maxsmd, ci = NULL) {
  data.table(model = lab, n = r$m$n, events = r$m$nevent, n_bypass_pos = sum(d$BYPASS_HIGH_CONFIDENCE_R),
             interaction_hr = exp(r$b), lcl = exp(r$b - 1.96 * r$se), ucl = exp(r$b + 1.96 * r$se),
             p = 2 * pnorm(-abs(r$b / r$se)),
             boot_lcl = if (is.null(ci)) NA_real_ else ci[1], boot_ucl = if (is.null(ci)) NA_real_ else ci[2],
             ess_antiEGFR = ess(w[d$ANTI_EGFR_TREAT_R == 1]), ess_bev = ess(w[d$ANTI_EGFR_TREAT_R == 0]),
             ess_bypass_pos = ess(w[d$BYPASS_HIGH_CONFIDENCE_R == 1]), max_abs_smd = maxsmd)
}
tabR3 <- rbind(
  row3("M1: frozen overlap weights (reference)", r1, dw, dw$OVERLAP_WEIGHT_R, max(abs(bal$smd_M1))),
  row3("M3: overlap weights from extended propensity model", r3, dw, dw$OW_EXT_R, max(abs(bal$smd_M3)), boot_ci),
  row3("M4: not previously exposed to both oxaliplatin and irinotecan; weights re-estimated", r4, early, early$OW_SUB_R, max_smd_M4))
tabR3[, boot_valid := c(NA, length(boot_ok), NA)]
fwrite(tabR3, file.path(OUT_DIR, "22_TableS13a_comparator_residual_confounding.csv"))
fwrite(x[, .(PATIENT_ID, STAGE4_AT_DX_R, STAGE4_SOURCE_R, IMAGING_CLASS_R, BACKBONE_EXPO_R, PRIOR_CLASS_N_R,
             PRIOR_CLASS_CAT_R, OUTSIDE_PRIOR_R, ECOG_CAT_R, ECOG90_R, ECOG180_R, COMBO_AT_INDEX_R)],
       file.path(OUT_DIR, "22_derived_variables_anti_EGFR_cohort.csv"))

# ------------------------------------------------------------------------------
# 6. Supplementary Figure S6: forest plot of rwPFS sensitivity analyses
# ------------------------------------------------------------------------------
fp <- tabR2[endpoint == "rwPFS"]
short <- c("Primary locked model (reference)",
           "+ ECOG (missing as category) + concurrent cytotoxic",
           "Prior backbone exposure, 3 categories",
           "Prior drug-class count, 0/1/>=2",
           "Stratified by prior backbone exposure",
           "Fully extended model",
           "Radiology-corroborated stage IV at diagnosis",
           "Model C + ECOG <=90 d (Table S5 model D)",
           "Same 90-day subset, without ECOG",
           "Model C + ECOG <=180 d",
           "Same 180-day subset, without ECOG",
           "No systemic therapy before MSK referral",
           "Not exposed to both oxaliplatin and irinotecan",
           "Exposed to both oxaliplatin and irinotecan")
stopifnot(length(short) == nrow(fp))
fp[, label := sprintf("%s (n=%d; %d bypass+)", short, n, n_bypass)]
fp[, label := factor(label, levels = rev(label))]
fp[, txt := sprintf("%.2f (%.2f-%.2f)", hr, lcl, ucl)]
g <- ggplot(fp, aes(x = hr, y = label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), width = 0.2, orientation = "y") +
  geom_point(shape = 15, size = 2.4) +
  geom_text(aes(x = max(fp$ucl) * 1.6, label = txt), hjust = 0, size = 2.8) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
  coord_cartesian(xlim = c(min(0.5, min(fp$lcl)), max(fp$ucl) * 4), clip = "off") +
  labs(x = "rwPFS HR, bypass-positive vs negative (log scale)", y = NULL) +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.margin = margin(5, 40, 5, 5))
ggsave(file.path(FIG_DIR, "FigureS6_revision_sensitivity_forest.tiff"), g,
       width = 174, height = 110, units = "mm", dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(FIG_DIR, "FigureS6_revision_sensitivity_forest.png"), g,
       width = 174, height = 110, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, "FigureS6_revision_sensitivity_forest.pdf"), g, width = 174, height = 110, units = "mm")

# ------------------------------------------------------------------------------
# 7. Word three-line tables (Tables S11-S14)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(flextable); library(officer) })
ci3 <- function(e, l, h) sprintf("%.2f (%.2f-%.2f)", e, l, h)
pv <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
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
t2 <- tabR2[, .(Analysis = analysis, Endpoint = endpoint, `N (events)` = sprintf("%d (%d)", n, events),
                `Bypass-positive, n` = n_bypass, `HR (95% CI)` = ci3(hr, lcl, ucl), P = pv(p),
                `PH global P` = ifelse(is.na(ph_global_p), "-", pv(ph_global_p)))]
t3 <- tabR3[, .(Model = model, `N (events)` = sprintf("%d (%d)", n, events), `Bypass-positive, n` = n_bypass_pos,
                `Interaction HR (95% CI)` = ci3(interaction_hr, lcl, ucl), P = pv(p),
                `Bootstrap 95% CI` = ifelse(is.na(boot_lcl), "-", sprintf("%.2f-%.2f", boot_lcl, boot_ucl)),
                `ESS anti-EGFR / bevacizumab / bypass-positive` = sprintf("%.1f / %.1f / %.1f", ess_antiEGFR, ess_bev, ess_bypass_pos),
                `Max |SMD|` = ifelse(max_abs_smd < 0.001, "<0.001", sprintf("%.3f", max_abs_smd)))]
t3b <- cells[, .(Weights = weights, Stratum = bypass, `Anti-EGFR, n (events)` = sprintf("%d (%d)", n_antiEGFR, events_antiEGFR),
                 `Bevacizumab, n (events)` = sprintf("%d (%d)", n_bev, events_bev),
                 `ESS anti-EGFR / bevacizumab` = sprintf("%.1f / %.1f", ess_antiEGFR, ess_bev),
                 `Max |SMD| within stratum` = sprintf("%.3f", max_abs_smd_within_stratum))]
lab_term <- c(BYPASS_HIGH_CONFIDENCE_R = "EGFR-bypass status", GLOBAL = "Global test", PANEL_OLD_F_R = "Panel generation",
              MET_SITE_COUNT_NUM_R = "Recent metastatic-site count", SAMPLE_METASTATIC_F_R = "Specimen type (metastatic vs primary)")
t2b <- rbind(
  ph_terms[term %in% c("BYPASS_HIGH_CONFIDENCE_R", "GLOBAL") | p < 0.05,
           .(Analysis = analysis, Endpoint = endpoint, Term = fifelse(term %in% names(lab_term), lab_term[term], term), `Chi-square (df)` = sprintf("%.2f (%d)", chisq, as.integer(df)), P = pv(p))],
  ph_fix[, .(Analysis = analysis, Endpoint = endpoint, Term = "PH-addressed refit: bypass HR (95% CI)",
             `Chi-square (df)` = ci3(hr, lcl, ucl), P = sprintf("%s (PH global P %s)", pv(p), pv(ph_global_p)))],
  ph_fix_terms[term %in% c("BYPASS_HIGH_CONFIDENCE_R", "GLOBAL") | p < 0.05,
               .(Analysis = analysis, Endpoint = endpoint, Term = paste0("Stratified refit: ", fifelse(term %in% names(lab_term), lab_term[term], term)),
                 `Chi-square (df)` = sprintf("%.2f (%d)", chisq, as.integer(df)), P = pv(p))])
lab_var <- c(GENDER_F_RMale = "Male sex", `SUBSITE_F_Rrectum/rectosigmoid` = "Rectum/rectosigmoid primary",
  LOG_STAGE4_TO_BIOLOGIC_R = "Log interval, stage IV diagnosis to biologic", PRIOR_MAJOR_N_F_R1 = "Prior major cytotoxic classes: 1",
  PRIOR_MAJOR_N_F_R2 = "Prior major cytotoxic classes: 2", ACTIVE_IRI_F_RYes = "Active irinotecan context",
  ACTIVE_OX_F_RYes = "Active oxaliplatin context", `PANEL_OLD_F_RIMPACT341/410` = "Older panel (IMPACT341/410)",
  SAMPLE_METASTATIC_F_RMetastatic = "Metastatic specimen", LOG_SPECIMEN_TO_BIOLOGIC_R = "Log interval, specimen to biologic",
  MET_SITE_COUNT_90D_R = "Recent metastatic-site count", ECOG_CAT_R1 = "ECOG 1", `ECOG_CAT_R2 or higher` = "ECOG 2 or higher",
  `ECOG_CAT_RNot recorded within 180 days` = "ECOG not recorded within 180 days", OUTSIDE_PRIOR_RYes = "Systemic therapy before MSK referral: yes",
  OUTSIDE_PRIOR_RUnknown = "Systemic therapy before MSK referral: unknown", `PRIOR_CLASS_CAT_ROne prior class` = "Prior drug classes: one",
  `PRIOR_CLASS_CAT_RTwo or more prior classes` = "Prior drug classes: two or more")
stopifnot(all(smd_pos$variable %in% names(lab_var)))
t3c <- smd_pos[, .(Covariate = unname(lab_var[variable]), Unweighted = sprintf("%.3f", smd_unweighted),
                   `M1 weights` = sprintf("%.3f", smd_M1), `M3 weights` = sprintf("%.3f", smd_M3))]
t4 <- tabEV[, .(Estimate = analysis, `HR (lower 95% CL)` = sprintf("%.2f (%.2f)", hr, lcl),
                `E-value, point estimate` = sprintf("%.2f", evalue_point),
                `E-value, confidence limit` = sprintf("%.2f", evalue_ci))]
doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "landscape")))
add_tab <- function(doc, cap, df, note) {
  doc <- body_add_par(doc, cap, style = "Normal")
  doc <- body_add_flextable(doc, three_line(df))
  doc <- body_add_par(doc, note, style = "Normal")
  body_add_par(doc, "", style = "Normal")
}
doc <- add_tab(doc, "Supplementary Table S11. Treatment-history proxies, stage at initial diagnosis, and related characteristics of the anti-EGFR cohort (post hoc)",
  tabR1, paste("Values are n (%) unless stated otherwise. Stage at initial diagnosis was taken from the first tumour-registry primary colorectal",
               "record (SEER summary stage, or AJCC clinical stage group). Radiology corroboration required a distant organ site (liver, lung,",
               "adrenal, bone, central nervous system, pleura, or reproductive organs) in the radiology-report-derived tumour-site timeline",
               "within 120 days before or after diagnosis and before the index biologic; intra-abdominal, lymph-node, and other sites were not",
               "counted, and the timeline contains only reports with an extracted tumour site, so absence of a record does not mean absence of",
               "imaging. Prior cytotoxic backbone exposure and the timeline-derived prior drug-class count are treatment-history proxies, not",
               "chart-verified lines of therapy; the drug-class count includes drugs started from 30 days before stage IV diagnosis to 28 days",
               "before the index biologic and counts drug classes rather than regimens. ECOG, Eastern Cooperative Oncology Group; IQR,",
               "interquartile range; MSK, Memorial Sloan Kettering."))
doc <- add_tab(doc, "Supplementary Table S12a. Sensitivity analyses addressing treatment history, stage at diagnosis, performance status, and baseline confounding (post hoc)",
  t2, paste0("All models compare bypass-positive with bypass-negative patients. Unless stated otherwise, models include sex, primary subsite,",
             " anti-EGFR agent, and log-transformed interval from stage IV diagnosis to anti-EGFR initiation, plus the covariates named in each row.",
             " Model C additionally includes prior exposure to 0/1/2 major cytotoxic classes, active irinotecan and oxaliplatin context, panel",
             " generation, specimen type, log specimen-to-index interval, and recent metastatic-site count; ECOG models use the nearest pre-index",
             " record within the stated window (0 vs >=1; complete case). Subgroup models by prior backbone exposure use the four core covariates.",
             " Interaction P for bypass x prior exposure to both oxaliplatin and irinotecan: rwPFS ", pv(int_p["RWPFS"]), "; OS ", pv(int_p["OS"]), ".",
             " PH global P, Grambsch-Therneau global test. CI, confidence interval; HR, hazard ratio; OS, overall survival; PH, proportional",
             " hazards; rwPFS, real-world progression-free survival."))
doc <- add_tab(doc, "Supplementary Table S12b. Proportional-hazards diagnostics for sensitivity models with global P < 0.05 (post hoc)",
  t2b, paste("Term-level Grambsch-Therneau tests (cox.zph) for every model in Table S12a whose global test gave P < 0.05 (bypass term, global test, and",
             "any term with P < 0.05 shown; all terms in 22_PH_term_diagnostics.csv), followed by",
             "refits in which the nonproportional adjustment covariate (panel generation) was moved into strata(), with their own term-level tests.",
             "In the 90-day OS subset the",
             "bypass term itself was nonproportional, so its HR is a time-averaged estimate and is not interpreted."))
doc <- add_tab(doc, "Supplementary Table S13a. Exploratory comparator analysis: residual confounding and treatment-history restriction (post hoc)",
  t3, paste("The comparator comprised the 319 propensity-score-complete patients; 3 anti-EGFR and 6 bevacizumab patients lacked a propensity covariate",
            "(recent metastatic-site count in 8; interval from stage IV diagnosis to biologic in 1).",
            "Interaction HR is the ratio of anti-EGFR-versus-bevacizumab rwPFS HRs in bypass-positive versus bypass-negative patients;",
            "confidence intervals use robust variance. M3 added ECOG category, systemic therapy before MSK referral, and the",
            "timeline-derived prior drug-class count to the 10 frozen propensity covariates (stage IV at initial diagnosis was constant and omitted); its bootstrap interval refits",
            "the propensity model, weights, and Cox model in each of 2,000 patient-level resamples (percentile method); all propensity and Cox",
            sprintf("models converged, and excluding the %d resamples with near-separation (fitted propensity score within 1e-6 of 0 or 1) gave %s. In M4, the frozen",
                    boot_diag$ps_near_separation, sub("-", " to ", boot_diag$ci_all_checks_passed)),
            "propensity model and overlap weights were re-estimated within the subgroup. Max |SMD| is computed over the covariates of the",
            "respective propensity model (for M1 and M3, the extended covariate set); overlap weights give exact mean balance on the covariates of",
            "the model from which they were estimated. M4 is reported for completeness; with 14 bypass-positive patients it is unstable and not interpreted. ESS, effective sample size; SMD, standardized mean difference."))
doc <- add_tab(doc, "Supplementary Table S13b. Comparator treatment-by-bypass cells: sample size, events, effective sample size, and within-stratum balance (post hoc)",
  t3b, paste("Within-stratum balance compares anti-EGFR with bevacizumab patients inside each bypass stratum over the extended covariate set;",
             "overlap weights balance covariates in the overall population, not necessarily within each stratum."))
doc <- add_tab(doc, "Supplementary Table S13c. Covariate balance within the bypass-positive comparator stratum (post hoc)",
  t3c, paste("Standardized mean differences (anti-EGFR minus bevacizumab) among the 37 bypass-positive comparator patients (23 anti-EGFR, 14",
             "bevacizumab), using the pooled unweighted standard deviation, before weighting and under the frozen (M1) and extended (M3) overlap",
             "weights; sorted by absolute M3 value. Propensity models did not include bypass status, so balance within this stratum is not targeted."))
doc <- add_tab(doc, "Supplementary Table S14. E-values for unmeasured confounding of the primary associations (post hoc)",
  t4, paste("The E-value is the minimum strength of association, on the risk-ratio scale, that an unmeasured confounder would need to have",
            "with both exposure and outcome to fully explain away the observed estimate (point) or to move the confidence limit to the null;",
            "HRs were converted with the common-outcome approximation. E-values do not exclude unmeasured confounding."))
print(doc, target = file.path(OUT_DIR, "22_revision_tables_three_line.docx"))

cat("\n==== Key results ====\n")
print(tabR2[, .(analysis = substr(analysis, 1, 60), endpoint, n, events, n_bypass, hr = round(hr, 2), lcl = round(lcl, 2),
                ucl = round(ucl, 2), p = signif(p, 2), ph = round(ph_global_p, 3))], nrows = 100)
cat("Interaction P (bypass x both backbones): rwPFS", signif(int_p["RWPFS"], 3), "; OS", signif(int_p["OS"], 3), "\n")
print(tabR1)
print(tabR3); print(cells); print(tabEV)
cat("\nPH term diagnostics (models with global P < 0.05):\n"); print(ph_terms[term %in% c("BYPASS_HIGH_CONFIDENCE_R", "GLOBAL") | p < 0.05, .(analysis = substr(analysis, 1, 45), endpoint, term, chisq = round(chisq, 2), p = signif(p, 2))], nrows = 100)
cat("\nPH-addressed refits:\n"); print(ph_fix[, .(analysis, n, events, hr = round(hr, 2), lcl = round(lcl, 2), ucl = round(ucl, 2), p = signif(p, 2), ph = round(ph_global_p, 3))])
cat("\nPH terms in stratified refits (bypass, global, P < 0.05):\n"); print(ph_fix_terms[term %in% c("BYPASS_HIGH_CONFIDENCE_R", "GLOBAL") | p < 0.05, .(analysis = substr(analysis, 1, 45), term, chisq = round(chisq, 2), p = round(p, 3))])
cat("\nBalance within bypass-positive comparator stratum:\n"); print(smd_pos, digits = 3)
cat("\nPropensity-score-incomplete patients excluded from comparator:\n"); print(ps_incomplete); print(ps_missing)
cat("Bootstrap resamples valid:", length(boot_ok), "of", B, "\n")
cat("\nBootstrap fit diagnostics:\n"); print(boot_diag)
cat("SMD computations with zero pooled SD:", smd_zero_sd_n, "\n")
cat("\nOutputs written to:", OUT_DIR, "\nFigure S6 written to:", FIG_DIR, "\n")
sessionInfo()
