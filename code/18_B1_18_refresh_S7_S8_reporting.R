# ============================================================
# B1-18 — Reporting-only refresh of Supplementary Tables S7/S8
# Project: A7 B-line anti-EGFR / exploratory bevacizumab comparator
#
# PURPOSE
#   1) Keep the original frozen S7 A-D and S8 A-E content.
#   2) Add the outcome-blind component-composition panel to S7.
#   3) Add B1-16 shared <=730-day comparator, contrasts, LOO, and PH
#      diagnostics to S8 so every number reported in the manuscript has
#      a displayed supplementary source.
#
# IMPORTANT
#   - REPORTING ONLY. No exposure redefinition.
#   - No new outcome model is fitted here.
#   - No figure is regenerated.
#   - Existing TableS7/TableS8 files are overwritten, per project rule.
#
# INPUTS
#   04_results/B1_16_specimen_timing_comparator_stress_test.rds
#   07_tables/TableS7_bevacizumab_comparator_feasibility.docx
#   07_tables/TableS8_exploratory_comparator_rwPFS.docx
#
# OUTPUTS (OVERWRITE)
#   07_tables/TableS7_bevacizumab_comparator_feasibility.docx
#   07_tables/TableS8_exploratory_comparator_rwPFS.docx
#   06_logs_and_audit/B1_18_refresh_S7_S8_reporting.txt
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory
results_dir  <- file.path(project_root, "04_results")
table_dir    <- file.path(project_root, "07_tables")
audit_dir    <- file.path(project_root, "06_logs_and_audit")

s7_file <- file.path(table_dir, "TableS7_bevacizumab_comparator_feasibility.docx")
s8_file <- file.path(table_dir, "TableS8_exploratory_comparator_rwPFS.docx")
b116_file <- file.path(results_dir, "B1_16_specimen_timing_comparator_stress_test.rds")
audit_file <- file.path(audit_dir, "B1_18_refresh_S7_S8_reporting.txt")

for (f in c(s7_file, s8_file, b116_file)) {
    if (!file.exists(f)) stop("Missing required input: ", f)
}

pkgs <- c("officer", "flextable", "data.table")
for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(officer)
library(flextable)
library(data.table)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
fmt_p <- function(p) {
    ifelse(
        is.na(p),
        NA_character_,
        ifelse(p < 0.001, "<.001", sprintf("%.3f", p))
    )
}

fmt_hr_ci <- function(hr, lo, hi) {
    sprintf("%.2f (%.2f–%.2f)", hr, lo, hi)
}

three_line_ft <- function(df, font_size = 8) {

    # Word-extracted tables may contain blank or duplicated visible headers.
    # flextable uses column names as internal col_keys, so use guaranteed-safe
    # internal keys and restore the original labels only for display.
    df <- as.data.frame(
        df,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (ncol(df) < 1L) {
        stop("three_line_ft(): input table has zero columns.")
    }

    display_labels <- names(df)

    if (is.null(display_labels)) {
        display_labels <- rep("", ncol(df))
    }

    display_labels[is.na(display_labels)] <- ""

    internal_keys <- paste0(
        ".B118_COL_",
        seq_len(ncol(df))
    )

    names(df) <- internal_keys

    ft <- flextable(
        df,
        col_keys = internal_keys
    )

    ft <- set_header_labels(
        ft,
        values = setNames(
            as.list(display_labels),
            internal_keys
        )
    )

    ft <- border_remove(ft)

    bd_top <- fp_border(
        color = "black",
        width = 1
    )

    bd_mid <- fp_border(
        color = "black",
        width = 0.6
    )

    ft <- hline_top(
        ft,
        border = bd_top,
        part = "header"
    )

    ft <- hline_bottom(
        ft,
        border = bd_mid,
        part = "header"
    )

    ft <- hline_bottom(
        ft,
        border = bd_top,
        part = "body"
    )

    ft <- bold(
        ft,
        part = "header"
    )

    ft <- fontsize(
        ft,
        size = font_size,
        part = "all"
    )

    ft <- font(
        ft,
        fontname = "Arial",
        part = "all"
    )

    ft <- align(
        ft,
        align = "left",
        part = "all"
    )

    ft <- valign(
        ft,
        valign = "center",
        part = "all"
    )

    ft <- padding(
        ft,
        padding.top = 1,
        padding.bottom = 1,
        padding.left = 2,
        padding.right = 2,
        part = "all"
    )

    ft <- autofit(ft)

    ft
}

# Read all existing Word tables into plain character data.frames.
# This lets us preserve the already-frozen S7 A-D / S8 A-E values while
# rebuilding the documents cleanly in the correct panel order.
extract_docx_tables <- function(path) {

    x <- read_docx(path)
    sm <- docx_summary(x)

    z <- sm[
        sm$content_type == "table cell",
        ,
        drop = FALSE
    ]

    if (nrow(z) == 0L) {
        stop("No tables found in: ", path)
    }

    z$text[is.na(z$text)] <- ""

    idx <- sort(
        unique(z$table_index)
    )

    out <- lapply(
        idx,
        function(k) {

            a <- z[
                z$table_index == k,
                ,
                drop = FALSE
            ]

            nr <- max(
                a$row_id,
                na.rm = TRUE
            )

            nc <- max(
                a$cell_id,
                na.rm = TRUE
            )

            mat <- matrix(
                "",
                nrow = nr,
                ncol = nc
            )

            # A Word table cell may contain multiple paragraphs. Preserve them
            # all rather than overwriting earlier text from the same cell.
            for (r in seq_len(nr)) {
                for (cc in seq_len(nc)) {

                    vals <- a$text[
                        a$row_id == r &
                        a$cell_id == cc
                    ]

                    vals <- trimws(
                        as.character(vals)
                    )

                    vals <- vals[
                        nzchar(vals)
                    ]

                    if (length(vals) > 0L) {
                        mat[r, cc] <- paste(
                            vals,
                            collapse = " "
                        )
                    }
                }
            }

            hdr <- as.character(
                mat[1, ]
            )

            hdr[is.na(hdr)] <- ""

            d <- as.data.frame(
                mat[-1, , drop = FALSE],
                stringsAsFactors = FALSE,
                check.names = FALSE
            )

            names(d) <- hdr

            d
        }
    )

    out
}

add_heading <- function(doc, txt, level = 1) {
    body_add_par(doc, txt, style = if (level == 1) "heading 1" else "heading 2")
}

add_note <- function(doc, txt) {
    body_add_par(doc, txt, style = "Normal")
}

landscape_section <- prop_section(
    page_size = page_size(orient = "landscape"),
    page_margins = page_mar(top = 0.45, bottom = 0.45, left = 0.45, right = 0.45)
)

# ------------------------------------------------------------
# Frozen B1-16 source
# ------------------------------------------------------------
b116 <- readRDS(b116_file)

if (is.null(b116$restricted_730$interactions) ||
    is.null(b116$restricted_730$contrasts) ||
    is.null(b116$bev_bypass_pos_loo$results) ||
    is.null(b116$outcome_blind_component_composition)) {
    stop("B1-16 RDS is missing required slots.")
}

int730 <- as.data.table(b116$restricted_730$interactions)
con730 <- as.data.table(b116$restricted_730$contrasts)
loo <- as.data.table(b116$bev_bypass_pos_loo$results)
arm730 <- as.data.table(b116$restricted_730$arm_summary)
comp <- as.data.table(b116$outcome_blind_component_composition)
ph730 <- as.data.table(b116$restricted_730$ph_diagnostics)

# Frozen hard anchors.
stopifnot(nrow(int730) == 2L)
m1 <- int730[MODEL_ID == "M1"]
m2 <- int730[MODEL_ID == "M2"]
if (nrow(m1) != 1L || nrow(m2) != 1L) stop("Could not identify B1-16 M1/M2.")
if (abs(m1$INTERACTION_HR - 3.195494) > 0.001) stop("B1-16 M1 anchor mismatch.")
if (abs(m2$INTERACTION_HR - 3.210642) > 0.001) stop("B1-16 M2 anchor mismatch.")
if (nrow(loo) != 14L || any(loo$STATUS != "SUCCESS")) stop("B1-16 LOO anchor mismatch.")
if (any(loo$INTERACTION_HR <= 1) || any(loo$LCL95 <= 1)) stop("Unexpected LOO direction/CI anchor.")

# ------------------------------------------------------------
# Prepare S7 panel E
# ------------------------------------------------------------
comp[, arm := fifelse(Treatment == "anti-EGFR", "Anti-EGFR", "Bevacizumab")]
anti_comp <- comp[arm == "Anti-EGFR", .(Component, anti_n = N, anti_pct = `Percent of HC+ carriers`)]
bev_comp  <- comp[arm == "Bevacizumab", .(Component, bev_n = N, bev_pct = `Percent of HC+ carriers`)]
s7e <- merge(anti_comp, bev_comp, by = "Component", all = TRUE, sort = FALSE)
s7e[, `Anti-EGFR, n` := anti_n]
s7e[, `% of HC+ (n=23)` := sprintf("%.1f%%", anti_pct)]
s7e[, `Bevacizumab, n` := bev_n]
s7e[, `% of HC+ (n=16)` := sprintf("%.1f%%", bev_pct)]
s7e <- s7e[, .(Component, `Anti-EGFR, n`, `% of HC+ (n=23)`, `Bevacizumab, n`, `% of HC+ (n=16)`)]

# ------------------------------------------------------------
# Prepare S8 panels F-I
# ------------------------------------------------------------
# F1: restricted arm summary
s8f1 <- copy(arm730)
setnames(s8f1, old = names(s8f1), new = names(s8f1))
# Handle column names as saved by B1-16.
tr_col <- grep("Treatment", names(s8f1), value = TRUE)[1]
n_col  <- grep("^N$", names(s8f1), value = TRUE)[1]
ev_col <- grep("Events", names(s8f1), value = TRUE, ignore.case = TRUE)[1]
bp_col <- grep("Bypass", names(s8f1), value = TRUE, ignore.case = TRUE)[1]
sp_col <- grep("Specimen age", names(s8f1), value = TRUE, ignore.case = TRUE)[1]
s8f1 <- s8f1[, .(
    Treatment = get(tr_col),
    N = get(n_col),
    Events = get(ev_col),
    `Bypass+, n` = get(bp_col),
    `Specimen age, days, median [IQR]` = get(sp_col)
)]

# F2: restricted interaction models
s8f2 <- int730[, .(
    Model = MODEL,
    N,
    Events = EVENTS,
    `Interaction HR (95% CI)` = fmt_hr_ci(INTERACTION_HR, LCL95, UCL95),
    `P value` = fmt_p(P)
)]

# G: restricted contrasts
s8g <- con730[, .(
    Model = MODEL_ID,
    Contrast,
    `HR (95% CI)` = fmt_hr_ci(HR, LCL95, UCL95),
    `P value` = fmt_p(P)
)]

# H: all 14 LOO rows
s8h <- loo[, .(
    `LOO replicate` = LOO,
    N,
    Events = EVENTS,
    `Interaction HR (95% CI)` = fmt_hr_ci(INTERACTION_HR, LCL95, UCL95)
)]

# I: original and restricted PH diagnostics.
# Original M1/M2 values are frozen in the B1-14B audit and are reporting
# anchors, not newly selected analyses:
#   M1 global .993505; interaction .802557
#   M2 global .268108; interaction .893388
get_ph_p <- function(ph, model_id, term_pattern) {
    hit <- ph[MODEL_ID == model_id & grepl(term_pattern, TERM, fixed = FALSE)]
    if (nrow(hit) != 1L) stop("PH term not uniquely found: ", model_id, " / ", term_pattern)
    hit$P[1]
}

s8i <- data.table(
    Population = c(
        "Original 319-patient",
        "Original 319-patient",
        "Shared <=730-day, n=260",
        "Shared <=730-day, n=260"
    ),
    Model = c("M1", "M2", "M1", "M2"),
    `Global PH P` = c(
        ".994",
        ".268",
        sprintf("%.3f", get_ph_p(ph730, "M1", "GLOBAL")),
        sprintf("%.3f", get_ph_p(ph730, "M2", "GLOBAL"))
    ),
    `Treatment x bypass PH P` = c(
        ".803",
        ".893",
        sprintf("%.3f", get_ph_p(ph730, "M1", "ANTI_EGFR_TREAT_R:BYPASS_HIGH_CONFIDENCE_R")),
        sprintf("%.3f", get_ph_p(ph730, "M2", "ANTI_EGFR_TREAT_R:BYPASS_HIGH_CONFIDENCE_R"))
    ),
    Interpretation = c(
        "No evidence of global or interaction nonproportionality",
        "Global and interaction tests non-significant; active-irinotecan term P=.050",
        "No evidence of global or interaction nonproportionality",
        "No evidence of global or interaction nonproportionality"
    )
)

# ------------------------------------------------------------
# Rebuild S7 from frozen existing A-D + new E
# ------------------------------------------------------------
s7_old <- extract_docx_tables(s7_file)
if (length(s7_old) < 4L) stop("Expected at least 4 existing tables in S7.")
s7_old <- s7_old[1:4]

cat("\nExisting S7 tables recovered successfully:\n")
for (i in seq_along(s7_old)) {
    cat(
        "  S7 table ", i,
        ": ", nrow(s7_old[[i]]), " rows x ",
        ncol(s7_old[[i]]), " columns\n",
        sep = ""
    )
}

s7_tmp <- tempfile(fileext = ".docx")
doc <- read_docx()
doc <- body_add_par(doc, "Supplementary Table S7. Outcome-blind bevacizumab comparator feasibility and propensity-overlap audit", style = "heading 1")

labels7 <- c(
    "A. Genomic positivity by treatment arm",
    "B. Outcome-blind treatment-context summary",
    "C. Propensity-score overlap diagnostics",
    "D. Predeclared feasibility gate"
)
for (i in 1:4) {
    doc <- body_add_par(doc, labels7[i], style = "heading 2")
    if (i == 3) doc <- body_add_par(doc, "Propensity-score common-support interval: 0.2472 to 0.8608.")
    doc <- body_add_flextable(doc, three_line_ft(s7_old[[i]], 8))
}
doc <- body_add_par(doc, "Gate decision: GO. An explicitly exploratory comparator analysis may be estimable. This does not establish exchangeability or causal comparability. Any next-stage analysis must retain an exploratory label and use overlap-aware adjustment.")

doc <- body_add_par(doc, "E. Outcome-blind high-confidence component composition by treatment arm", style = "heading 2")
doc <- body_add_flextable(doc, three_line_ft(as.data.frame(s7e), 8))
doc <- body_add_par(doc, "Patients may carry more than one qualifying component; percentages therefore need not sum to 100%. This composition audit was outcome-blind and descriptive.")

doc <- body_add_par(doc, paste0(
    "Notes: B1-14 was completed before comparator outcomes were accessed. No rwPFS, OS, death, progression, hazard ratio, treatment effect, or treatment-by-bypass interaction was analyzed. ",
    "The frozen high-confidence genomic definition was unchanged. The outcome-blind propensity model was: I(FIRST_BIOLOGIC_R == \"anti-EGFR\") ~ GENDER_F_R + SUBSITE_F_R + LOG_STAGE4_TO_BIOLOGIC_R + PRIOR_MAJOR_N_F_R + ACTIVE_IRI_F_R + ACTIVE_OX_F_R + PANEL_OLD_F_R + SAMPLE_METASTATIC_F_R + LOG_SPECIMEN_TO_BIOLOGIC_R + MET_SITE_COUNT_90D_R. ",
    "Overlap-weight ESS denotes effective sample size under overlap weighting. GO indicates feasibility for an explicitly exploratory comparator analysis under the predeclared support criteria; it does not establish exchangeability, causal comparability, treatment predictiveness, or external validation. Common-support coverage was assessed as a feasibility diagnostic; the subsequent overlap-weighted outcome analysis retained the full propensity-score-complete population and did not trim to the empirical common-support interval."
))
doc <- body_end_section_landscape(doc)
print(doc, target = s7_tmp)
file.copy(s7_tmp, s7_file, overwrite = TRUE)
unlink(s7_tmp)

# ------------------------------------------------------------
# Rebuild S8 from frozen existing A-E + new F-I
# ------------------------------------------------------------
s8_old <- extract_docx_tables(s8_file)
if (length(s8_old) < 5L) stop("Expected at least 5 existing tables in S8.")
s8_old <- s8_old[1:5]

cat("\nExisting S8 tables recovered successfully:\n")
for (i in seq_along(s8_old)) {
    cat(
        "  S8 table ", i,
        ": ", nrow(s8_old[[i]]), " rows x ",
        ncol(s8_old[[i]]), " columns\n",
        sep = ""
    )
}

s8_tmp <- tempfile(fileext = ".docx")
doc <- read_docx()
doc <- body_add_par(doc, "Supplementary Table S8. Exploratory overlap-weighted bevacizumab comparator analysis of rwPFS", style = "heading 1")
labels8 <- c(
    "A. Treatment-by-bypass endpoint counts",
    "B. Treatment-by-bypass interaction models",
    "C. Contrasts derived from each interaction model",
    "D. Supportive within-arm bypass associations using shared covariates",
    "E. Propensity-score balance"
)
for (i in 1:5) {
    doc <- body_add_par(doc, labels8[i], style = "heading 2")
    doc <- body_add_flextable(doc, three_line_ft(s8_old[[i]], 8))
}

doc <- body_add_par(doc, "F. Shared <=730-day restricted comparator population and interaction models", style = "heading 2")
doc <- body_add_flextable(doc, three_line_ft(as.data.frame(s8f1), 8))
doc <- body_add_flextable(doc, three_line_ft(as.data.frame(s8f2), 8))

doc <- body_add_par(doc, "G. Contrasts from the shared <=730-day restricted M1/M2 models", style = "heading 2")
doc <- body_add_flextable(doc, three_line_ft(as.data.frame(s8g), 8))
doc <- body_add_par(doc, "The restricted comparator within-arm HRs are overlap-weighted estimates for the comparator overlap population and use covariate sets different from the anti-EGFR-only Model E/PH-addressed analyses in Tables S5-S6; these estimates are not interchangeable.")

doc <- body_add_par(doc, "H. Bevacizumab bypass-positive leave-one-out full-pipeline influence analysis", style = "heading 2")
doc <- body_add_flextable(doc, three_line_ft(as.data.frame(s8h), 8))
doc <- body_add_par(doc, sprintf(
    "Summary: 14/14 refits succeeded; interaction HR range, %.2f–%.2f; 14/14 estimates were >1 and 14/14 lower 95%% confidence limits were >1. Each replicate refit the propensity model, recalculated overlap weights, and refit M1.",
    min(loo$INTERACTION_HR), max(loo$INTERACTION_HR)
))

doc <- body_add_par(doc, "I. Proportional-hazards diagnostics for original and shared <=730-day comparator models", style = "heading 2")
doc <- body_add_flextable(doc, three_line_ft(as.data.frame(s8i), 8))

doc <- body_add_par(doc, paste0(
    "Notes: This is an exploratory comparator specificity analysis, not a primary causal comparative-effectiveness analysis. M1 used propensity-score overlap weights fixed in the outcome-blind B1-14 workflow; M2 additionally adjusted for sex, frozen primary subsite, Stage IV-to-biologic interval, prior exposure to 0/1/2 major cytotoxic classes, and active irinotecan/oxaliplatin context. Robust sandwich variance was used for weighted Cox models. ",
    "The interaction HR is the ratio of the anti-EGFR-versus-bevacizumab HR in bypass-positive patients to the corresponding treatment HR in bypass-negative patients. Maximum absolute overlap-weighted SMD before outcome access was 0.000 after rounding. Panels F-I are post-hoc (post-outcome) diagnostics from B1-16 and do not redefine the primary exposure or the original comparator model. ",
    "The shared <=730-day analysis refit the unchanged propensity model and overlap weights after applying the same specimen-age restriction to both arms. Within-bevacizumab bypass HRs were imprecise and compatible with no association; values below 1 must not be interpreted as evidence of a protective effect. The restricted overlap-weighted within-arm HRs target a different population and use different adjustment sets from the anti-EGFR-only Model E/PH-addressed estimates in Tables S5-S6 and are not directly interchangeable. This analysis remains observational, exploratory, and noncausal."
))
doc <- body_end_section_landscape(doc)
print(doc, target = s8_tmp)
file.copy(s8_tmp, s8_file, overwrite = TRUE)
unlink(s8_tmp)

# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------
audit_lines <- c(
    "B1-18 REPORTING-ONLY S7/S8 REFRESH",
    "===================================",
    "",
    "NO NEW OUTCOME MODEL WAS FIT.",
    "NO FIGURE WAS REGENERATED.",
    "",
    sprintf("B1-16 restricted M1 interaction HR: %.6f (%.6f–%.6f), P=%.6g",
            m1$INTERACTION_HR, m1$LCL95, m1$UCL95, m1$P),
    sprintf("B1-16 restricted M2 interaction HR: %.6f (%.6f–%.6f), P=%.6g",
            m2$INTERACTION_HR, m2$LCL95, m2$UCL95, m2$P),
    sprintf("LOO successful: %d/%d", sum(loo$STATUS == "SUCCESS"), nrow(loo)),
    sprintf("LOO interaction HR range: %.6f–%.6f", min(loo$INTERACTION_HR), max(loo$INTERACTION_HR)),
    sprintf("LOO lower 95%% CI >1: %d/%d", sum(loo$LCL95 > 1), nrow(loo)),
    "",
    "Original comparator PH reporting anchors from frozen B1-14B audit:",
    "M1 global P=.993505; treatment-by-bypass PH P=.802557",
    "M2 global P=.268108; treatment-by-bypass PH P=.893388",
    "",
    paste("Overwrote:", s7_file),
    paste("Overwrote:", s8_file)
)
writeLines(audit_lines, audit_file)

cat("\nB1-18 completed.\n")
cat("Overwrote:\n", s7_file, "\n", s8_file, "\n", sep = "")
cat("Audit:\n", audit_file, "\n", sep = "")
