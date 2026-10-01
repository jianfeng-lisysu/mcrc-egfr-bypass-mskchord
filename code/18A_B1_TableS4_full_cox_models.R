# ============================================================
# A7 / B-line
# 18A_B1_TableS4_full_cox_models.R
#
# PURPOSE
#   Re-run the frozen B1-07 analysis locally, identify the exact
#   primary adjusted rwPFS and secondary adjusted OS Cox models
#   from their hard anchors, and export the COMPLETE coefficient
#   table as Supplementary Table S4.
#
# IMPORTANT
#   - No new analysis is introduced.
#   - No covariate is added/removed.
#   - The B1-07 frozen model is re-used exactly.
#   - The script HARD-STOPS if the reproduced bypass HRs do not
#     match the frozen B1-07 anchors.
#   - Same output filename is overwritten on rerun.
#
# INPUT
#   02_code/07_B1_07_first_unblinded_analysis.R
#
# OUTPUTS
#   07_tables/TableS4_full_multivariable_cox_models.docx
#   06_logs_and_audit/TableS4_full_multivariable_cox_models_audit.txt
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

source_file <- file.path(
    project_root,
    "02_code",
    "07_B1_07_first_unblinded_analysis.R"
)

table_dir <- file.path(project_root, "07_tables")
audit_dir <- file.path(project_root, "06_logs_and_audit")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

word_file <- file.path(
    table_dir,
    "TableS4_full_multivariable_cox_models.docx"
)

audit_file <- file.path(
    audit_dir,
    "TableS4_full_multivariable_cox_models_audit.txt"
)

if (!file.exists(source_file)) {
    stop(
        "Missing B1-07 source script:\n",
        source_file
    )
}

required_packages <- c(
    "survival",
    "data.table",
    "officer",
    "flextable"
)

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

library(survival)
library(data.table)
library(officer)
library(flextable)

cat("\n============================================================\n")
cat("B1-18A TABLE S4: FULL FROZEN COX MODELS\n")
cat("============================================================\n\n")


# ============================================================
# 1. Re-run B1-07 in an isolated environment
# ============================================================

b107_env <- new.env(parent = globalenv())

sys.source(
    source_file,
    envir = b107_env
)

cat("B1-07 source completed.\n")


# ============================================================
# 2. Recursively collect every coxph object produced by B1-07
# ============================================================

collect_coxph <- function(x, path = "root", depth = 0L) {

    out <- list()

    if (inherits(x, "coxph")) {
        out[[path]] <- x
        return(out)
    }

    if (depth >= 6L) {
        return(out)
    }

    if (is.list(x) && !inherits(x, c("data.frame", "data.table"))) {

        nms <- names(x)

        for (i in seq_along(x)) {

            child_name <- if (
                !is.null(nms) &&
                length(nms) >= i &&
                nzchar(nms[i])
            ) {
                nms[i]
            } else {
                paste0("[[", i, "]]")
            }

            child_path <- paste0(
                path,
                "$",
                child_name
            )

            out <- c(
                out,
                collect_coxph(
                    x[[i]],
                    child_path,
                    depth + 1L
                )
            )
        }
    }

    out
}

cox_objects <- list()

object_names <- ls(
    envir = b107_env,
    all.names = TRUE
)

for (nm in object_names) {

    obj <- get(
        nm,
        envir = b107_env
    )

    cox_objects <- c(
        cox_objects,
        collect_coxph(
            obj,
            path = nm
        )
    )
}

if (length(cox_objects) == 0L) {
    stop(
        "No coxph objects were found after sourcing B1-07."
    )
}

cat(
    "Number of coxph objects found: ",
    length(cox_objects),
    "\n",
    sep = ""
)


# ============================================================
# 3. Identify the exact frozen adjusted rwPFS and OS models
#    using HARD numerical anchors
# ============================================================

get_bypass_hr <- function(fit) {

    b <- coef(fit)

    if (is.null(b)) {
        return(NA_real_)
    }

    term_names <- names(b)

    idx <- which(
        term_names ==
            "BYPASS_HIGH_CONFIDENCE_R"
    )

    if (length(idx) != 1L) {
        return(NA_real_)
    }

    exp(
        unname(
            b[idx]
        )
    )
}

model_inventory <- rbindlist(
    lapply(
        names(cox_objects),
        function(nm) {

            fit <- cox_objects[[nm]]

            data.table(
                OBJECT = nm,
                N = ifelse(
                    is.null(fit$n),
                    NA_integer_,
                    as.integer(fit$n)
                ),
                EVENTS = ifelse(
                    is.null(fit$nevent),
                    NA_integer_,
                    as.integer(fit$nevent)
                ),
                BYPASS_HR =
                    get_bypass_hr(fit),
                N_COEFFICIENTS =
                    length(coef(fit))
            )
        }
    ),
    fill = TRUE
)

rw_candidates <- model_inventory[
    N == 191L &
    EVENTS == 183L &
    is.finite(BYPASS_HR) &
    abs(BYPASS_HR - 2.334372) < 0.0001
]

os_candidates <- model_inventory[
    N == 191L &
    EVENTS == 129L &
    is.finite(BYPASS_HR) &
    abs(BYPASS_HR - 1.474230) < 0.0001
]

if (nrow(rw_candidates) != 1L) {
    cat("\nCandidate rwPFS models:\n")
    print(rw_candidates)

    stop(
        "Could not uniquely identify the frozen adjusted rwPFS model."
    )
}

if (nrow(os_candidates) != 1L) {
    cat("\nCandidate OS models:\n")
    print(os_candidates)

    stop(
        "Could not uniquely identify the frozen adjusted OS model."
    )
}

rw_name <- rw_candidates$OBJECT[1]
os_name <- os_candidates$OBJECT[1]

fit_rw <- cox_objects[[rw_name]]
fit_os <- cox_objects[[os_name]]

rw_anchor_hr <- get_bypass_hr(fit_rw)
os_anchor_hr <- get_bypass_hr(fit_os)

if (abs(rw_anchor_hr - 2.334372) >= 0.0001) {
    stop(
        "rwPFS hard anchor failed."
    )
}

if (abs(os_anchor_hr - 1.474230) >= 0.0001) {
    stop(
        "OS hard anchor failed."
    )
}

cat(
    "rwPFS frozen model: ",
    rw_name,
    "\n",
    sep = ""
)

cat(
    "OS frozen model: ",
    os_name,
    "\n",
    sep = ""
)

cat(
    "rwPFS bypass HR = ",
    sprintf("%.6f", rw_anchor_hr),
    " [PASS]\n",
    sep = ""
)

cat(
    "OS bypass HR = ",
    sprintf("%.6f", os_anchor_hr),
    " [PASS]\n",
    sep = ""
)


# ============================================================
# 4. Extract COMPLETE coefficient tables
# ============================================================

extract_full_model <- function(
    fit,
    outcome_label
) {

    s <- summary(fit)

    cf <- as.data.table(
        s$coefficients,
        keep.rownames = "TERM"
    )

    ci <- as.data.table(
        s$conf.int,
        keep.rownames = "TERM"
    )

    setnames(
        cf,
        old = names(cf),
        new = make.names(
            names(cf),
            unique = TRUE
        )
    )

    setnames(
        ci,
        old = names(ci),
        new = make.names(
            names(ci),
            unique = TRUE
        )
    )

    # Coefficient table column names can vary slightly by survival version.
    coef_col <- grep(
        "^coef$",
        names(cf),
        value = TRUE
    )

    p_col <- grep(
        "Pr",
        names(cf),
        value = TRUE
    )

    hr_col <- grep(
        "^exp.coef",
        names(ci),
        value = TRUE
    )

    lower_col <- grep(
        "lower",
        names(ci),
        value = TRUE
    )

    upper_col <- grep(
        "upper",
        names(ci),
        value = TRUE
    )

    if (
        length(coef_col) != 1L ||
        length(p_col) != 1L ||
        length(hr_col) != 1L ||
        length(lower_col) != 1L ||
        length(upper_col) != 1L
    ) {
        stop(
            "Unexpected coxph summary column structure."
        )
    }

    out <- merge(
        cf[
            ,
            .(
                TERM,
                BETA =
                    get(coef_col),
                P =
                    get(p_col)
            )
        ],
        ci[
            ,
            .(
                TERM,
                HR =
                    get(hr_col),
                LCL95 =
                    get(lower_col),
                UCL95 =
                    get(upper_col)
            )
        ],
        by = "TERM",
        all = TRUE,
        sort = FALSE
    )

    out[
        ,
        OUTCOME := outcome_label
    ]

    out
}

rw_full <- extract_full_model(
    fit_rw,
    "rwPFS"
)

os_full <- extract_full_model(
    fit_os,
    "OS"
)


# ============================================================
# 5. Harmonize rows across rwPFS and OS
# ============================================================

pretty_term <- function(x) {

    z <- x

    z <- sub(
        "^BYPASS_HIGH_CONFIDENCE_R$",
        "High-confidence EGFR-bypass positive",
        z
    )

    z <- sub(
        "^GENDER_F_R",
        "Sex: ",
        z
    )

    z <- sub(
        "^SUBSITE_F_R",
        "Primary subsite: ",
        z
    )

    z <- sub(
        "^ANTI_EGFR_AGENT_F_R",
        "Anti-EGFR agent: ",
        z
    )

    z <- sub(
        "^LOG_STAGE4_TO_ANTI_EGFR_R$",
        "log(Stage IV-to-anti-EGFR interval)",
        z
    )

    z <- sub(
        "^PRIOR_HEAVY_F_R",
        "Prior oxaliplatin + irinotecan exposure: ",
        z
    )

    z
}

all_terms <- unique(
    c(
        rw_full$TERM,
        os_full$TERM
    )
)

combined <- data.table(
    TERM = all_terms
)

combined <- merge(
    combined,
    rw_full[
        ,
        .(
            TERM,
            RW_BETA = BETA,
            RW_HR = HR,
            RW_LCL95 = LCL95,
            RW_UCL95 = UCL95,
            RW_P = P
        )
    ],
    by = "TERM",
    all.x = TRUE,
    sort = FALSE
)

combined <- merge(
    combined,
    os_full[
        ,
        .(
            TERM,
            OS_BETA = BETA,
            OS_HR = HR,
            OS_LCL95 = LCL95,
            OS_UCL95 = UCL95,
            OS_P = P
        )
    ],
    by = "TERM",
    all.x = TRUE,
    sort = FALSE
)

combined[
    ,
    COVARIATE := vapply(
        TERM,
        pretty_term,
        character(1)
    )
]


# ============================================================
# 6. Formatting helpers
# ============================================================

fmt_beta <- function(x) {

    ifelse(
        is.finite(x),
        sprintf("%.3f", x),
        "—"
    )
}

fmt_hr_ci <- function(
    hr,
    lo,
    hi
) {

    ifelse(
        is.finite(hr) &
        is.finite(lo) &
        is.finite(hi),
        sprintf(
            "%.2f (%.2f–%.2f)",
            hr,
            lo,
            hi
        ),
        "—"
    )
}

fmt_p <- function(x) {

    ifelse(
        !is.finite(x),
        "—",
        ifelse(
            x < 0.001,
            "<0.001",
            sprintf("%.3f", x)
        )
    )
}

table_data <- combined[
    ,
    .(
        Covariate = COVARIATE,
        `rwPFS β` =
            fmt_beta(RW_BETA),
        `rwPFS HR (95% CI)` =
            fmt_hr_ci(
                RW_HR,
                RW_LCL95,
                RW_UCL95
            ),
        `rwPFS P` =
            fmt_p(RW_P),
        `OS β` =
            fmt_beta(OS_BETA),
        `OS HR (95% CI)` =
            fmt_hr_ci(
                OS_HR,
                OS_LCL95,
                OS_UCL95
            ),
        `OS P` =
            fmt_p(OS_P)
    )
]


# ============================================================
# 7. Generate publication-style Word table
# ============================================================

ft <- flextable(
    table_data
)

ft <- theme_booktabs(ft)

ft <- font(
    ft,
    fontname = "Arial",
    part = "all"
)

ft <- fontsize(
    ft,
    size = 9,
    part = "all"
)

ft <- bold(
    ft,
    part = "header"
)

ft <- align(
    ft,
    j = 1,
    align = "left",
    part = "all"
)

ft <- align(
    ft,
    j = 2:ncol(table_data),
    align = "center",
    part = "all"
)

ft <- valign(
    ft,
    valign = "center",
    part = "all"
)

ft <- padding(
    ft,
    padding.top = 2,
    padding.bottom = 2,
    padding.left = 3,
    padding.right = 3,
    part = "all"
)

ft <- width(
    ft,
    j = 1,
    width = 2.7
)

ft <- width(
    ft,
    j = 2,
    width = 0.75
)

ft <- width(
    ft,
    j = 3,
    width = 1.45
)

ft <- width(
    ft,
    j = 4,
    width = 0.70
)

ft <- width(
    ft,
    j = 5,
    width = 0.75
)

ft <- width(
    ft,
    j = 6,
    width = 1.45
)

ft <- width(
    ft,
    j = 7,
    width = 0.70
)

ft <- autofit(ft)

doc <- read_docx()

doc <- body_add_par(
    doc,
    "Supplementary Table S4. Full multivariable Cox model coefficients for the primary rwPFS and secondary overall-survival analyses",
    style = "heading 1"
)

doc <- body_add_par(
    doc,
    paste0(
        "Primary rwPFS model: N = ",
        fit_rw$n,
        ", events = ",
        fit_rw$nevent,
        ". Secondary OS model: N = ",
        fit_os$n,
        ", deaths = ",
        fit_os$nevent,
        "."
    )
)

doc <- body_add_flextable(
    doc,
    value = ft
)

doc <- body_add_par(
    doc,
    paste(
        "Notes:",
        "This table reports every coefficient from the exact frozen B1-07 adjusted Cox models.",
        "No covariate was added, removed, or recoded for Table S4.",
        "Factor coefficients are shown relative to the reference levels used in the frozen B1-07 model.",
        "For log(Stage IV-to-anti-EGFR interval), the HR is per 1-unit increase on the natural-log scale.",
        "rwPFS denotes real-world progression-free survival; OS, overall survival; HR, hazard ratio; CI, confidence interval."
    )
)

sec <- prop_section(
    page_size = page_size(
        orient = "landscape"
    ),
    page_margins = page_mar(
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

print(
    doc,
    target = word_file
)


# ============================================================
# 8. Audit
# ============================================================

audit_lines <- c(
    "A7 B1-18A TABLE S4 AUDIT",
    "========================",
    "",
    paste0(
        "B1-07 source: ",
        source_file
    ),
    "",
    paste0(
        "rwPFS model object: ",
        rw_name
    ),
    paste0(
        "rwPFS N/events: ",
        fit_rw$n,
        "/",
        fit_rw$nevent
    ),
    paste0(
        "rwPFS bypass HR reproduced: ",
        sprintf("%.8f", rw_anchor_hr),
        " [expected 2.334372; PASS]"
    ),
    "",
    paste0(
        "OS model object: ",
        os_name
    ),
    paste0(
        "OS N/events: ",
        fit_os$n,
        "/",
        fit_os$nevent
    ),
    paste0(
        "OS bypass HR reproduced: ",
        sprintf("%.8f", os_anchor_hr),
        " [expected 1.474230; PASS]"
    ),
    "",
    "FULL rwPFS COEFFICIENTS",
    "-----------------------",
    capture.output(
        print(rw_full)
    ),
    "",
    "FULL OS COEFFICIENTS",
    "--------------------",
    capture.output(
        print(os_full)
    ),
    "",
    paste0(
        "Saved Word table: ",
        word_file
    )
)

writeLines(
    audit_lines,
    audit_file
)

cat("\n============================================================\n")
cat("TABLE S4 COMPLETE\n")
cat("============================================================\n")

cat(
    "\nSaved Word table:\n",
    word_file,
    "\n",
    sep = ""
)

cat(
    "\nSaved audit:\n",
    audit_file,
    "\n",
    sep = ""
)
