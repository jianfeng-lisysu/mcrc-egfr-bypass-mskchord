# ============================================================
# A7 / B-line
# 12_B1_12_finalize_manuscript_tables.R
#
# PURPOSE
#   Final manuscript-table cleanup only.
#   NO new primary analysis.
#   NO exposure redefinition.
#
# OVERWRITES:
#   07_tables/Table1_baseline_characteristics.docx
#   07_tables/Table2_primary_secondary_results.docx
#   07_tables/TableS1_bypass_component_definitions.docx
#   07_tables/TableS2_robustness_results.docx
#   07_tables/TableS3_TTNTD_validation.docx
#
# KEY FIXES
#   - remove automatic "1." numbering from titles
#   - clean manuscript-facing terminology
#   - Table 1: remove redundant Level column
#   - Table S1: landscape layout, add observed carrier counts
#   - Table S2: remove internal _AMP/_LOF/_HOTSPOT labels
#               and include available locked sensitivity outputs
#   - Table S3: replace internal "PH-addressed" wording with
#               manuscript-facing "Adjusted Cox, stratified"
#   - all files remain Word three-line tables
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

intermediate_dir <- file.path(project_root, "03_intermediate")
results_dir      <- file.path(project_root, "04_results")
table_dir        <- file.path(project_root, "07_tables")
audit_dir        <- file.path(project_root, "06_logs_and_audit")

if (!dir.exists(table_dir)) {
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(audit_dir)) {
    dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
}

locked_file <- file.path(
    intermediate_dir,
    "A7_B1_locked_analysis.rds"
)

b107_file <- file.path(
    results_dir,
    "B1_07_primary_models.rds"
)

b108_file <- file.path(
    results_dir,
    "B1_08_robustness_validation.rds"
)

b110_file <- file.path(
    results_dir,
    "B1_10_refined_TTNTD_PH_validation.rds"
)

for (f in c(
    locked_file,
    b107_file,
    b108_file,
    b110_file
)) {
    if (!file.exists(f)) {
        stop(
            "Missing required input:\n",
            f
        )
    }
}

packages <- c(
    "data.table",
    "officer",
    "flextable"
)

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

library(data.table)
library(officer)
library(flextable)

cat("\n============================================================\n")
cat("B1-12 FINALIZE MANUSCRIPT TABLES\n")
cat("============================================================\n\n")


# ============================================================
# 1. Read locked objects
# ============================================================

dt   <- as.data.table(readRDS(locked_file))
b107 <- readRDS(b107_file)
b108 <- readRDS(b108_file)
b110 <- readRDS(b110_file)

if (
    nrow(dt) != 191L ||
    dt[BYPASS_HIGH_CONFIDENCE_R == 1L, .N] != 23L
) {
    stop("Locked cohort anchors failed.")
}


# ============================================================
# 2. General helpers
# ============================================================

fmt_num <- function(x, digits = 2) {

    vapply(
        x,
        function(z) {

            if (
                length(z) == 0L ||
                is.na(z) ||
                !is.finite(z)
            ) {
                return("NA")
            }

            formatC(
                z,
                format = "f",
                digits = digits
            )
        },
        character(1)
    )
}

fmt_p <- function(x) {

    vapply(
        x,
        function(z) {

            if (
                length(z) == 0L ||
                is.na(z) ||
                !is.finite(z)
            ) {
                return("NA")
            }

            if (z < 0.001) {
                return("<0.001")
            }

            formatC(
                z,
                format = "f",
                digits = 3
            )
        },
        character(1)
    )
}

median_iqr <- function(x, digits = 1) {

    z <- suppressWarnings(
        as.numeric(x)
    )

    z <- z[
        is.finite(z)
    ]

    if (length(z) == 0L) {
        return("NA")
    }

    q <- quantile(
        z,
        probs = c(
            0.25,
            0.50,
            0.75
        ),
        na.rm = TRUE,
        names = FALSE
    )

    paste0(
        fmt_num(q[2], digits),
        " [",
        fmt_num(q[1], digits),
        ", ",
        fmt_num(q[3], digits),
        "]"
    )
}

n_pct <- function(n, den, digits = 1) {

    paste0(
        n,
        " (",
        fmt_num(
            100 * n / den,
            digits
        ),
        "%)"
    )
}

get_first_existing <- function(
    x,
    candidates
) {

    hit <- candidates[
        candidates %in%
            names(x)
    ]

    if (length(hit) == 0L) {
        return(NULL)
    }

    hit[1]
}

pretty_analysis_label <- function(x) {

    z <- as.character(x)

    replacements <- c(
        "ERBB2_AMP_R" =
            "ERBB2 amplification",

        "ERBB2_AMP" =
            "ERBB2 amplification",

        "ERBB2_ACTIVATING_R" =
            "ERBB2 activating mutation",

        "ERBB2_ACTIVATING" =
            "ERBB2 activating mutation",

        "MET_AMP_R" =
            "MET amplification",

        "MET_AMP" =
            "MET amplification",

        "PIK3CA_EX20_R" =
            "PIK3CA exon 20",

        "PIK3CA_EX20" =
            "PIK3CA exon 20",

        "PTEN_LOF_R" =
            "PTEN loss-of-function",

        "PTEN_LOF" =
            "PTEN loss-of-function",

        "KRAS_AMP_R" =
            "KRAS amplification",

        "KRAS_AMP" =
            "KRAS amplification",

        "MAP2K1_HOTSPOT_R" =
            "MAP2K1 hotspot",

        "MAP2K1_HOTSPOT" =
            "MAP2K1 hotspot",

        "NF1_LOF_R" =
            "NF1 loss-of-function",

        "NF1_LOF" =
            "NF1 loss-of-function",

        "BYPASS_PRESSING_STRICT_R" =
            "Strict PRESSING-informed composite",

        "BYPASS_HIGH_CONFIDENCE_NO_D67E_R" =
            "High-confidence composite excluding MAP2K1 D67E",

        "BYPASS_ORIGINAL_BROAD_R" =
            "Original broad composite"
    )

    for (old in names(replacements)) {

        z <- gsub(
            old,
            replacements[[old]],
            z,
            fixed = TRUE
        )
    }

    z <- gsub(
        "_",
        " ",
        z,
        fixed = TRUE
    )

    z
}


# ============================================================
# 3. Flextable / Word helpers
# ============================================================

make_three_line_ft <- function(
    df,
    font_size = 9,
    col_widths = NULL,
    align_first = "left"
) {

    ft <- flextable(
        df
    )

    ft <- border_remove(
        ft
    )

    border_top <- fp_border(
        color = "black",
        width = 1.0
    )

    border_mid <- fp_border(
        color = "black",
        width = 0.7
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
        align = align_first,
        part = "all"
    )

    ft <- valign(
        ft,
        valign = "center",
        part = "all"
    )

    ft <- padding(
        ft,
        padding.top = 3,
        padding.bottom = 3,
        padding.left = 3,
        padding.right = 3,
        part = "all"
    )

    ft <- set_table_properties(
        ft,
        layout = "fixed",
        width = 1
    )

    if (!is.null(col_widths)) {

        if (
            length(col_widths) !=
                ncol(df)
        ) {
            stop(
                "col_widths length does not match table columns."
            )
        }

        for (j in seq_along(col_widths)) {

            ft <- width(
                ft,
                j = j,
                width = col_widths[j]
            )
        }
    }

    ft
}


write_table_docx <- function(
    df,
    title,
    filename,
    footnotes = NULL,
    orientation = c(
        "portrait",
        "landscape"
    ),
    font_size = 9,
    col_widths = NULL
) {

    orientation <- match.arg(
        orientation
    )

    doc <- read_docx()

    section <- prop_section(
        page_size =
            page_size(
                orient = orientation
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
        section
    )

    # Important: use normal paragraph, not Heading 1,
    # to prevent automatic "1." numbering.
    doc <- body_add_fpar(
        doc,
        fpar(
            ftext(
                title,
                fp_text(
                    font.family = "Arial",
                    font.size = 11,
                    bold = TRUE
                )
            )
        )
    )

    doc <- body_add_par(
        doc,
        ""
    )

    ft <- make_three_line_ft(
        df,
        font_size = font_size,
        col_widths = col_widths
    )

    doc <- body_add_flextable(
        doc,
        ft
    )

    if (!is.null(footnotes)) {

        doc <- body_add_par(
            doc,
            ""
        )

        for (note in footnotes) {

            doc <- body_add_fpar(
                doc,
                fpar(
                    ftext(
                        note,
                        fp_text(
                            font.family = "Arial",
                            font.size = 8
                        )
                    )
                )
            )
        }
    }

    print(
        doc,
        target =
            file.path(
                table_dir,
                filename
            )
    )
}


# ============================================================
# 4. TABLE 1
# ============================================================

dt[
    ,
    BYPASS_GROUP_TABLE :=
        fifelse(
            BYPASS_HIGH_CONFIDENCE_R == 0L,
            "Bypass−",
            "Bypass+"
        )
]

group_levels <- c(
    "Bypass−",
    "Bypass+"
)

group_n <- dt[
    ,
    .N,
    by = BYPASS_GROUP_TABLE
]

get_group_n <- function(g) {

    group_n[
        BYPASS_GROUP_TABLE == g,
        N
    ]
}

binary_row <- function(
    label,
    condition_fun
) {

    vals <- vapply(
        group_levels,
        function(g) {

            d0 <- dt[
                BYPASS_GROUP_TABLE == g
            ]

            n0 <- sum(
                condition_fun(d0),
                na.rm = TRUE
            )

            n_pct(
                n0,
                nrow(d0)
            )
        },
        character(1)
    )

    data.table(
        Characteristic = label,
        `Bypass− (n=168)` = vals[1],
        `Bypass+ (n=23)` = vals[2]
    )
}

continuous_row <- function(
    label,
    variable
) {

    vals <- vapply(
        group_levels,
        function(g) {

            d0 <- dt[
                BYPASS_GROUP_TABLE == g
            ]

            median_iqr(
                d0[[variable]],
                digits = 1
            )
        },
        character(1)
    )

    data.table(
        Characteristic = label,
        `Bypass− (n=168)` = vals[1],
        `Bypass+ (n=23)` = vals[2]
    )
}

table1 <- rbindlist(
    list(
        binary_row(
            "Male sex",
            function(d0) d0$GENDER == "Male"
        ),

        binary_row(
            "Rectum/rectosigmoid primary",
            function(d0)
                d0$SUBSITE_FROZEN_R ==
                    "rectum/rectosigmoid"
        ),

        binary_row(
            "Panitumumab",
            function(d0)
                d0$ANTI_EGFR_AGENT_R ==
                    "PANITUMUMAB"
        ),

        binary_row(
            "Prior exposure to both oxaliplatin and irinotecan",
            function(d0)
                d0$PRIOR_HEAVY_R == 1L
        ),

        binary_row(
            "Metastatic tumor sample",
            function(d0)
                grepl(
                    "metast",
                    as.character(
                        d0$SAMPLE_TYPE
                    ),
                    ignore.case = TRUE
                )
        ),

        binary_row(
            "MSK-IMPACT 341/410 panel",
            function(d0)
                d0$PANEL_OLD_R == 1L
        ),

        continuous_row(
            "Stage IV diagnosis to first anti-EGFR, days",
            "STAGE4_TO_ANTI_EGFR_DAYS_R"
        ),

        continuous_row(
            "Specimen acquisition to first anti-EGFR, days",
            "SPECIMEN_TO_ANTI_EGFR_DAYS_R"
        ),

        continuous_row(
            "Sequencing to first anti-EGFR, days",
            "SEQUENCING_TO_ANTI_EGFR_DAYS_R"
        ),

        continuous_row(
            "Recent metastatic-site count",
            "MET_SITE_COUNT_90D_R"
        )
    ),
    fill = TRUE
)

site_avail <- dt[
    ,
    .(
        available =
            sum(
                is.finite(
                    as.numeric(
                        MET_SITE_COUNT_90D_R
                    )
                )
            )
    ),
    by =
        BYPASS_GROUP_TABLE
]

table1_foot <- c(
    "Categorical variables are presented as n (%); continuous variables are presented as median [IQR].",
    "Table 1 is descriptive; no between-group significance testing was performed.",
    "Age at the anti-EGFR index is unavailable in MSK-CHORD and is therefore not reported."
)

if (
    any(
        site_avail$available <
            c(
                get_group_n("Bypass−"),
                get_group_n("Bypass+")
            )
    )
) {

    table1_foot <- c(
        table1_foot,
        paste0(
            "Recent metastatic-site count was available for ",
            site_avail[
                BYPASS_GROUP_TABLE == "Bypass−",
                available
            ],
            "/168 bypass− and ",
            site_avail[
                BYPASS_GROUP_TABLE == "Bypass+",
                available
            ],
            "/23 bypass+ patients."
        )
    )
}

write_table_docx(
    table1,
    "Table 1. Baseline characteristics of the locked anti-EGFR cohort",
    "Table1_baseline_characteristics.docx",
    footnotes =
        table1_foot,
    orientation =
        "portrait",
    font_size =
        8.8,
    col_widths =
        c(
            3.8,
            1.75,
            1.75
        )
)


# ============================================================
# 5. TABLE 2
# ============================================================

primary_results <- as.data.table(
    b107$primary_results
)

table2_a <- primary_results[
    ,
    .(
        Endpoint = OUTCOME,
        Analysis = ANALYSIS,
        N,
        Events = EVENTS,

        `HR (95% CI)` =
            paste0(
                fmt_num(
                    HR,
                    2
                ),
                " (",
                fmt_num(
                    LCL95,
                    2
                ),
                "–",
                fmt_num(
                    UCL95,
                    2
                ),
                ")"
            ),

        `P value` =
            fmt_p(P)
    )
]

tt_models <- as.data.table(
    b110$model_results
)

tt_main <- tt_models[
    ENDPOINT %chin%
        c(
            "Conservative CRC TTNTD30",
            "Conservative CRC TTNTD60"
        )
]

tt_main[
    ,
    Endpoint :=
        fcase(
            ENDPOINT ==
                "Conservative CRC TTNTD30",
            "CRC-compatible TTNTD (30-day grace)",

            ENDPOINT ==
                "Conservative CRC TTNTD60",
            "CRC-compatible TTNTD (60-day grace)",

            default =
                ENDPOINT
        )
]

tt_main[
    ,
    Analysis :=
        fcase(
            MODEL ==
                "Standard adjusted",
            "Adjusted Cox",

            MODEL ==
                "PH-addressed stratified",
            "Adjusted Cox, stratified",

            default =
                MODEL
        )
]

table2_b <- tt_main[
    ,
    .(
        Endpoint,
        Analysis,
        N,
        Events = EVENTS,

        `HR (95% CI)` =
            paste0(
                fmt_num(
                    HR,
                    2
                ),
                " (",
                fmt_num(
                    LCL95,
                    2
                ),
                "–",
                fmt_num(
                    UCL95,
                    2
                ),
                ")"
            ),

        `P value` =
            fmt_p(P)
    )
]

table2 <- rbindlist(
    list(
        table2_a,
        table2_b
    ),
    fill = TRUE
)

write_table_docx(
    table2,
    "Table 2. Primary, secondary, and treatment-based outcome analyses",
    "Table2_primary_secondary_results.docx",
    footnotes = c(
        "Hazard ratios >1 indicate worse outcomes for patients with high-confidence EGFR-bypass alterations.",
        "rwPFS was the primary endpoint; OS was secondary.",
        "TTNTD denotes time to next treatment or death. The 30- and 60-day values refer to the prespecified regimen-assembly grace windows used for the treatment-based endpoint.",
        "Adjusted Cox models used the locked primary covariate set. Stratified Cox models stratified by primary subsite and prior heavy cytotoxic exposure."
    ),
    orientation =
        "portrait",
    font_size =
        8.6,
    col_widths =
        c(
            2.15,
            1.75,
            0.55,
            0.65,
            1.55,
            0.75
        )
)


# ============================================================
# 6. TABLE S1
# ============================================================

count_flag <- function(
    candidates
) {

    hit <- candidates[
        candidates %in%
            names(dt)
    ]

    if (length(hit) == 0L) {
        return(NA_integer_)
    }

    sum(
        dt[[hit[1]]] == 1L,
        na.rm = TRUE
    )
}

braf_candidates <- grep(
    "BRAF.*III|BRAF.*CLASS",
    names(dt),
    value = TRUE,
    ignore.case = TRUE
)

braf_count <- if (
    length(braf_candidates) > 0L
) {
    sum(
        dt[[braf_candidates[1]]] == 1L,
        na.rm = TRUE
    )
} else {
    NA_integer_
}

tableS1 <- data.table(
    `Definition tier` = c(
        rep(
            "Literature-anchored strict core",
            7
        ),
        rep(
            "Primary high-confidence extension",
            3
        ),
        "Sensitivity only",
        "Sensitivity only"
    ),

    `Component / rule` = c(
        "ERBB2 amplification",
        "High-confidence activating ERBB2 mutation",
        "MET amplification",
        "PIK3CA exon 20 / kinase-domain hotspot",
        "PTEN clear loss-of-function",
        "AKT1 E17K",
        "Canonical ALK/ROS1/NTRK1-3/RET fusion",
        "KRAS amplification",
        "MAP2K1 activating/resistance hotspot",
        "NF1 clear loss-of-function",
        "Original broad composite",
        "BRAF class III exclusion sensitivity"
    ),

    `Observed carriers, n` = c(
        count_flag("ERBB2_AMP_R"),
        count_flag("ERBB2_ACTIVATING_R"),
        count_flag("MET_AMP_R"),
        count_flag("PIK3CA_EX20_R"),
        count_flag("PTEN_LOF_R"),
        count_flag("AKT1_E17K_R"),
        count_flag("CANONICAL_PRESSING_FUSION_R"),
        count_flag("KRAS_AMP_R"),
        count_flag("MAP2K1_HOTSPOT_R"),
        count_flag("NF1_LOF_R"),
        count_flag("BYPASS_ORIGINAL_BROAD_R"),
        braf_count
    ),

    `Role in analysis` = c(
        rep(
            "Included in strict and primary high-confidence definitions",
            7
        ),
        rep(
            "Included in primary high-confidence definition",
            3
        ),
        "Broad-composite sensitivity analysis only",
        "Exclusion sensitivity; not classified as bypass solely by BRAF"
    )
)

write_table_docx(
    tableS1,
    "Supplementary Table S1. Locked EGFR-bypass component definitions",
    "TableS1_bypass_component_definitions.docx",
    footnotes = c(
        "The primary high-confidence exposure definition was locked before the first formal outcome analysis.",
        "Component counts may overlap because individual patients can carry more than one qualifying alteration.",
        "The strict core was literature anchored and PRESSING informed; it was not intended as an exact reproduction of the original PRESSING panel.",
        "Component-level effects are not interpreted as gene-specific efficacy estimates because individual carrier counts are sparse."
    ),
    orientation =
        "landscape",
    font_size =
        8.3,
    col_widths =
        c(
            2.25,
            3.15,
            1.25,
            4.0
        )
)


# ============================================================
# 7. TABLE S2
#    rwPFS robustness / sensitivity ONLY
# ============================================================

pretty_component <- c(
    ERBB2_AMP_R =
        "ERBB2 amplification",

    ERBB2_ACTIVATING_R =
        "ERBB2 activating mutation",

    MET_AMP_R =
        "MET amplification",

    PIK3CA_EX20_R =
        "PIK3CA exon 20",

    PTEN_LOF_R =
        "PTEN loss-of-function",

    KRAS_AMP_R =
        "KRAS amplification",

    MAP2K1_HOTSPOT_R =
        "MAP2K1 hotspot",

    NF1_LOF_R =
        "NF1 loss-of-function"
)

fmt_hr_ci <- function(
    hr,
    lo,
    hi
) {

    paste0(
        fmt_num(
            hr,
            2
        ),
        " (",
        fmt_num(
            lo,
            2
        ),
        "–",
        fmt_num(
            hi,
            2
        ),
        ")"
    )
}

# ------------------------------------------------------------
# 7.1 Primary
# ------------------------------------------------------------

primary_s2 <- as.data.table(
    b108$primary_result
)

s2_primary <- primary_s2[
    ,
    .(
        Section =
            "Primary",

        Analysis =
            "Primary adjusted rwPFS",

        N,

        Events =
            EVENTS,

        `Exposed, n` =
            BYPASS_POS,

        `HR (95% CI)` =
            fmt_hr_ci(
                HR,
                LCL95,
                UCL95
            ),

        `P value` =
            fmt_p(P)
    )
]


# ------------------------------------------------------------
# 7.2 Exposure-definition sensitivity
#     IMPORTANT: rwPFS rows ONLY
# ------------------------------------------------------------

exp_sens <- as.data.table(
    b107$exposure_sensitivity
)

if (
    !"OUTCOME" %in%
        names(exp_sens) ||
    !"EXPOSURE" %in%
        names(exp_sens)
) {
    stop(
        "B1-07 exposure_sensitivity structure mismatch."
    )
}

exp_sens <- exp_sens[
    OUTCOME == "rwPFS"
]

exp_count_map <- c(
    "PRESSING-like strict" =
        dt[
            BYPASS_PRESSING_STRICT_R == 1L,
            .N
        ],

    "High-confidence excluding MAP2K1 D67E" =
        dt[
            BYPASS_HIGH_CONFIDENCE_NO_D67E_R == 1L,
            .N
        ],

    "Original broad composite" =
        dt[
            BYPASS_ORIGINAL_BROAD_R == 1L,
            .N
        ]
)

exp_label_map <- c(
    "PRESSING-like strict" =
        "Strict literature-anchored composite",

    "High-confidence excluding MAP2K1 D67E" =
        "High-confidence composite excluding MAP2K1 D67E",

    "Original broad composite" =
        "Original broad composite"
)

exp_sens[
    ,
    ANALYSIS_LABEL :=
        unname(
            exp_label_map[
                EXPOSURE
            ]
        )
]

if (
    any(
        is.na(
            exp_sens$ANALYSIS_LABEL
        )
    )
) {
    stop(
        "Unknown exposure label found in B1-07 exposure sensitivity."
    )
}

exp_sens[
    ,
    EXPOSED_N :=
        as.integer(
            exp_count_map[
                EXPOSURE
            ]
        )
]

s2_exposure <- exp_sens[
    ,
    .(
        Section =
            "Exposure-definition sensitivity",

        Analysis =
            ANALYSIS_LABEL,

        N,

        Events =
            EVENTS,

        `Exposed, n` =
            EXPOSED_N,

        `HR (95% CI)` =
            fmt_hr_ci(
                HR,
                LCL95,
                UCL95
            ),

        `P value` =
            fmt_p(P)
    )
]


# ------------------------------------------------------------
# 7.3 Cohort-restriction sensitivity
#     IMPORTANT: rwPFS rows ONLY
# ------------------------------------------------------------

restriction <- as.data.table(
    b107$restriction_sensitivity
)

if (
    !"OUTCOME" %in%
        names(restriction) ||
    !"ANALYSIS" %in%
        names(restriction)
) {
    stop(
        "B1-07 restriction_sensitivity structure mismatch."
    )
}

restriction <- restriction[
    OUTCOME == "rwPFS"
]

restriction[
    ,
    EXPOSED_N :=
        fcase(

            ANALYSIS ==
                "Restricted: exclude BRAF class III",

            dt[
                BRAF_CLASS3_R == 0L &
                    BYPASS_HIGH_CONFIDENCE_R == 1L,
                .N
            ],

            ANALYSIS ==
                "Restricted: specimen <=730 days",

            dt[
                SPECIMEN_WITHIN_730D_R == 1L &
                    BYPASS_HIGH_CONFIDENCE_R == 1L,
                .N
            ],

            ANALYSIS ==
                "Restricted: exclude sidedness conflict",

            dt[
                SIDE_CONFLICT_R == 0L &
                    BYPASS_HIGH_CONFIDENCE_R == 1L,
                .N
            ],

            ANALYSIS ==
                "Restricted: exclude BRAF class III + specimen <=730 days",

            dt[
                BRAF_CLASS3_R == 0L &
                    SPECIMEN_WITHIN_730D_R == 1L &
                    BYPASS_HIGH_CONFIDENCE_R == 1L,
                .N
            ],

            default =
                NA_integer_
        )
]

restriction[
    ,
    ANALYSIS_LABEL :=
        fcase(

            ANALYSIS ==
                "Restricted: exclude BRAF class III",

            "Exclude BRAF class III",

            ANALYSIS ==
                "Restricted: specimen <=730 days",

            "Specimen acquired within 730 days before anti-EGFR",

            ANALYSIS ==
                "Restricted: exclude sidedness conflict",

            "Exclude sidedness-conflict case",

            ANALYSIS ==
                "Restricted: exclude BRAF class III + specimen <=730 days",

            "Exclude BRAF class III and restrict specimen to <=730 days",

            default =
                ANALYSIS
        )
]

s2_restriction <- restriction[
    ,
    .(
        Section =
            "Cohort-restriction sensitivity",

        Analysis =
            ANALYSIS_LABEL,

        N,

        Events =
            EVENTS,

        `Exposed, n` =
            EXPOSED_N,

        `HR (95% CI)` =
            fmt_hr_ci(
                HR,
                LCL95,
                UCL95
            ),

        `P value` =
            fmt_p(P)
    )
]


# ------------------------------------------------------------
# 7.4 Leave-one-component-out composite
# ------------------------------------------------------------

loo <- as.data.table(
    b108$leave_one_component_out
)

if (
    !"COMPONENT" %in%
        names(loo)
) {
    stop(
        "B1-08 leave_one_component_out missing COMPONENT."
    )
}

loo[
    ,
    COMPONENT_LABEL :=
        unname(
            pretty_component[
                COMPONENT
            ]
        )
]

if (
    any(
        is.na(
            loo$COMPONENT_LABEL
        )
    )
) {
    stop(
        "Unknown component in leave-one-component-out results."
    )
}

s2_loo <- loo[
    ,
    .(
        Section =
            "Leave-one-component-out",

        Analysis =
            paste0(
                "Leave out ",
                COMPONENT_LABEL
            ),

        N,

        Events =
            EVENTS,

        `Exposed, n` =
            BYPASS_POS,

        `HR (95% CI)` =
            fmt_hr_ci(
                HR,
                LCL95,
                UCL95
            ),

        `P value` =
            fmt_p(P)
    )
]


# ------------------------------------------------------------
# 7.5 Exclude carriers of each component
# ------------------------------------------------------------

carrier_excl <- as.data.table(
    b108$carrier_exclusion
)

if (
    !"COMPONENT" %in%
        names(carrier_excl)
) {
    stop(
        "B1-08 carrier_exclusion missing COMPONENT."
    )
}

carrier_excl[
    ,
    COMPONENT_LABEL :=
        unname(
            pretty_component[
                COMPONENT
            ]
        )
]

if (
    any(
        is.na(
            carrier_excl$COMPONENT_LABEL
        )
    )
) {
    stop(
        "Unknown component in carrier-exclusion results."
    )
}

s2_carrier <- carrier_excl[
    ,
    .(
        Section =
            "Carrier-exclusion sensitivity",

        Analysis =
            paste0(
                "Exclude carriers of ",
                COMPONENT_LABEL
            ),

        N,

        Events =
            EVENTS,

        `Exposed, n` =
            BYPASS_POS,

        `HR (95% CI)` =
            fmt_hr_ci(
                HR,
                LCL95,
                UCL95
            ),

        `P value` =
            fmt_p(P)
    )
]


# ------------------------------------------------------------
# 7.6 Landmark sensitivity
# ------------------------------------------------------------

landmark <- as.data.table(
    b108$landmark
)

s2_landmark <- landmark[
    ,
    .(
        Section =
            "Landmark sensitivity",

        Analysis =
            ANALYSIS,

        N,

        Events =
            EVENTS,

        `Exposed, n` =
            BYPASS_POS,

        `HR (95% CI)` =
            fmt_hr_ci(
                HR,
                LCL95,
                UCL95
            ),

        `P value` =
            fmt_p(P)
    )
]


# ------------------------------------------------------------
# 7.7 Combine
# ------------------------------------------------------------

tableS2 <- rbindlist(
    list(
        s2_primary,
        s2_exposure,
        s2_restriction,
        s2_loo,
        s2_carrier,
        s2_landmark
    ),
    fill = TRUE
)

# Hard audit: S2 must contain rwPFS-event counts only.
# No OS row with 129 events may survive.
if (
    any(
        tableS2$Events == 129L,
        na.rm = TRUE
    )
) {
    stop(
        "S2 hard audit failed: an OS row (129 events) entered the rwPFS table."
    )
}

# Required exposure-definition anchors.
required_exp_rows <- c(
    "Strict literature-anchored composite",
    "High-confidence composite excluding MAP2K1 D67E",
    "Original broad composite"
)

if (
    !all(
        required_exp_rows %in%
            tableS2$Analysis
    )
) {
    stop(
        "S2 hard audit failed: exposure-definition labels are incomplete."
    )
}

# Required leave-one-component-out labels.
required_loo <- paste0(
    "Leave out ",
    c(
        "ERBB2 amplification",
        "ERBB2 activating mutation",
        "MET amplification",
        "PIK3CA exon 20",
        "PTEN loss-of-function",
        "KRAS amplification",
        "MAP2K1 hotspot",
        "NF1 loss-of-function"
    )
)

if (
    !all(
        required_loo %in%
            tableS2$Analysis
    )
) {
    stop(
        "S2 hard audit failed: leave-one-component-out labels are incomplete."
    )
}

boot_hr <- as.numeric(
    b108$bootstrap_hr
)

boot_hr <- boot_hr[
    is.finite(
        boot_hr
    )
]

bootstrap_note <- paste0(
    "Patient-level bootstrap: ",
    length(
        boot_hr
    ),
    "/2,000 successful iterations; median HR ",
    sprintf(
        "%.2f",
        median(
            boot_hr
        )
    ),
    "; empirical 95% interval ",
    sprintf(
        "%.2f",
        quantile(
            boot_hr,
            0.025
        )
    ),
    "–",
    sprintf(
        "%.2f",
        quantile(
            boot_hr,
            0.975
        )
    ),
    "; ",
    sprintf(
        "%.1f%%",
        100 *
            mean(
                boot_hr > 1
            )
    ),
    " of bootstrap HRs were >1."
)

write_table_docx(
    tableS2,
    "Supplementary Table S2. Robustness and sensitivity analyses for rwPFS",
    "TableS2_robustness_results.docx",
    footnotes = c(
        "All rows refer to rwPFS. No OS results are included in this table.",
        "Exposed, n denotes the number classified as bypass-positive under the exposure definition or restriction used in the corresponding row.",
        "Leave-one-component-out analyses reconstructed the composite after removing the specified alteration while retaining all other qualifying components.",
        "Carrier-exclusion analyses removed all carriers of the specified alteration and then re-estimated the association using the locked high-confidence composite among the remaining patients.",
        "Landmark analyses were post-unblinding data-quality sensitivity analyses and did not replace the primary analysis.",
        bootstrap_note
    ),
    orientation =
        "landscape",
    font_size =
        7.6,
    col_widths =
        c(
            1.85,
            3.75,
            0.50,
            0.60,
            0.75,
            1.45,
            0.70
        )
)


# ============================================================
# 8. TABLE S3
# ============================================================

s3 <- as.data.table(
    b110$model_results
)

s3[
    ,
    Endpoint :=
        fcase(
            ENDPOINT ==
                "Conservative CRC TTNTD30",
            "CRC-compatible TTNTD (30-day grace)",

            ENDPOINT ==
                "Conservative CRC TTNTD60",
            "CRC-compatible TTNTD (60-day grace)",

            ENDPOINT ==
                "Broad cancer-directed TTNTD30",
            "Broad cancer-directed TTNTD (30-day grace)",

            ENDPOINT ==
                "Broad cancer-directed TTNTD60",
            "Broad cancer-directed TTNTD (60-day grace)",

            default =
                ENDPOINT
        )
]

s3[
    ,
    Model :=
        fcase(
            MODEL ==
                "Standard adjusted",
            "Adjusted Cox",

            MODEL ==
                "PH-addressed stratified",
            "Adjusted Cox, stratified",

            default =
                MODEL
        )
]

tableS3 <- s3[
    ,
    .(
        Endpoint,
        Model,
        N,
        Events = EVENTS,

        `HR (95% CI)` =
            paste0(
                fmt_num(
                    HR,
                    2
                ),
                " (",
                fmt_num(
                    LCL95,
                    2
                ),
                "–",
                fmt_num(
                    UCL95,
                    2
                ),
                ")"
            ),

        `P value` =
            fmt_p(P)
    )
]

write_table_docx(
    tableS3,
    "Supplementary Table S3. Refined treatment-based TTNTD validation",
    "TableS3_TTNTD_validation.docx",
    footnotes = c(
        "TTNTD denotes time to next treatment or death.",
        "CRC-compatible endpoints counted only initiation of a new CRC-compatible systemic/targeted treatment or death after the regimen-assembly grace period.",
        "Broad cancer-directed endpoints excluded prespecified clearly non-CRC/supportive/endocrine treatments but otherwise allowed new cancer-directed therapy or death.",
        "Adjusted Cox models used the locked primary covariate set; stratified models stratified by primary subsite and prior heavy cytotoxic exposure."
    ),
    orientation =
        "portrait",
    font_size =
        8.3,
    col_widths =
        c(
            2.65,
            1.85,
            0.55,
            0.65,
            1.55,
            0.75
        )
)


# ============================================================
# 9. Audit
# ============================================================

expected <- file.path(
    table_dir,
    c(
        "Table1_baseline_characteristics.docx",
        "Table2_primary_secondary_results.docx",
        "TableS1_bypass_component_definitions.docx",
        "TableS2_robustness_results.docx",
        "TableS3_TTNTD_validation.docx"
    )
)

missing <- expected[
    !file.exists(
        expected
    )
]

audit_file <- file.path(
    audit_dir,
    "B1_12_final_manuscript_tables_audit.txt"
)

writeLines(
    c(
        "B1-12 FINAL MANUSCRIPT TABLES AUDIT",
        "===================================",
        "",
        paste0(
            "Locked cohort N = ",
            nrow(dt)
        ),
        paste0(
            "Bypass- = ",
            dt[
                BYPASS_HIGH_CONFIDENCE_R == 0L,
                .N
            ]
        ),
        paste0(
            "Bypass+ = ",
            dt[
                BYPASS_HIGH_CONFIDENCE_R == 1L,
                .N
            ]
        ),
        "",
        "Table titles use unnumbered manuscript-facing paragraphs.",
        "Table S1 uses landscape orientation to prevent clipping.",
        "Internal R-code labels were removed from manuscript-facing tables.",
        "Table S2 contains rwPFS rows only; OS rows are explicitly excluded.",
        "Table S2 exposure/component labels are manuscript-facing and non-generic.",
        "No new primary analysis was introduced.",
        "",
        paste0(
            "Missing outputs = ",
            length(missing)
        ),
        if (
            length(missing) > 0L
        ) {
            missing
        } else {
            "All expected table outputs exist."
        }
    ),
    audit_file
)

cat("\n============================================================\n")
cat("B1-12 COMPLETE\n")
cat("============================================================\n")

cat(
    "Missing outputs: ",
    length(missing),
    "\n",
    sep = ""
)

cat(
    "\nAudit:\n",
    audit_file,
    "\n"
)
