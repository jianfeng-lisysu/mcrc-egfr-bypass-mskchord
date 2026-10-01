# ============================================================
# A7 / B-line
# 11B_B1_11_generate_manuscript_outputs_colored.R
#
# PURPOSE
#   Generate the first publication-ready manuscript package
#   from the LOCKED B-line analyses.
#
# NO NEW PRIMARY ANALYSIS IS INTRODUCED HERE.
#
# INPUTS
#   03_intermediate/A7_B1_locked_analysis.rds
#   04_results/B1_07_primary_models.rds
#   04_results/B1_08_robustness_validation.rds
#   04_results/B1_10_refined_TTNTD_PH_validation.rds
#
# MAIN OUTPUTS
#   05_figures/Figure1_study_design_and_cohort_flow.png/tiff
#   05_figures/Figure2_high_confidence_bypass_landscape.png/tiff
#   05_figures/Figure3_primary_rwPFS_OS.png/tiff
#   05_figures/Figure4_robustness_validation.png/tiff
#   05_figures/Figure5_TTNTD_validation.png/tiff
#
#   07_tables/Table1_baseline_characteristics.docx
#   07_tables/Table2_primary_secondary_results.docx
#
# SUPPLEMENTARY TABLES
#   07_tables/TableS1_bypass_component_definitions.docx
#   07_tables/TableS2_robustness_results.docx
#   07_tables/TableS3_TTNTD_validation.docx
#
# AUDIT
#   06_logs_and_audit/B1_11_manuscript_outputs_audit.txt
#
# FIGURE RULES
#   - overwrite same file names
#   - PNG + TIFF
#   - TIFF 600 dpi
#   - no subfolders in 05_figures
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory

intermediate_dir <- file.path(project_root, "03_intermediate")
results_dir <- file.path(project_root, "04_results")
fig_dir <- file.path(project_root, "05_figures")
audit_dir <- file.path(project_root, "06_logs_and_audit")
table_dir <- file.path(project_root, "07_tables")

for (d in c(fig_dir, audit_dir, table_dir)) {
    if (!dir.exists(d)) {
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
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
        stop("Missing required input:\n", f)
    }
}

packages <- c(
    "data.table",
    "survival",
    "ggplot2",
    "patchwork",
    "scales",
    "officer",
    "flextable"
)

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

library(data.table)
library(survival)
library(ggplot2)
library(patchwork)
library(scales)
library(officer)
library(flextable)

# ============================================================
# Publication palette
# Okabe-Ito-inspired; consistent semantic colors across figures.
# ============================================================

COL_NEG <- "#0072B2"       # bypass-negative: blue
COL_POS <- "#D55E00"       # bypass-positive: vermillion
COL_TEAL <- "#009E73"      # secondary validation / follow-up
COL_ORANGE <- "#E69F00"    # alternate treatment definition
COL_SKY <- "#56B4E9"
COL_PURPLE <- "#CC79A7"
COL_DARK <- "#333333"
COL_LIGHT <- "#E9EEF2"
COL_LIGHT_BLUE <- "#DCEEF8"
COL_LIGHT_ORANGE <- "#FCE8D5"
COL_LIGHT_GREEN <- "#DDF2E9"

cat("\n============================================================\n")
cat("B1-11B COLORED MANUSCRIPT OUTPUTS - STALE-FUNCTION-SAFE\n")
cat("============================================================\n\n")


# ============================================================
# 1. Read locked results
# ============================================================

dt <- as.data.table(readRDS(locked_file))
b107 <- readRDS(b107_file)
b108 <- readRDS(b108_file)
b110 <- readRDS(b110_file)

if (
    nrow(dt) != 191L ||
    dt[BYPASS_HIGH_CONFIDENCE_R == 1L, .N] != 23L
) {
    stop("Locked cohort anchor failed.")
}

if (!all(c(
    "primary_results",
    "exposure_sensitivity",
    "restriction_sensitivity",
    "km_rwpfs",
    "km_os"
) %in% names(b107))) {
    stop("B1-07 model object structure mismatch.")
}

if (!all(c(
    "primary_result",
    "bootstrap_summary",
    "component_counts",
    "leave_one_component_out",
    "carrier_exclusion",
    "landmark",
    "rmst"
) %in% names(b108))) {
    stop("B1-08 object structure mismatch.")
}

if (!all(c(
    "patient_data",
    "model_results",
    "event_summary"
) %in% names(b110))) {
    stop("B1-10 object structure mismatch.")
}

ttntd_dt <- as.data.table(b110$patient_data)


# ============================================================
# 2. General helpers
# ============================================================

safe_chr <- function(x) {
    z <- as.character(x)
    z[is.na(z)] <- ""
    z
}

b111_fmt_num <- function(x, digits = 2) {
    vapply(
        x,
        function(z) {
            if (length(z) == 0L || is.na(z) || !is.finite(z)) {
                return("NA")
            }
            formatC(z, format = "f", digits = digits)
        },
        character(1)
    )
}

b111_fmt_p <- function(x) {
    vapply(
        x,
        function(z) {
            if (length(z) == 0L || is.na(z) || !is.finite(z)) {
                return("NA")
            }
            if (z < 0.001) {
                return("<0.001")
            }
            formatC(z, format = "f", digits = 3)
        },
        character(1)
    )
}

median_iqr_text <- function(x, digits = 1) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
        return("NA")
    }

    q <- quantile(
        x,
        c(0.25, 0.50, 0.75),
        na.rm = TRUE,
        names = FALSE
    )

    paste0(
        b111_fmt_num(q[2], digits),
        " [",
        b111_fmt_num(q[1], digits),
        ", ",
        b111_fmt_num(q[3], digits),
        "]"
    )
}

n_pct_text <- function(n, d, digits = 1) {
    paste0(
        n,
        " (",
        b111_fmt_num(100 * n / d, digits),
        "%)"
    )
}

theme_pub <- function(base_size = 10) {
    theme_bw(base_size = base_size) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey95"),
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold", size = base_size + 1),
            plot.subtitle = element_text(size = base_size),
            axis.title = element_text(face = "bold"),
            legend.title = element_text(face = "bold"),
            legend.position = "top",
            plot.margin = margin(6, 6, 6, 6)
        )
}

save_pub_figure <- function(
    plot,
    basename,
    width,
    height
) {
    png_file <- file.path(
        fig_dir,
        paste0(basename, ".png")
    )

    tiff_file <- file.path(
        fig_dir,
        paste0(basename, ".tiff")
    )

    ggsave(
        png_file,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = 300,
        bg = "white"
    )

    ggsave(
        tiff_file,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = 600,
        compression = "lzw",
        bg = "white"
    )

    invisible(
        c(
            png_file,
            tiff_file
        )
    )
}

three_line_ft <- function(df) {
    ft <- flextable(df)

    ft <- theme_booktabs(ft)

    ft <- border_remove(ft)

    border_black <- fp_border(
        color = "black",
        width = 1
    )

    border_mid <- fp_border(
        color = "black",
        width = 0.75
    )

    ft <- hline_top(
        ft,
        border = border_black,
        part = "header"
    )

    ft <- hline_bottom(
        ft,
        border = border_mid,
        part = "header"
    )

    ft <- hline_bottom(
        ft,
        border = border_black,
        part = "body"
    )

    ft <- fontsize(
        ft,
        size = 9,
        part = "all"
    )

    ft <- font(
        ft,
        fontname = "Arial",
        part = "all"
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

    ft <- bold(
        ft,
        part = "header"
    )

    ft <- autofit(ft)

    ft
}

write_word_table <- function(
    df,
    title,
    filename,
    footnotes = NULL
) {
    doc <- read_docx()

    doc <- body_add_par(
        doc,
        title,
        style = "heading 1"
    )

    ft <- three_line_ft(df)

    doc <- body_add_flextable(
        doc,
        ft
    )

    if (!is.null(footnotes)) {
        for (note in footnotes) {
            doc <- body_add_par(
                doc,
                note,
                style = "Normal"
            )
        }
    }

    print(
        doc,
        target = file.path(
            table_dir,
            filename
        )
    )
}

survfit_to_df <- function(
    sf,
    group_label
) {
    s <- summary(sf)

    data.table(
        time = s$time,
        surv = s$surv,
        lower = s$lower,
        upper = s$upper,
        strata = as.character(s$strata)
    )[
        ,
        group := fifelse(
            grepl("=1$", strata),
            "High-confidence bypass+",
            "High-confidence bypass−"
        )
    ][
        ,
        group := factor(
            group,
            levels = c(
                "High-confidence bypass−",
                "High-confidence bypass+"
            )
        )
    ]
}

km_plot <- function(
    sf,
    title,
    xlab,
    x_limit = NULL
) {
    d <- survfit_to_df(sf, "group")

    p <- ggplot(
        d,
        aes(
            x = time,
            y = surv,
            linetype = group,
            color = group
        )
    ) +
        geom_step(
            linewidth = 0.95
        ) +
        scale_linetype_manual(
            values = c(
                "High-confidence bypass−" = "solid",
                "High-confidence bypass+" = "dashed"
            )
        ) +
        scale_color_manual(
            values = c(
                "High-confidence bypass−" = COL_NEG,
                "High-confidence bypass+" = COL_POS
            )
        ) +
        scale_y_continuous(
            labels = percent_format(
                accuracy = 1
            ),
            limits = c(0, 1)
        ) +
        labs(
            title = title,
            x = xlab,
            y = "Survival probability",
            linetype = NULL,
            color = NULL
        ) +
        theme_pub(10)

    if (!is.null(x_limit)) {
        p <- p +
            coord_cartesian(
                xlim = c(0, x_limit)
            )
    }

    p
}

forest_plot <- function(
    d,
    title,
    label_col = "LABEL"
) {
    d <- copy(d)

    d[
        ,
        LABEL := factor(
            get(label_col),
            levels = rev(
                unique(
                    get(label_col)
                )
            )
        )
    ]

    ggplot(
        d,
        aes(
            x = HR,
            y = LABEL
        )
    ) +
        geom_vline(
            xintercept = 1,
            linetype = "dotted"
        ) +
        geom_errorbarh(
            aes(
                xmin = LCL95,
                xmax = UCL95
            ),
            height = 0.15,
            color = COL_POS,
            linewidth = 0.7
        ) +
        geom_point(
            size = 2.5,
            color = COL_POS
        ) +
        scale_x_log10(
            breaks = c(
                0.5,
                1,
                2,
                4,
                8
            )
        ) +
        labs(
            title = title,
            x = "Hazard ratio (log scale)",
            y = NULL
        ) +
        theme_pub(9) +
        theme(
            legend.position = "none"
        )
}


# ============================================================
# 3. FIGURE 1: study design + cohort flow
# ============================================================

flow_nodes <- data.table(
    x = 1,
    y = c(5, 4, 3, 2, 1),
    label = c(
        "MSK-CHORD colorectal cancer\nn = 5,543",
        "Legacy genomic eligibility screen\nleft-sided stage IV, MSS,\nno detected KRAS/NRAS mutation,\nBRAF V600E-negative\nn = 803",
        "First relevant biologic = anti-EGFR\nn = 230",
        "Sequencing completed on/before\nfirst anti-EGFR exposure\nn = 191",
        "Frozen analysis cohort\nhigh-confidence bypass+ n = 23\nbypass− n = 168"
    ),
    fill_color = c(
        COL_LIGHT,
        COL_LIGHT_BLUE,
        COL_LIGHT_ORANGE,
        COL_LIGHT_BLUE,
        COL_LIGHT_GREEN
    ),
    border_color = c(
        COL_DARK,
        COL_NEG,
        COL_POS,
        COL_NEG,
        COL_TEAL
    )
)

p1a <- ggplot() +
    geom_segment(
        data = data.table(
            x = 1,
            xend = 1,
            y = c(4.65, 3.65, 2.65, 1.65),
            yend = c(4.35, 3.35, 2.35, 1.35)
        ),
        aes(
            x = x,
            xend = xend,
            y = y,
            yend = yend
        ),
        arrow = arrow(
            length = unit(
                0.08,
                "inches"
            )
        ),
        linewidth = 0.6
    ) +
    geom_label(
        data = flow_nodes,
        aes(
            x = x,
            y = y,
            label = label,
            fill = fill_color,
            color = border_color
        ),
        size = 3.0,
        label.size = 0.5,
        label.padding = unit(
            0.18,
            "lines"
        )
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    coord_cartesian(
        xlim = c(0.25, 1.75),
        ylim = c(0.55, 5.45)
    ) +
    labs(
        title = "A  Cohort derivation"
    ) +
    theme_void(base_size = 10) +
    theme(
        plot.title = element_text(
            face = "bold"
        )
    )

timeline <- data.table(
    x = c(
        0,
        1.5,
        3,
        4.5,
        6
    ),
    label = c(
        "Stage IV\nregistry diagnosis",
        "Pre-treatment\ntumor specimen",
        "MSK-IMPACT\nsequencing",
        "First anti-EGFR\n(index)",
        "rwPFS / OS /\nTTNTD follow-up"
    ),
    point_color = c(
        COL_DARK,
        COL_NEG,
        COL_NEG,
        COL_POS,
        COL_TEAL
    )
)

p1b <- ggplot() +
    geom_segment(
        aes(
            x = 0,
            xend = 6,
            y = 1,
            yend = 1
        ),
        linewidth = 0.8
    ) +
    geom_point(
        data = timeline,
        aes(
            x = x,
            y = 1,
            color = point_color
        ),
        size = 3.0
    ) +
    scale_color_identity() +
    geom_text(
        data = timeline,
        aes(
            x = x,
            y = 0.83,
            label = label
        ),
        size = 3,
        vjust = 1
    ) +
    annotate(
        "rect",
        xmin = 0,
        xmax = 4.5,
        ymin = 1.12,
        ymax = 1.33,
        fill = COL_LIGHT_BLUE,
        color = COL_NEG,
        linewidth = 0.35
    ) +
    annotate(
        "text",
        x = 2.25,
        y = 1.225,
        label = "Genomic ascertainment must precede anti-EGFR index",
        size = 3.1
    ) +
    annotate(
        "text",
        x = 3.0,
        y = 1.57,
        label = "Frozen exposure: high-confidence EGFR-bypass alterations",
        size = 3.25,
        fontface = "bold"
    ) +
    coord_cartesian(
        xlim = c(-0.4, 6.4),
        ylim = c(0.05, 1.8)
    ) +
    labs(
        title = "B  Study time zero and endpoint framework"
    ) +
    theme_void(base_size = 10) +
    theme(
        plot.title = element_text(
            face = "bold"
        )
    )

fig1 <- p1a / p1b +
    plot_layout(
        heights = c(1.45, 0.75)
    )

save_pub_figure(
    fig1,
    "Figure1_study_design_and_cohort_flow",
    width = 8.2,
    height = 9.2
)


# ============================================================
# 4. FIGURE 2: high-confidence bypass genomic landscape
# ============================================================

component_vars <- c(
    "ERBB2_AMP_R",
    "ERBB2_ACTIVATING_R",
    "MET_AMP_R",
    "PIK3CA_EX20_R",
    "PTEN_LOF_R",
    "KRAS_AMP_R",
    "MAP2K1_HOTSPOT_R",
    "NF1_LOF_R"
)

component_labels <- c(
    "ERBB2 amplification",
    "ERBB2 activating mutation",
    "MET amplification",
    "PIK3CA exon 20",
    "PTEN loss-of-function",
    "KRAS amplification",
    "MAP2K1 hotspot",
    "NF1 loss-of-function"
)

names(component_labels) <- component_vars

comp_counts <- data.table(
    component = component_vars,
    label = component_labels
)[
    ,
    n := vapply(
        component,
        function(v) {
            sum(
                dt[[v]] == 1L,
                na.rm = TRUE
            )
        },
        integer(1)
    )
][
    ,
    prevalence := n / nrow(dt)
]

comp_counts[
    ,
    label := factor(
        label,
        levels = rev(
            label
        )
    )
]

p2a <- ggplot(
    comp_counts,
    aes(
        x = prevalence,
        y = label
    )
) +
    geom_col(
        width = 0.65,
        fill = COL_NEG
    ) +
    geom_text(
        aes(
            label = paste0(
                n,
                " (",
                percent(
                    prevalence,
                    accuracy = 0.1
                ),
                ")"
            )
        ),
        hjust = -0.08,
        size = 3
    ) +
    scale_x_continuous(
        labels = percent_format(
            accuracy = 1
        ),
        expand = expansion(
            mult = c(0, 0.25)
        )
    ) +
    labs(
        title = "A  Prevalence of high-confidence bypass components",
        x = "Prevalence in the locked cohort",
        y = NULL
    ) +
    theme_pub(9) +
    theme(
        legend.position = "none"
    )

carriers <- dt[
    BYPASS_HIGH_CONFIDENCE_R == 1L
]

carrier_scores <- rowSums(
    as.data.frame(
        carriers[
            ,
            ..component_vars
        ]
    ),
    na.rm = TRUE
)

carriers[
    ,
    CARRIER_SCORE_R := carrier_scores
]

setorder(
    carriers,
    -CARRIER_SCORE_R,
    -ERBB2_AMP_R,
    -ERBB2_ACTIVATING_R,
    -PIK3CA_EX20_R,
    -KRAS_AMP_R
)

carriers[
    ,
    patient_order := seq_len(.N)
]

heat <- rbindlist(
    lapply(
        component_vars,
        function(v) {
            data.table(
                patient_order =
                    carriers$patient_order,
                patient =
                    paste0(
                        "P",
                        sprintf(
                            "%02d",
                            carriers$patient_order
                        )
                    ),
                component =
                    component_labels[[v]],
                present =
                    carriers[[v]]
            )
        }
    )
)

heat[
    ,
    patient := factor(
        patient,
        levels = paste0(
            "P",
            sprintf(
                "%02d",
                seq_len(
                    nrow(carriers)
                )
            )
        )
    )
]

heat[
    ,
    component := factor(
        component,
        levels = rev(
            component_labels
        )
    )
]

p2b <- ggplot(
    heat,
    aes(
        x = patient,
        y = component,
        fill = factor(
            present
        )
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.25
    ) +
    scale_fill_manual(
        values = c(
            "0" = COL_LIGHT,
            "1" = COL_POS
        ),
        labels = c(
            "Absent",
            "Present"
        )
    ) +
    labs(
        title = "B  Component landscape among 23 high-confidence bypass+ patients",
        x = "High-confidence bypass+ patients",
        y = NULL,
        fill = NULL
    ) +
    theme_pub(8.5) +
    theme(
        axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1,
            size = 6.5
        ),
        legend.position = "top"
    )

fig2 <- p2a / p2b +
    plot_layout(
        heights = c(0.85, 1.35)
    )

save_pub_figure(
    fig2,
    "Figure2_high_confidence_bypass_landscape",
    width = 8.6,
    height = 8.6
)


# ============================================================
# 5. FIGURE 3: primary rwPFS and secondary OS
# ============================================================

p3a <- km_plot(
    b107$km_rwpfs,
    "A  Real-world progression-free survival",
    "Days from first anti-EGFR exposure",
    x_limit = 800
)

p3b <- km_plot(
    b107$km_os,
    "B  Overall survival",
    "Days from first anti-EGFR exposure",
    x_limit = 1400
)

primary_results <- as.data.table(
    b107$primary_results
)

forest_main <- primary_results[
    ANALYSIS %chin%
        c(
            "Crude",
            "Primary adjusted"
        )
][
    ,
    LABEL := paste0(
        OUTCOME,
        " — ",
        ANALYSIS
    )
]

p3c <- forest_plot(
    forest_main,
    "C  Hazard-ratio estimates",
    "LABEL"
)

rmst <- as.data.table(
    b108$rmst
)

rmst_long <- rbindlist(
    list(
        rmst[
            ,
            .(
                TAU_DAYS,
                group =
                    "High-confidence bypass−",
                rmst =
                    RMST_BYPASS_NEG_DAYS
            )
        ],
        rmst[
            ,
            .(
                TAU_DAYS,
                group =
                    "High-confidence bypass+",
                rmst =
                    RMST_BYPASS_POS_DAYS
            )
        ]
    )
)

rmst_long[
    ,
    rmst_label := sprintf("%.1f", rmst)
]

rmst_long[
    ,
    TAU_LABEL := factor(
        paste0(
            TAU_DAYS,
            "-day horizon"
        ),
        levels = c(
            "180-day horizon",
            "365-day horizon"
        )
    )
]

p3d <- ggplot(
    rmst_long,
    aes(
        x = TAU_LABEL,
        y = rmst,
        fill = group
    )
) +
    geom_col(
        position = position_dodge(
            width = 0.75
        ),
        width = 0.65
    ) +
    geom_text(
        aes(label = rmst_label),
        position = position_dodge(
            width = 0.75
        ),
        vjust = -0.3,
        size = 3
    ) +
    scale_fill_manual(
        values = c(
            "High-confidence bypass−" = COL_NEG,
            "High-confidence bypass+" = COL_POS
        )
    ) +
    labs(
        title = "D  Restricted mean rwPFS",
        x = NULL,
        y = "RMST, days",
        fill = NULL
    ) +
    theme_pub(9) +
    theme(
        legend.position = "top"
    )

fig3 <- (p3a | p3b) /
    (p3c | p3d) +
    plot_layout(
        heights = c(1.0, 0.95)
    )

save_pub_figure(
    fig3,
    "Figure3_primary_rwPFS_OS",
    width = 11.2,
    height = 8.4
)


# ============================================================
# 6. FIGURE 4: robustness
# ============================================================

robust <- rbindlist(
    list(
        as.data.table(
            b108$primary_result
        )[
            ,
            LABEL :=
                "Primary adjusted"
        ],

        as.data.table(
            b108$leave_one_component_out
        )[
            ,
            LABEL := paste0(
                "Leave out ",
                gsub(
                    "_R$",
                    "",
                    COMPONENT
                )
            )
        ],

        as.data.table(
            b108$landmark
        )[
            ,
            LABEL := ANALYSIS
        ]
    ),
    fill = TRUE
)

robust <- robust[
    is.finite(HR) &
        is.finite(LCL95) &
        is.finite(UCL95)
]

p4a <- forest_plot(
    robust,
    "A  Robustness of adjusted rwPFS association",
    "LABEL"
)

boot_values <- as.numeric(
    b108$bootstrap_hr
)

boot_df <- data.table(
    HR = boot_values
)

primary_hr <- as.data.table(
    b108$primary_result
)$HR[1]

p4b <- ggplot(
    boot_df,
    aes(
        x = HR
    )
) +
    geom_histogram(
        bins = 35,
        boundary = 1,
        fill = COL_NEG,
        color = "white"
    ) +
    geom_vline(
        xintercept = 1,
        linetype = "dotted",
        linewidth = 0.7
    ) +
    geom_vline(
        xintercept = primary_hr,
        linetype = "dashed",
        linewidth = 0.9,
        color = COL_POS
    ) +
    labs(
        title = "B  Patient-level bootstrap distribution",
        subtitle = paste0(
            "2,000/2,000 successful; median HR ",
            b111_fmt_num(
                median(
                    boot_values
                ),
                2
            )
        ),
        x = "Adjusted rwPFS hazard ratio",
        y = "Bootstrap iterations"
    ) +
    theme_pub(9)

component_counts <- as.data.table(
    b108$component_counts
)

component_counts[
    ,
    COMPONENT_LABEL := gsub(
        "_R$",
        "",
        COMPONENT
    )
]

component_counts[
    ,
    COMPONENT_LABEL := factor(
        COMPONENT_LABEL,
        levels = rev(
            COMPONENT_LABEL
        )
    )
]

p4c <- ggplot(
    component_counts,
    aes(
        x = N_CARRIERS,
        y = COMPONENT_LABEL
    )
) +
    geom_col(
        width = 0.65,
        fill = COL_NEG
    ) +
    geom_text(
        aes(
            label = N_CARRIERS
        ),
        hjust = -0.15,
        size = 3
    ) +
    scale_x_continuous(
        expand = expansion(
            mult = c(0, 0.2)
        )
    ) +
    labs(
        title = "C  Carrier counts",
        x = "Patients",
        y = NULL
    ) +
    theme_pub(9)

rmst_diff <- as.data.table(
    b108$rmst
)

rmst_diff[
    ,
    diff_label := sprintf("%.1f", DIFFERENCE_POS_MINUS_NEG_DAYS)
]

rmst_diff[
    ,
    HORIZON := factor(
        paste0(
            TAU_DAYS,
            " days"
        ),
        levels = c(
            "180 days",
            "365 days"
        )
    )
]

p4d <- ggplot(
    rmst_diff,
    aes(
        x = HORIZON,
        y = DIFFERENCE_POS_MINUS_NEG_DAYS
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = "dotted"
    ) +
    geom_col(
        width = 0.55,
        fill = COL_POS
    ) +
    geom_text(
        aes(label = diff_label),
        vjust = 1.4,
        color = "white",
        size = 3.2
    ) +
    labs(
        title = "D  Absolute rwPFS deficit",
        x = "RMST horizon",
        y = "Bypass+ minus bypass−, days"
    ) +
    theme_pub(9)

fig4 <- (p4a | p4b) /
    (p4c | p4d)

save_pub_figure(
    fig4,
    "Figure4_robustness_validation",
    width = 11.3,
    height = 9.0
)


# ============================================================
# 7. FIGURE 5: refined TTNTD validation
# ============================================================

tt_models <- as.data.table(
    b110$model_results
)

tt_models[
    ,
    LABEL := paste0(
        gsub(
            "Conservative CRC ",
            "Conservative ",
            ENDPOINT
        ),
        " — ",
        MODEL
    )
]

p5a <- forest_plot(
    tt_models,
    "A  Refined treatment-based validation",
    "LABEL"
)

# Build KM curves for conservative CRC TTNTD30 from patient data.
tt_sf <- survfit(
    Surv(
        CRC30_DAYS_R,
        CRC30_EVENT_R
    ) ~
        BYPASS_HIGH_CONFIDENCE_R,
    data = ttntd_dt,
    conf.type = "log-log"
)

p5b <- km_plot(
    tt_sf,
    "B  Conservative CRC TTNTD30",
    "Days from first anti-EGFR exposure",
    x_limit = 1200
)

tt_event <- as.data.table(
    b110$event_summary
)

tt_event[
    ,
    EVENT_RATE := EVENTS / N
]

tt_event[
    ,
    FILTER_TYPE_R := fifelse(
        grepl(
            "^Conservative",
            ENDPOINT
        ),
        "Conservative CRC-compatible",
        "Broad cancer-directed"
    )
]

p5c <- ggplot(
    tt_event,
    aes(
        x = ENDPOINT,
        y = EVENT_RATE,
        fill = FILTER_TYPE_R
    )
) +
    geom_col(
        width = 0.6
    ) +
    scale_fill_manual(
        values = c(
            "Conservative CRC-compatible" = COL_TEAL,
            "Broad cancer-directed" = COL_ORANGE
        )
    ) +
    geom_text(
        aes(
            label = paste0(
                EVENTS,
                "/",
                N
            )
        ),
        vjust = -0.35,
        size = 3
    ) +
    scale_y_continuous(
        labels = percent_format(
            accuracy = 1
        ),
        limits = c(
            0,
            min(
                1,
                max(
                    tt_event$EVENT_RATE
                ) +
                    0.12
            )
        )
    ) +
    labs(
        title = "C  Event completeness across TTNTD definitions",
        x = NULL,
        y = "Event proportion",
        fill = NULL
    ) +
    theme_pub(8.5) +
    theme(
        axis.text.x = element_text(
            angle = 25,
            hjust = 1
        )
    )

tt_summary <- ttntd_dt[
    ,
    .(
        N = .N,
        EVENTS =
            sum(
                CRC30_EVENT_R,
                na.rm = TRUE
            ),
        MEDIAN_DAYS =
            median(
                CRC30_DAYS_R,
                na.rm = TRUE
            )
    ),
    by = BYPASS_HIGH_CONFIDENCE_R
][
    ,
    GROUP := fifelse(
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        "High-confidence bypass+",
        "High-confidence bypass−"
    )
][
    ,
    median_label := sprintf("%.1f", MEDIAN_DAYS)
]

p5d <- ggplot(
    tt_summary,
    aes(
        x = GROUP,
        y = MEDIAN_DAYS,
        fill = GROUP
    )
) +
    geom_col(
        width = 0.55
    ) +
    scale_fill_manual(
        values = c(
            "High-confidence bypass−" = COL_NEG,
            "High-confidence bypass+" = COL_POS
        )
    ) +
    geom_text(
        aes(label = median_label),
        vjust = -0.3,
        size = 3.2
    ) +
    labs(
        title = "D  Raw median conservative TTNTD30",
        x = NULL,
        y = "Median days"
    ) +
    theme_pub(9) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(
            angle = 15,
            hjust = 1
        )
    )

fig5 <- (p5a | p5b) /
    (p5c | p5d)

save_pub_figure(
    fig5,
    "Figure5_TTNTD_validation",
    width = 11.2,
    height = 8.7
)


# ============================================================
# 8. TABLE 1: baseline characteristics
# ============================================================

dt[
    ,
    BYPASS_GROUP_R := fifelse(
        BYPASS_HIGH_CONFIDENCE_R == 1L,
        "High-confidence bypass+",
        "High-confidence bypass−"
    )
]

group_levels <- c(
    "High-confidence bypass−",
    "High-confidence bypass+"
)

n_by_group <- dt[
    ,
    .N,
    by = BYPASS_GROUP_R
]

get_group_n <- function(g) {
    n_by_group[
        BYPASS_GROUP_R == g,
        N
    ]
}

table1_rows <- list()

add_cat_row <- function(
    variable,
    level,
    condition
) {
    counts <- vapply(
        group_levels,
        function(g) {
            d <- dt[
                BYPASS_GROUP_R == g
            ]
            sum(
                condition(d),
                na.rm = TRUE
            )
        },
        integer(1)
    )

    dens <- vapply(
        group_levels,
        get_group_n,
        integer(1)
    )

    data.table(
        Characteristic = variable,
        Level = level,
        `Bypass−` =
            n_pct_text(
                counts[1],
                dens[1]
            ),
        `Bypass+` =
            n_pct_text(
                counts[2],
                dens[2]
            )
    )
}

add_cont_row <- function(
    variable,
    extractor,
    digits = 1
) {
    vals <- vapply(
        group_levels,
        function(g) {
            d <- dt[
                BYPASS_GROUP_R == g
            ]
            median_iqr_text(
                extractor(d),
                digits = digits
            )
        },
        character(1)
    )

    data.table(
        Characteristic = variable,
        Level = "Median [IQR]",
        `Bypass−` = vals[1],
        `Bypass+` = vals[2]
    )
}

table1_rows[[length(table1_rows) + 1L]] <-
    add_cat_row(
        "Sex",
        "Male",
        function(d) d$GENDER == "Male"
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cat_row(
        "Primary subsite",
        "Rectum/rectosigmoid",
        function(d) d$SUBSITE_FROZEN_R == "rectum/rectosigmoid"
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cat_row(
        "Anti-EGFR agent",
        "Panitumumab",
        function(d) d$ANTI_EGFR_AGENT_R == "PANITUMUMAB"
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cat_row(
        "Prior cytotoxic exposure",
        "Both oxaliplatin + irinotecan",
        function(d) d$PRIOR_HEAVY_R == 1L
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cat_row(
        "Sample type",
        "Metastasis",
        function(d) grepl(
            "metast",
            safe_chr(
                d$SAMPLE_TYPE
            ),
            ignore.case = TRUE
        )
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cat_row(
        "Gene panel",
        "IMPACT341/410",
        function(d) d$PANEL_OLD_R == 1L
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cont_row(
        "Stage IV diagnosis to anti-EGFR, days",
        function(d) d$STAGE4_TO_ANTI_EGFR_DAYS_R,
        digits = 1
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cont_row(
        "Specimen acquisition to anti-EGFR, days",
        function(d) d$SPECIMEN_TO_ANTI_EGFR_DAYS_R,
        digits = 1
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cont_row(
        "Sequencing to anti-EGFR, days",
        function(d) d$SEQUENCING_TO_ANTI_EGFR_DAYS_R,
        digits = 1
    )

table1_rows[[length(table1_rows) + 1L]] <-
    add_cont_row(
        "Recent metastatic-site count",
        function(d) d$MET_SITE_COUNT_90D_R,
        digits = 1
    )

table1 <- rbindlist(
    table1_rows,
    fill = TRUE
)

setnames(
    table1,
    c(
        "Bypass−",
        "Bypass+"
    ),
    c(
        paste0(
            "Bypass− (n=",
            get_group_n(
                "High-confidence bypass−"
            ),
            ")"
        ),
        paste0(
            "Bypass+ (n=",
            get_group_n(
                "High-confidence bypass+"
            ),
            ")"
        )
    )
)

write_word_table(
    table1,
    "Table 1. Baseline characteristics of the locked anti-EGFR cohort",
    "Table1_baseline_characteristics.docx",
    footnotes = c(
        "Values are n (%) unless otherwise specified.",
        "No significance testing is shown in Table 1; the purpose is descriptive baseline characterization.",
        "Age at anti-EGFR index is unavailable in MSK-CHORD and is therefore not reported."
    )
)


# ============================================================
# 9. TABLE 2: primary + secondary results
# ============================================================

table2_main <- as.data.table(
    b107$primary_results
)[
    ,
    .(
        Endpoint = OUTCOME,
        Analysis = ANALYSIS,
        N,
        Events = EVENTS,
        `HR (95% CI)` = paste0(
            b111_fmt_num(
                HR,
                2
            ),
            " (",
            b111_fmt_num(
                LCL95,
                2
            ),
            "–",
            b111_fmt_num(
                UCL95,
                2
            ),
            ")"
        ),
        `P value` =
            vapply(
                P,
                fmt_p,
                character(1)
            )
    )
]

tt_primary <- as.data.table(
    b110$model_results
)[
    ENDPOINT %chin%
        c(
            "Conservative CRC TTNTD30",
            "Conservative CRC TTNTD60"
        )
][
    ,
    .(
        Endpoint = ENDPOINT,
        Analysis = MODEL,
        N,
        Events = EVENTS,
        `HR (95% CI)` = paste0(
            b111_fmt_num(
                HR,
                2
            ),
            " (",
            b111_fmt_num(
                LCL95,
                2
            ),
            "–",
            b111_fmt_num(
                UCL95,
                2
            ),
            ")"
        ),
        `P value` =
            vapply(
                P,
                fmt_p,
                character(1)
            )
    )
]

table2 <- rbindlist(
    list(
        table2_main,
        tt_primary
    ),
    fill = TRUE
)

write_word_table(
    table2,
    "Table 2. Primary, secondary, and orthogonal treatment-based outcome analyses",
    "Table2_primary_secondary_results.docx",
    footnotes = c(
        "Hazard ratios >1 indicate worse outcomes for patients with high-confidence EGFR-bypass alterations.",
        "rwPFS is the primary endpoint; OS is secondary.",
        "TTNTD analyses are independent treatment-based validation endpoints and are not external validation."
    )
)


# ============================================================
# 10. TABLE S1: bypass definitions / components
# ============================================================

s1 <- data.table(
    `Exposure tier` = c(
        rep(
            "PRESSING-like strict core",
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

    `Role in manuscript` = c(
        rep(
            "Included in strict and primary high-confidence definitions",
            7
        ),
        rep(
            "Included in primary high-confidence definition",
            3
        ),
        "Broad sensitivity only; not primary",
        "Not classified as bypass solely by BRAF; exclusion sensitivity"
    )
)

write_word_table(
    s1,
    "Supplementary Table S1. Frozen EGFR-bypass component definitions",
    "TableS1_bypass_component_definitions.docx",
    footnotes = c(
        "The primary high-confidence exposure definition was frozen before the first formal outcome analysis.",
        "Component-level effects are not interpreted as gene-specific efficacy estimates because individual counts are sparse."
    )
)


# ============================================================
# 11. TABLE S2: robustness
# ============================================================

s2_primary <- as.data.table(
    b108$primary_result
)[
    ,
    .(
        Analysis = "Primary adjusted rwPFS",
        N,
        Events = EVENTS,
        `Bypass+` = BYPASS_POS,
        `HR (95% CI)` = paste0(
            b111_fmt_num(
                HR,
                2
            ),
            " (",
            b111_fmt_num(
                LCL95,
                2
            ),
            "–",
            b111_fmt_num(
                UCL95,
                2
            ),
            ")"
        ),
        `P value` =
            vapply(
                P,
                fmt_p,
                character(1)
            )
    )
]

s2_loo <- as.data.table(
    b108$leave_one_component_out
)[
    ,
    .(
        Analysis = paste0(
            "Leave out ",
            gsub(
                "_R$",
                "",
                COMPONENT
            )
        ),
        N,
        Events = EVENTS,
        `Bypass+` = BYPASS_POS,
        `HR (95% CI)` = paste0(
            b111_fmt_num(
                HR,
                2
            ),
            " (",
            b111_fmt_num(
                LCL95,
                2
            ),
            "–",
            b111_fmt_num(
                UCL95,
                2
            ),
            ")"
        ),
        `P value` =
            vapply(
                P,
                fmt_p,
                character(1)
            )
    )
]

s2_landmark <- as.data.table(
    b108$landmark
)[
    ,
    .(
        Analysis = ANALYSIS,
        N,
        Events = EVENTS,
        `Bypass+` = BYPASS_POS,
        `HR (95% CI)` = paste0(
            b111_fmt_num(
                HR,
                2
            ),
            " (",
            b111_fmt_num(
                LCL95,
                2
            ),
            "–",
            b111_fmt_num(
                UCL95,
                2
            ),
            ")"
        ),
        `P value` =
            vapply(
                P,
                fmt_p,
                character(1)
            )
    )
]

s2 <- rbindlist(
    list(
        s2_primary,
        s2_loo,
        s2_landmark
    ),
    fill = TRUE
)

write_word_table(
    s2,
    "Supplementary Table S2. Post-unblinding robustness analyses for rwPFS",
    "TableS2_robustness_results.docx",
    footnotes = c(
        "These analyses were explicitly post-unblinding robustness checks and did not redefine the frozen primary exposure."
    )
)


# ============================================================
# 12. TABLE S3: refined TTNTD
# ============================================================

s3 <- as.data.table(
    b110$model_results
)[
    ,
    .(
        Endpoint = ENDPOINT,
        Model = MODEL,
        N,
        Events = EVENTS,
        `HR (95% CI)` = paste0(
            b111_fmt_num(
                HR,
                2
            ),
            " (",
            b111_fmt_num(
                LCL95,
                2
            ),
            "–",
            b111_fmt_num(
                UCL95,
                2
            ),
            ")"
        ),
        `P value` =
            vapply(
                P,
                fmt_p,
                character(1)
            )
    )
]

write_word_table(
    s3,
    "Supplementary Table S3. Refined treatment-based TTNTD validation",
    "TableS3_TTNTD_validation.docx",
    footnotes = c(
        "Conservative endpoints count only new CRC-compatible systemic/targeted treatment or death.",
        "PH-addressed models stratify by primary subsite and prior heavy cytotoxic exposure."
    )
)


# ============================================================
# 13. Audit
# ============================================================

audit_file <- file.path(
    audit_dir,
    "B1_11_manuscript_outputs_audit.txt"
)

generated <- c(
    file.path(
        fig_dir,
        paste0(
            "Figure",
            1:5,
            c(
                "_study_design_and_cohort_flow",
                "_high_confidence_bypass_landscape",
                "_primary_rwPFS_OS",
                "_robustness_validation",
                "_TTNTD_validation"
            ),
            ".png"
        )
    ),
    file.path(
        fig_dir,
        paste0(
            "Figure",
            1:5,
            c(
                "_study_design_and_cohort_flow",
                "_high_confidence_bypass_landscape",
                "_primary_rwPFS_OS",
                "_robustness_validation",
                "_TTNTD_validation"
            ),
            ".tiff"
        )
    ),
    file.path(
        table_dir,
        c(
            "Table1_baseline_characteristics.docx",
            "Table2_primary_secondary_results.docx",
            "TableS1_bypass_component_definitions.docx",
            "TableS2_robustness_results.docx",
            "TableS3_TTNTD_validation.docx"
        )
    )
)

missing_outputs <- generated[
    !file.exists(
        generated
    )
]

lines <- c(
    "B1-11 MANUSCRIPT OUTPUTS AUDIT",
    "==============================",
    "",
    paste0(
        "Locked cohort: ",
        nrow(dt)
    ),
    paste0(
        "High-confidence bypass+: ",
        dt[
            BYPASS_HIGH_CONFIDENCE_R == 1L,
            .N
        ]
    ),
    paste0(
        "High-confidence bypass-: ",
        dt[
            BYPASS_HIGH_CONFIDENCE_R == 0L,
            .N
        ]
    ),
    "",
    "Generated files:",
    generated,
    "",
    paste0(
        "Missing outputs: ",
        length(
            missing_outputs
        )
    ),
    if (
        length(
            missing_outputs
        ) > 0L
    ) {
        missing_outputs
    } else {
        "All expected manuscript outputs exist."
    },
    "",
    "No new primary analysis was introduced in B1-11."
)

writeLines(
    lines,
    audit_file
)

cat("\n============================================================\n")
cat("B1-11 COMPLETE\n")
cat("============================================================\n")
cat(
    "Expected manuscript files: ",
    length(generated),
    "\n",
    sep = ""
)
cat(
    "Missing outputs: ",
    length(missing_outputs),
    "\n",
    sep = ""
)
cat(
    "\nAudit:\n",
    audit_file,
    "\n"
)
