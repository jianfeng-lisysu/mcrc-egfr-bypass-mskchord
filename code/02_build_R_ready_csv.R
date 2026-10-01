# ============================================================
# A7 - 02_build_R_ready_csv.R
# Corrected for the ACTUAL MSK-CHORD 2024 file structure
#
# R-ready CSV output:
#   A7_MSK_CHORD_CRC_clean_R.csv
#
# All outputs are written directly to 01_raw_data.
# Raw TXT/TAR files are never modified.
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "/path/to/A7"  # set to your local project directory
raw_dir <- file.path(project_root, "01_raw_data")
audit_dir <- file.path(project_root, "06_logs_and_audit")

if (!dir.exists(raw_dir)) stop("Missing raw directory: ", raw_dir)
if (!dir.exists(audit_dir)) {
    dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}
library(data.table)

cat("\n============================================================\n")
cat("A7 BUILD R-READY MSK-CHORD CSV\n")
cat("Corrected MSK-CHORD parser\n")
cat("============================================================\n\n")


# ============================================================
# Helpers
# ============================================================

min_or_na <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    min(x)
}

max_or_na <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    max(x)
}

row_min_na <- function(a, b) {
    out <- pmin(a, b, na.rm = TRUE)
    out[is.infinite(out)] <- NA_real_
    out
}

safe_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
}

collapse_unique <- function(x, sep = "|") {
    x <- safe_chr(x)
    x <- x[nzchar(x)]
    x <- unique(x)
    if (length(x) == 0L) return(NA_character_)
    paste(x, collapse = sep)
}


# ============================================================
# PART B. Locate and temporarily extract MSK-CHORD
# ============================================================

expected_tar <- file.path(
    raw_dir,
    "msk_chord_2024.tar"
)

if (file.exists(expected_tar)) {
    msk_tar <- expected_tar
} else {
    tar_candidates <- list.files(
        raw_dir,
        pattern = "\\.(tar|tar\\.gz|tgz)$",
        full.names = TRUE,
        ignore.case = TRUE
    )

    if (length(tar_candidates) == 0L) {
        stop("No MSK-CHORD TAR archive found in: ", raw_dir)
    }

    msk_tar <- tar_candidates[
        which.max(file.info(tar_candidates)$size)
    ]
}

cat("MSK archive:\n", msk_tar, "\n\n")

tmp_dir <- tempfile("A7_MSK_CHORD_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

untar(
    msk_tar,
    exdir = tmp_dir
)

all_msk_files <- list.files(
    tmp_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
)

find_msk <- function(filename) {
    hit <- all_msk_files[
        basename(all_msk_files) == filename
    ]

    if (length(hit) == 0L) {
        stop("Missing MSK file: ", filename)
    }

    hit[1]
}


# ============================================================
# PART C. Read ACTUAL MSK-CHORD files
#
# data_clinical_patient.txt and data_clinical_sample.txt
# contain four cBioPortal metadata rows before the real header.
# ============================================================

cat("Reading MSK clinical files...\n")

patient <- fread(
    find_msk("data_clinical_patient.txt"),
    sep = "\t",
    skip = 4,
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

sample <- fread(
    find_msk("data_clinical_sample.txt"),
    sep = "\t",
    skip = 4,
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

required_patient <- c(
    "PATIENT_ID",
    "STAGE_HIGHEST_RECORDED",
    "OS_MONTHS",
    "OS_STATUS"
)

required_sample <- c(
    "SAMPLE_ID",
    "PATIENT_ID",
    "CANCER_TYPE",
    "PRIMARY_SITE",
    "CANCER_TYPE_DETAILED",
    "GENE_PANEL",
    "MSI_TYPE"
)

if (!all(required_patient %in% names(patient))) {
    stop(
        "Patient clinical header still not parsed correctly. Missing: ",
        paste(
            setdiff(required_patient, names(patient)),
            collapse = ", "
        )
    )
}

if (!all(required_sample %in% names(sample))) {
    stop(
        "Sample clinical header still not parsed correctly. Missing: ",
        paste(
            setdiff(required_sample, names(sample)),
            collapse = ", "
        )
    )
}

sample_crc <- sample[
    CANCER_TYPE == "Colorectal Cancer"
]

crc_ids <- unique(sample_crc$PATIENT_ID)
crc_samples <- unique(sample_crc$SAMPLE_ID)

patient_crc <- patient[
    PATIENT_ID %chin% crc_ids
]

cat(
    "CRC:",
    nrow(sample_crc),
    "sample rows;",
    uniqueN(sample_crc$PATIENT_ID),
    "patients\n"
)

if (nrow(sample_crc) != 5543L ||
    uniqueN(sample_crc$PATIENT_ID) != 5543L) {
    stop(
        "CRC anchor failed. Expected 5,543 patients / sample rows; got ",
        nrow(sample_crc),
        " rows and ",
        uniqueN(sample_crc$PATIENT_ID),
        " patients."
    )
}

if (uniqueN(sample_crc$SAMPLE_ID) != 5543L) {
    stop(
        "Expected exactly one public genomic sample per CRC patient."
    )
}

# One patient-level row, retaining original clinical/sample variables.
msk <- merge(
    sample_crc,
    patient_crc,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_SAMPLE", "_PATIENT")
)

msk[, OS_EVENT_R := fifelse(
    OS_STATUS == "1:DECEASED",
    1L,
    fifelse(OS_STATUS == "0:LIVING", 0L, NA_integer_)
)]

msk[, OS_DAYS_APPROX_R := OS_MONTHS * 30.4375]

msk[, MUCINOUS_R := as.integer(
    grepl(
        "mucinous",
        CANCER_TYPE_DETAILED,
        ignore.case = TRUE
    )
)]


# ============================================================
# PART D. Diagnosis and sidedness
# ============================================================

diagnosis <- fread(
    find_msk("data_timeline_diagnosis.txt"),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

diagnosis <- diagnosis[
    PATIENT_ID %chin% crc_ids
]

left_dx_codes <- c(
    "C185", "C186", "C187", "C199", "C209"
)

right_dx_codes <- c(
    "C180", "C181", "C182", "C183", "C184"
)

dx_side_one <- function(x) {

    txt <- toupper(
        paste(
            x[!is.na(x)],
            collapse = " "
        )
    )

    if (!nzchar(txt)) return("unknown")

    matches <- gregexpr(
        "C18[0-9]|C19[0-9]|C20[0-9]",
        txt,
        perl = TRUE
    )

    codes <- regmatches(
        txt,
        matches
    )[[1]]

    if (length(codes) == 0L) return("unknown")

    has_left <- any(codes %in% left_dx_codes)
    has_right <- any(codes %in% right_dx_codes)

    if (has_left && !has_right) return("left")
    if (has_right && !has_left) return("right")
    if (has_left && has_right) return("conflict")

    "unknown"
}

diagnosis_sum <- diagnosis[
    order(START_DATE),
    .(
        DX_SIDE_R = dx_side_one(DX_DESCRIPTION),

        DIAGNOSIS_FIRST_DAY_R = min_or_na(START_DATE),
        DIAGNOSIS_LAST_DAY_R = max_or_na(START_DATE),

        DIAGNOSIS_SEQUENCE_R = paste(
            paste(
                START_DATE,
                SUBTYPE,
                SOURCE,
                DX_DESCRIPTION,
                STAGE_CDM_DERIVED,
                SUMMARY,
                sep = "~"
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]

msk <- merge(
    msk,
    diagnosis_sum,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

left_structured <- c(
    "rectum",
    "sigmoid colon",
    "rectosigmoid colon",
    "descending colon",
    "left colon"
)

right_structured <- c(
    "cecum",
    "ascending colon",
    "transverse colon",
    "hepatic flexure"
)

primary_site_lower <- tolower(
    trimws(
        safe_chr(msk$PRIMARY_SITE)
    )
)

msk[, STRUCTURED_SIDE_R := fcase(
    primary_site_lower %chin% left_structured,
    "left",

    primary_site_lower %chin% right_structured,
    "right",

    default = "unknown"
)]

msk[is.na(DX_SIDE_R), DX_SIDE_R := "unknown"]

# Frozen audit rule:
# structured PRIMARY_SITE has priority when definitive;
# diagnosis ICD-O code is salvage only when structured site
# is uninformative.
msk[, PRIMARY_SIDE_R := fifelse(
    STRUCTURED_SIDE_R %chin% c("left", "right"),
    STRUCTURED_SIDE_R,
    fifelse(
        DX_SIDE_R %chin% c("left", "right"),
        DX_SIDE_R,
        "unknown"
    )
)]


# ============================================================
# PART E. Sequencing and specimen acquisition timing
# ============================================================

specimen <- fread(
    find_msk("data_timeline_specimen.txt"),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

specimen <- specimen[
    PATIENT_ID %chin% crc_ids
]

specimen_sum <- specimen[
    order(START_DATE),
    .(
        SEQUENCING_FIRST_DAY_R = min_or_na(START_DATE),
        SEQUENCING_LAST_DAY_R = max_or_na(START_DATE),
        SEQUENCING_SEQUENCE_R = paste(
            paste(
                START_DATE,
                SAMPLE_ID,
                sep = "~"
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]

specimen_surgery <- fread(
    find_msk("data_timeline_specimen_surgery.txt"),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

specimen_surgery <- specimen_surgery[
    PATIENT_ID %chin% crc_ids
]

acq_sum <- specimen_surgery[
    order(START_DATE),
    .(
        SPECIMEN_ACQUISITION_FIRST_DAY_R =
            min_or_na(START_DATE),

        SPECIMEN_ACQUISITION_LAST_DAY_R =
            max_or_na(START_DATE),

        SPECIMEN_ACQUISITION_SEQUENCE_R = paste(
            paste(
                START_DATE,
                SAMPLE_ID,
                SEQ_DATE,
                sep = "~"
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]

msk <- merge(
    msk,
    specimen_sum,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

msk <- merge(
    msk,
    acq_sum,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)


# ============================================================
# PART F. Treatment timeline
# ============================================================

treatment <- fread(
    find_msk("data_timeline_treatment.txt"),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

treatment <- treatment[
    PATIENT_ID %chin% crc_ids
]

treatment[, AGENT_UPPER_R := toupper(
    safe_chr(AGENT)
)]

treatment_sum <- treatment[
    order(START_DATE, STOP_DATE),
    .(
        TREATMENT_EVENT_N_R = .N,

        TREATMENT_FIRST_DAY_R =
            min_or_na(START_DATE),

        TREATMENT_LAST_START_DAY_R =
            max_or_na(START_DATE),

        TREATMENT_SEQUENCE_R = paste(
            paste(
                START_DATE,
                STOP_DATE,
                SUBTYPE,
                AGENT,
                RX_INVESTIGATIVE,
                sep = "~"
            ),
            collapse = "|"
        ),

        CETUXIMAB_FIRST_DAY_R =
            min_or_na(
                START_DATE[
                    AGENT_UPPER_R == "CETUXIMAB"
                ]
            ),

        PANITUMUMAB_FIRST_DAY_R =
            min_or_na(
                START_DATE[
                    AGENT_UPPER_R == "PANITUMUMAB"
                ]
            ),

        BEVACIZUMAB_FIRST_DAY_R =
            min_or_na(
                START_DATE[
                    AGENT_UPPER_R == "BEVACIZUMAB"
                ]
            ),

        OXALIPLATIN_FIRST_DAY_R =
            min_or_na(
                START_DATE[
                    AGENT_UPPER_R == "OXALIPLATIN"
                ]
            ),

        IRINOTECAN_FIRST_DAY_R =
            min_or_na(
                START_DATE[
                    AGENT_UPPER_R == "IRINOTECAN"
                ]
            ),

        FLUOROURACIL_FIRST_DAY_R =
            min_or_na(
                START_DATE[
                    AGENT_UPPER_R == "FLUOROURACIL"
                ]
            ),

        CAPECITABINE_FIRST_DAY_R =
            min_or_na(
                START_DATE[
                    AGENT_UPPER_R == "CAPECITABINE"
                ]
            )
    ),
    by = PATIENT_ID
]

msk <- merge(
    msk,
    treatment_sum,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

msk[, ANTI_EGFR_FIRST_DAY_R := row_min_na(
    CETUXIMAB_FIRST_DAY_R,
    PANITUMUMAB_FIRST_DAY_R
)]

msk[, FIRST_BIOLOGIC_DAY_R := row_min_na(
    ANTI_EGFR_FIRST_DAY_R,
    BEVACIZUMAB_FIRST_DAY_R
)]

msk[, FIRST_BIOLOGIC_R := fifelse(
    is.na(FIRST_BIOLOGIC_DAY_R),
    NA_character_,
    fifelse(
        !is.na(ANTI_EGFR_FIRST_DAY_R) &
            ANTI_EGFR_FIRST_DAY_R == FIRST_BIOLOGIC_DAY_R,
        "anti-EGFR",
        "bevacizumab"
    )
)]


# ============================================================
# PART G. Progression
# ============================================================

progression <- fread(
    find_msk("data_timeline_progression.txt"),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

progression <- progression[
    PATIENT_ID %chin% crc_ids
]

progression_sum <- progression[
    order(START_DATE),
    .(
        PROGRESSION_EVENT_N_R = .N,
        PROGRESSION_FIRST_DAY_R =
            min_or_na(START_DATE),
        PROGRESSION_LAST_DAY_R =
            max_or_na(START_DATE),

        PROGRESSION_SEQUENCE_R = paste(
            paste(
                START_DATE,
                PROGRESSION,
                NLP_PROGRESSION_PROBABILITY,
                PROCEDURE_TYPE,
                sep = "~"
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]

msk <- merge(
    msk,
    progression_sum,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)


# ============================================================
# PART H. Surgery and tumor-site timelines
# ============================================================

surgery <- fread(
    find_msk("data_timeline_surgery.txt"),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

surgery <- surgery[
    PATIENT_ID %chin% crc_ids
]

surgery_sum <- surgery[
    order(START_DATE),
    .(
        SURGERY_EVENT_N_R = .N,
        SURGERY_FIRST_DAY_R =
            min_or_na(START_DATE),
        SURGERY_LAST_DAY_R =
            max_or_na(START_DATE),

        SURGERY_SEQUENCE_R = paste(
            paste(
                START_DATE,
                SUBTYPE,
                sep = "~"
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]

tumor_sites <- fread(
    find_msk("data_timeline_tumor_sites.txt"),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

tumor_sites <- tumor_sites[
    PATIENT_ID %chin% crc_ids
]

tumor_site_sum <- tumor_sites[
    order(START_DATE),
    .(
        TUMOR_SITE_EVENT_N_R = .N,
        TUMOR_SITE_FIRST_DAY_R =
            min_or_na(START_DATE),
        TUMOR_SITE_LAST_DAY_R =
            max_or_na(START_DATE),

        TUMOR_SITE_SEQUENCE_R = paste(
            paste(
                START_DATE,
                SOURCE,
                SOURCE_SPECIFIC,
                TUMOR_SITE,
                CHEST,
                ABDOMEN,
                PELVIS,
                HEAD,
                OTHER,
                sep = "~"
            ),
            collapse = "|"
        )
    ),
    by = PATIENT_ID
]

msk <- merge(
    msk,
    surgery_sum,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

msk <- merge(
    msk,
    tumor_site_sum,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)


# ============================================================
# PART I. Mutation data
# ============================================================

mutation_file <- find_msk(
    "data_mutations.txt"
)

cat("Reading key mutation columns...\n")

mut <- fread(
    mutation_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    select = c(
        "Hugo_Symbol",
        "Tumor_Sample_Barcode",
        "HGVSp",
        "HGVSp_Short",
        "Variant_Classification"
    ),
    showProgress = TRUE
)

mut <- mut[
    Tumor_Sample_Barcode %chin% crc_samples
]

sample_map <- sample_crc[
    ,
    .(
        SAMPLE_ID,
        PATIENT_ID
    )
]

mut[
    ,
    PATIENT_ID := sample_map$PATIENT_ID[
        match(
            Tumor_Sample_Barcode,
            sample_map$SAMPLE_ID
        )
    ]
]

key_genes <- c(
    "KRAS",
    "NRAS",
    "BRAF",
    "MAP2K1",
    "PIK3CA",
    "GNAS",
    "NF1"
)

mut_key <- mut[
    Hugo_Symbol %chin% key_genes
]

mut_flag_long <- unique(
    mut_key[
        ,
        .(
            PATIENT_ID,
            GENE = Hugo_Symbol,
            FLAG = 1L
        )
    ]
)

if (nrow(mut_flag_long) > 0L) {

    mut_flags <- dcast(
        mut_flag_long,
        PATIENT_ID ~ GENE,
        value.var = "FLAG",
        fill = 0L
    )

    for (g in key_genes) {

        if (!g %in% names(mut_flags)) {
            mut_flags[, (g) := 0L]
        }

        setnames(
            mut_flags,
            g,
            paste0(g, "_MUTATION_R")
        )
    }

} else {

    mut_flags <- data.table(
        PATIENT_ID = character()
    )

    for (g in key_genes) {
        mut_flags[, (paste0(g, "_MUTATION_R")) := integer()]
    }
}

mut_details <- mut_key[
    ,
    .(
        MUTATION_DETAIL = collapse_unique(
            paste(
                HGVSp_Short,
                Variant_Classification,
                sep = ":"
            ),
            sep = ";"
        )
    ),
    by = .(
        PATIENT_ID,
        Hugo_Symbol
    )
]

if (nrow(mut_details) > 0L) {

    mut_details_wide <- dcast(
        mut_details,
        PATIENT_ID ~ Hugo_Symbol,
        value.var = "MUTATION_DETAIL"
    )

    detail_genes <- setdiff(
        names(mut_details_wide),
        "PATIENT_ID"
    )

    setnames(
        mut_details_wide,
        detail_genes,
        paste0(
            detail_genes,
            "_VARIANTS_R"
        )
    )

} else {

    mut_details_wide <- data.table(
        PATIENT_ID = character()
    )
}

braf_v600e <- unique(
    mut[
        Hugo_Symbol == "BRAF" &
            (
                grepl(
                    "V600E",
                    safe_chr(HGVSp_Short),
                    ignore.case = TRUE
                ) |
                grepl(
                    "V600E",
                    safe_chr(HGVSp),
                    ignore.case = TRUE
                )
            ),
        .(
            PATIENT_ID,
            BRAF_V600E_R = 1L
        )
    ]
)

msk <- merge(
    msk,
    mut_flags,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

msk <- merge(
    msk,
    mut_details_wide,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

msk <- merge(
    msk,
    braf_v600e,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

for (g in key_genes) {
    v <- paste0(g, "_MUTATION_R")

    if (!v %in% names(msk)) {
        msk[, (v) := 0L]
    }

    msk[
        is.na(get(v)),
        (v) := 0L
    ]
}

if (!"BRAF_V600E_R" %in% names(msk)) {
    msk[, BRAF_V600E_R := 0L]
}
msk[
    is.na(BRAF_V600E_R),
    BRAF_V600E_R := 0L
]


# ============================================================
# PART J. CNA data
#
# ACTUAL data_cna.txt is a WIDE matrix:
#   Hugo_Symbol x sample columns
# ============================================================

cna_file <- find_msk(
    "data_cna.txt"
)

cat("Reading ERBB2/KRAS CNA from wide matrix...\n")

cna_header <- names(
    fread(
        cna_file,
        sep = "\t",
        header = TRUE,
        nrows = 0,
        check.names = FALSE,
        showProgress = FALSE
    )
)

crc_cna_samples <- intersect(
    crc_samples,
    cna_header
)

if (length(crc_cna_samples) == 0L) {
    stop("No CRC sample IDs found in CNA matrix header.")
}

cna_small <- fread(
    cna_file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    select = c(
        "Hugo_Symbol",
        crc_cna_samples
    ),
    showProgress = TRUE
)

cna_small <- cna_small[
    Hugo_Symbol %chin% c(
        "ERBB2",
        "KRAS"
    )
]

amp_long_list <- lapply(
    c(
        "ERBB2",
        "KRAS"
    ),
    function(g) {

        row_g <- cna_small[
            Hugo_Symbol == g
        ]

        if (nrow(row_g) == 0L) {

            return(
                data.table(
                    SAMPLE_ID = crc_cna_samples,
                    GENE = g,
                    AMP = 0L
                )
            )
        }

        vals <- suppressWarnings(
            as.numeric(
                unlist(
                    row_g[
                        1,
                        ..crc_cna_samples
                    ],
                    use.names = FALSE
                )
            )
        )

        data.table(
            SAMPLE_ID = crc_cna_samples,
            GENE = g,
            AMP = as.integer(
                vals == 2
            )
        )
    }
)

amp_long <- rbindlist(
    amp_long_list
)

amp_wide <- dcast(
    amp_long,
    SAMPLE_ID ~ GENE,
    value.var = "AMP",
    fill = 0L
)

setnames(
    amp_wide,
    "ERBB2",
    "ERBB2_AMP_R"
)

setnames(
    amp_wide,
    "KRAS",
    "KRAS_AMP_R"
)

amp_wide[
    ,
    PATIENT_ID := sample_map$PATIENT_ID[
        match(
            SAMPLE_ID,
            sample_map$SAMPLE_ID
        )
    ]
]

amp_patient <- amp_wide[
    ,
    .(
        ERBB2_AMP_R = max(
            ERBB2_AMP_R,
            na.rm = TRUE
        ),
        KRAS_AMP_R = max(
            KRAS_AMP_R,
            na.rm = TRUE
        )
    ),
    by = PATIENT_ID
]

msk <- merge(
    msk,
    amp_patient,
    by = "PATIENT_ID",
    all.x = TRUE,
    sort = FALSE
)

msk[
    is.na(ERBB2_AMP_R),
    ERBB2_AMP_R := 0L
]

msk[
    is.na(KRAS_AMP_R),
    KRAS_AMP_R := 0L
]


# ============================================================
# PART K. Panel metadata
# ============================================================

panel_matrix <- fread(
    find_msk("data_gene_panel_matrix.txt"),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    showProgress = FALSE
)

panel_crc <- panel_matrix[
    SAMPLE_ID %chin% crc_samples,
    .(
        SAMPLE_ID,
        MUTATION_PANEL_R = mutations,
        CNA_PANEL_R = cna,
        STRUCTURAL_VARIANT_PANEL_R =
            structural_variants
    )
]

msk <- merge(
    msk,
    panel_crc,
    by = "SAMPLE_ID",
    all.x = TRUE,
    sort = FALSE
)


# ============================================================
# PART L. Original bypass composite + old feasibility anchors
# ============================================================

msk[, BYPASS_ORIGINAL_R := as.integer(
    ERBB2_AMP_R == 1L |
        KRAS_AMP_R == 1L |
        MAP2K1_MUTATION_R == 1L |
        PIK3CA_MUTATION_R == 1L |
        GNAS_MUTATION_R == 1L |
        NF1_MUTATION_R == 1L
)]

msk[, ANY_RAS_MUTATION_R := as.integer(
    KRAS_MUTATION_R == 1L |
        NRAS_MUTATION_R == 1L
)]

msk[, OLD_ELIGIBLE_803_R := as.integer(
    PRIMARY_SIDE_R == "left" &
        STAGE_HIGHEST_RECORDED == "Stage 4" &
        MSI_TYPE == "Stable" &
        ANY_RAS_MUTATION_R == 0L &
        BRAF_V600E_R == 0L
)]

msk[, OLD_BIOLOGIC_458_R := as.integer(
    OLD_ELIGIBLE_803_R == 1L &
        !is.na(FIRST_BIOLOGIC_R)
)]


# ============================================================
# PART M. Final hard audit BEFORE writing MSK CSV
# ============================================================

n_crc <- nrow(msk)

n_eligible <- msk[
    OLD_ELIGIBLE_803_R == 1L,
    .N
]

n_biologic <- msk[
    OLD_BIOLOGIC_458_R == 1L,
    .N
]

n_egfr <- msk[
    OLD_BIOLOGIC_458_R == 1L &
        FIRST_BIOLOGIC_R == "anti-EGFR",
    .N
]

n_bev <- msk[
    OLD_BIOLOGIC_458_R == 1L &
        FIRST_BIOLOGIC_R == "bevacizumab",
    .N
]

n_egfr_bp_pos <- msk[
    OLD_BIOLOGIC_458_R == 1L &
        FIRST_BIOLOGIC_R == "anti-EGFR" &
        BYPASS_ORIGINAL_R == 1L,
    .N
]

n_egfr_bp_neg <- msk[
    OLD_BIOLOGIC_458_R == 1L &
        FIRST_BIOLOGIC_R == "anti-EGFR" &
        BYPASS_ORIGINAL_R == 0L,
    .N
]

n_bev_bp_pos <- msk[
    OLD_BIOLOGIC_458_R == 1L &
        FIRST_BIOLOGIC_R == "bevacizumab" &
        BYPASS_ORIGINAL_R == 1L,
    .N
]

n_bev_bp_neg <- msk[
    OLD_BIOLOGIC_458_R == 1L &
        FIRST_BIOLOGIC_R == "bevacizumab" &
        BYPASS_ORIGINAL_R == 0L,
    .N
]

cat("\n============================================================\n")
cat("MSK HARD AUDIT\n")
cat("============================================================\n")
cat("CRC patients                         :", n_crc, "\n")
cat("Old eligible                         :", n_eligible, "\n")
cat("Old biologic-treated                 :", n_biologic, "\n")
cat("anti-EGFR first / BEV first          :", n_egfr, "/", n_bev, "\n")
cat(
    "anti-EGFR bypass- / bypass+          :",
    n_egfr_bp_neg,
    "/",
    n_egfr_bp_pos,
    "\n"
)
cat(
    "BEV bypass- / bypass+                :",
    n_bev_bp_neg,
    "/",
    n_bev_bp_pos,
    "\n"
)

anchor_ok <- (
    n_crc == 5543L &&
        n_eligible == 803L &&
        n_biologic == 458L &&
        n_egfr == 230L &&
        n_bev == 228L &&
        n_egfr_bp_neg == 189L &&
        n_egfr_bp_pos == 41L &&
        n_bev_bp_neg == 176L &&
        n_bev_bp_pos == 52L
)

if (!anchor_ok) {

    stop(
        paste0(
            "\nHARD AUDIT FAILED.\n",
            "The cleaned MSK table does not reproduce the raw-data anchors.\n",
            "Do NOT use the CSV for analysis."
        )
    )
}

cat("\nPASS: all MSK raw-data anchors reproduced exactly.\n")


# ============================================================
# PART N. Write the ONE unified MSK patient-level CSV
# ============================================================

setorder(
    msk,
    PATIENT_ID
)

msk_out <- file.path(
    raw_dir,
    "A7_MSK_CHORD_CRC_clean_R.csv"
)

fwrite(
    msk,
    msk_out,
    na = ""
)

cat(
    "\nMSK written:",
    nrow(msk),
    "rows x",
    ncol(msk),
    "columns\n"
)


# ============================================================
# PART O. Audit report
# ============================================================

audit_file <- file.path(
    audit_dir,
    "02_build_two_clean_R_csv_audit.txt"
)

audit_lines <- c(
    "A7 R-READY MSK-CHORD CSV BUILD AUDIT",
    "===================================",
    "",
    paste0(
        "MSK output: ",
        msk_out
    ),
    paste0(
        "MSK CRC patients: ",
        n_crc
    ),
    paste0(
        "MSK columns: ",
        ncol(msk)
    ),
    "",
    "Hard anchors:",
    paste0(
        "CRC = ",
        n_crc,
        " [expected 5543]"
    ),
    paste0(
        "Old eligible = ",
        n_eligible,
        " [expected 803]"
    ),
    paste0(
        "Old biologic = ",
        n_biologic,
        " [expected 458]"
    ),
    paste0(
        "anti-EGFR / BEV = ",
        n_egfr,
        " / ",
        n_bev,
        " [expected 230 / 228]"
    ),
    paste0(
        "anti-EGFR bypass- / bypass+ = ",
        n_egfr_bp_neg,
        " / ",
        n_egfr_bp_pos,
        " [expected 189 / 41]"
    ),
    paste0(
        "BEV bypass- / bypass+ = ",
        n_bev_bp_neg,
        " / ",
        n_bev_bp_pos,
        " [expected 176 / 52]"
    ),
    "",
    paste0(
        "Hard audit PASS: ",
        anchor_ok
    ),
    "",
    "R-ready CSV file used:",
    "A7_MSK_CHORD_CRC_clean_R.csv",
    "",
    "No treatment-effect model was run."
)

writeLines(
    audit_lines,
    audit_file
)

unlink(
    tmp_dir,
    recursive = TRUE,
    force = TRUE
)

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
cat("A7_MSK_CHORD_CRC_clean_R.csv\n")
cat("\nAudit:\n", audit_file, "\n")
