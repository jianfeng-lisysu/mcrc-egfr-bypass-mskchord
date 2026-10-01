# Variant-level bypass alterations and outcomes after EGFR-targeted therapy in molecularly selected metastatic colorectal cancer

Analysis code for the manuscript *Variant-Level Bypass Alterations and Outcomes After Epidermal Growth Factor Receptor–Targeted Therapy in Molecularly Selected Metastatic Colorectal Cancer: An Exploratory Real-World Clinicogenomic Study*.

## What this repository contains

R scripts for cohort derivation, the locked genomic exposure definition, the primary and supportive outcome analyses, the exploratory bevacizumab comparator, and all sensitivity and robustness analyses reported in the paper and its supplementary materials.

## What this repository does not contain

No patient-level data of any kind. The analysis uses the publicly released 2024 MSK-CHORD dataset, distributed through cBioPortal under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International licence. That licence does not permit redistribution of derived patient-level files, so no raw data, no intermediate patient-level tables and no patient-level audit files are included here.

To reproduce the analyses, download the source dataset yourself:

- Study: `msk_chord_2024`
- URL: https://www.cbioportal.org/study/summary?id=msk_chord_2024
- Source publication: Jee J, Fong C, Pichotta K, et al. Automated real-world data integration improves cancer outcome prediction. *Nature* 2024;636:728–736. doi:10.1038/s41586-024-08167-5

Then set `project_root` (in `14C_*`, `root`) at the top of each script, currently the placeholder `"/path/to/A7"`, to your own directory and place the downloaded files in `01_raw_data/`.

## Script order

Scripts are numbered in execution order. The numbering also records the sequence in which the study was conducted, which matters for two claims made in the paper.

| Script | Purpose |
|---|---|
| `02_build_R_ready_csv.R` | Build the analysis-ready dataset from the cBioPortal release; define genomic flags, including the original sensitivity composite `BYPASS_ORIGINAL_R` |
| `03_*`, `04_*`, `05_*` | Timeline reconstruction, first-relevant-biologic assignment, feasibility audits |
| `06_B1_06_freeze_bypass_covariates.R` | **Locks the high-confidence EGFR-bypass definition and the covariate set.** Run before any outcome model |
| `07_B1_07_first_unblinded_analysis.R` | First formal outcome analysis: primary rwPFS and secondary OS |
| `08B_*`, `09_*`, `10_*` | Robustness analyses and the supportive treatment-based TTNTD endpoint ("validation" in these file names refers to endpoint auditing, not external validation) |
| `11B_*`, `11C_*`, `12_*` | Manuscript figures and tables |
| `13_*`, `13A_*`, `13C_*` | Post-hoc treatment-context adjustment, proportional-hazards handling, specimen-timing diagnostics |
| `14_*` | Outcome-blind bevacizumab comparator feasibility gate |
| `14B_*`, `14C_*` | Exploratory bevacizumab comparator and full estimation-pipeline bootstrap |
| `16_*`, `16A_*` | Shared specimen-recency stress test; revised Figure 4 |
| `18_*`, `18A_*` | Supplementary Tables S4, S7 and S8 |
| `19_*` | Figure 5 relabelling (no model refitted) |

Two points that the numbering documents:

1. The genomic exposure was constructed and frozen in `06_*` before any survival model was run in `07_*`. This is what the paper means by "locked", as distinct from prospectively preregistered; there is no registered prospective protocol.
2. The bevacizumab comparator passed an outcome-blind feasibility gate in `14_*` before comparator outcomes were analysed.

## Notes on this release

- Only the final version of each script is included; earlier iterations that were superseded during figure and table refinement are not part of this archive.
- An unrelated data-preparation section for a separate project was removed from `02_build_R_ready_csv.R`; no reported analysis depends on it.
- Absolute local paths were replaced by a placeholder. Descriptive labels in some console and audit outputs were reworded for consistency with the manuscript; no model, cohort definition, or numerical output was changed.
- Table titles and notes in the submitted manuscript were subsequently copy-edited, so file names and captions produced by the scripts may differ slightly from the published versions. Numerical content is unchanged.

## Environment

- R version 4.5.2
- `survival` 3.8.6
- `data.table` 1.18.4
- `survRM2` 1.0.4

## Licence

Code: MIT (see `LICENSE`). The source dataset is licensed separately by its distributor; see above.

## Citation

If you use this code, please cite the paper. A DOI for this archive is issued by Zenodo and shown on the repository page.
