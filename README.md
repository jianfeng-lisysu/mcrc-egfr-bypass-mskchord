# Variant-level EGFR-bypass alterations and outcomes after anti-EGFR therapy in molecularly selected metastatic colorectal cancer

Analysis code for the manuscript *Variant-Level EGFR-Bypass Alterations and Outcomes After Anti-EGFR Therapy in Molecularly Selected Metastatic Colorectal Cancer*.

## What this repository contains

R scripts for cohort derivation, the locked genomic exposure definition, the primary and supportive outcome analyses, the exploratory bevacizumab comparator, and all sensitivity and robustness analyses reported in the paper and its supplementary materials, and (from v1.1.0) the literature-based external consistency analysis (Supplementary Table S10, Supplementary Figure S4), and (from v1.2.0) the reviewer-requested post hoc analyses for the first revision (Supplementary Tables S11–S14, Supplementary Figure S6).

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
| `21_literature_consistency_analysis.R` | Literature-based external consistency analysis (added in v1.1.0, post hoc): published HRs pooled with DerSimonian-Laird random effects; writes CSV results, three-line Word tables and Supplementary Figure S4. Uses only published summary estimates; needs no patient-level data |
| `22_revision_line_denovo_confounding.R` | Reviewer-requested post hoc analyses for the first revision (added in v1.2.0): registry stage at initial diagnosis with radiology corroboration; prior-treatment-exposure proxies (cytotoxic backbone exposure and a timeline-derived drug-class count; not chart-verified lines of therapy); ECOG and treatment-history sensitivity models; proportional-hazards term diagnostics and panel-stratified refits; extended-propensity overlap weights with a 2,000-resample full estimation-pipeline bootstrap, treatment-by-bypass cell diagnostics and within-stratum balance; E-values. Reads the locked outputs of `07_*` and `14B_*` and the timeline files of the source release; writes Supplementary Tables S11–S14 and Supplementary Figure S6. The locked cohort, exposure definition, endpoints and primary analysis are not changed |

Two points that the numbering documents:

1. The genomic exposure was constructed and frozen in `06_*` before any survival model was run in `07_*`. This is what the paper means by "locked", as distinct from prospectively preregistered; there is no registered prospective protocol.
2. The bevacizumab comparator passed an outcome-blind feasibility gate in `14_*` before comparator outcomes were analysed.

## Version history

- v1.2.0 (2026-11-24): added `22_revision_line_denovo_confounding.R` (post hoc analyses requested during peer review); no existing script or result changed.
- v1.1.0 (2026-10-23): added `21_literature_consistency_analysis.R`; no existing script or result changed.
- v1.0.0 (2026-10-01): initial release.

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
- For `21_*` only: `metafor` 5.0-1, `ggplot2` 4.0.3, `patchwork` 1.3.2, `flextable`, `officer`
- For `22_*`: `ggplot2`, `flextable`, `officer`; the bootstrap (2,000 resamples, seed 20261124) takes about a minute

## Licence

Code: MIT (see `LICENSE`). The source dataset is licensed separately by its distributor; see above.

## Citation

If you use this code, please cite the paper. A DOI for this archive is issued by Zenodo and shown on the repository page.
