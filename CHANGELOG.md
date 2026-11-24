# Changelog

## v1.2.0 – Major Revision 1 (2026-11-24)

Added `code/22_revision_line_denovo_confounding.R`, containing post hoc sensitivity analyses requested during peer review at BMC Gastroenterology:

- Prior treatment exposure assessed with a cytotoxic-backbone proxy and a timeline-derived drug-class proxy (not chart-verified lines of therapy)
- Registry-based stage IV classification at initial diagnosis, with radiology corroboration within 120 days
- Additional ECOG-adjusted complete-case sensitivity analyses (90- and 180-day windows, with and without ECOG in the same patients)
- Proportional-hazards term diagnostics and panel-stratified Cox refits
- Extended-propensity overlap weighting with a 2,000-resample full estimation-pipeline bootstrap
- Treatment-by-bypass cell sample sizes, effective sample sizes and within-stratum covariate balance
- E-values for the primary associations
- Supplementary Tables S11–S14 and Supplementary Figure S6

The original locked cohort, exposure definition, endpoints and primary analysis are unchanged; no earlier script was modified.

## v1.1.0 (2026-10-23)

Added `code/21_literature_consistency_analysis.R` (literature-based external consistency analysis). No existing script or result changed.

## v1.0.0 (2026-10-01)

Initial release.
