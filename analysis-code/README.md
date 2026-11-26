# analysis-code

Analyses, results, and figures for the main study.

## Data

- `all_withNorway_noTest_countryFill_raw.csv`: cleaned full dataset, copied here by `../data-preparation/prepare_dataset.R`.
- `codebook_full_data.xlsx`: codebook for the full dataset.
- `GCI.csv`: country-level index data used in the analyses.
- `TM_WORLD_BORDERS_SIMPL-0.3.*`: world borders shapefile used for maps.

## Analysis files

- `analysis_code.Rmd` / `.html`: planned analysis code.
- `results.qmd` (rendered to `.docx`, `.html`, `.pdf`): results, robustness checks, and exploratory analyses.
- `measurement_invariance_visualizemi.Rmd` / `.html`: measurement invariance analyses and visualizations.
- `scrap_descriptives.R`: scratch descriptive tables.

## Figures

- `fig3_001.R`, `fig3_002.R`: alternative drafts of the main-text Figure 3. They read the cached models and data created by `results.qmd`, so render that first.
- `fig3_002_*.png` / `.pdf`: the output of those scripts (flags / no flags, options A and B).
