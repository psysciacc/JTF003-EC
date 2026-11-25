# pilots

Pilot studies run before the main data collection.

- `pilot1/`: Qualtrics survey (`.qsf`), pilot data (`.csv`), codebook (`codebook_pilot1.xlsx`), and analysis script (`pilot_analyses.R`).
- `pilot2/`: Qualtrics survey (`.qsf`), pilot data (`ECFullPilot2.csv`), codebook (`codebook_pilot2.xlsx`), and a report (`report_pilot_2.qmd` / `.pdf`).

The pilot codebooks reuse the rows of `../analysis-code/codebook_full_data_updated.xlsx` for matching variable names and add pilot-specific columns (observed ranges, unique counts, and checks against the main codebook's options). Unlike the main dataset, the pilot CSVs store numeric codes, not text labels.
