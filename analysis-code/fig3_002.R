# Alternative ggplot layouts using the exact ml_theme() from results.qmd.
# Run from analysis-code: Rscript fig3_002.R
# Requires fig3_001.R and the report caches. Produces both A and B in PNG/PDF,
# each with country flags and without flags. Flag-free versions retain country
# points and their faint intervals. Panel B uses condition facets for clarity.
# Optional: FIG3_OUTPUT_DIR sets the destination; FIG3_FLAG_DIR reuses flag assets.
local({
  previous <- options(fig3.data_only = TRUE)
  on.exit(options(previous), add = TRUE)
  source("fig3_001.R", local = TRUE)
  if (!use_flag_images)
    message("Flag images unavailable: the flagged variant will use country dots.")
  for (flags in c(TRUE, FALSE)) {
    suffix <- if (flags) "flags" else "no_flags"
    figures <- make_figures(flags = flags, report_style = TRUE, facets = TRUE)
    save_figures(figures, paste0("fig3_002_", suffix), report_style = TRUE)
  }
})
