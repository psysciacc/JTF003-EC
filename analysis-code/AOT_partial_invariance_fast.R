# Targeted partial measurement-invariance analysis for AOT.
# Put this file alongside AOT_partial_invariance.R and the two input RDS files.
# Run in that directory, for example:
#   source("AOT_partial_invariance_fast.R")
#   result <- run_aot_partial_fast()                 # exploratory: top 2 scalar residual items
#   result <- run_aot_partial_fast(intercept_items = c("AOT_3", "AOT_9"))
#   result <- run_aot_partial_fast(loading_items = "AOT_4",
#                                 intercept_items = "AOT_9")
# Original exhaustive search is NOT called.

source("AOT_partial_invariance.R")

# Work out only the models needed, separately releasing metric vs scalar constraints.
aot_fast_specs <- function(loading_items, intercept_items) {
  wording <- paste0("NEG=~", aot_negative)
  # In R, paste0(prefix, character(0)) can return a spurious string.
  loadings <- if (length(loading_items)) paste0("AOT=~", loading_items) else character()
  intercepts <- if (length(intercept_items)) paste0(intercept_items, "~1") else character()
  specs <- list(
    list(name = "Configural", equal = character(), partial = character()),
    list(name = "Metric", equal = "loadings", partial = wording),
    list(name = "Scalar", equal = c("loadings", "intercepts"), partial = wording)
  )
  if (length(loading_items)) {
    specs[[length(specs) + 1L]] <- list(
      name = "Partial metric", equal = "loadings",
      partial = c(wording, loadings))
  }
  if (length(intercept_items) || length(loading_items)) {
    specs[[length(specs) + 1L]] <- list(
      name = "Partial scalar", equal = c("loadings", "intercepts"),
      partial = c(wording, loadings, intercepts))
  }
  specs
}

# One CFA at a time. Cache each successful or failed attempt and distinguish
# optimiser runtime from fit-index calculations.
aot_fast_fit <- function(ctx, spec, max_iterations, cache_dir) {
  key <- digest::digest(list(ctx$signature, aot_model_syntax(aot_items),
    spec, max_iterations, R.version.string,
    as.character(utils::packageVersion("lavaan"))))
  prefix <- file.path(cache_dir, paste0(gsub(" ", "_", spec$name), "-", key))
  fit_path <- paste0(prefix, "-model.rds")
  summary_path <- paste0(prefix, "-summary.rds")
  if (file.exists(summary_path)) {
    message("CACHED: ", spec$name)
    return(readRDS(summary_path))
  }

  message(format(Sys.time(), "%H:%M:%S"), " START ", spec$name,
    " [", paste(spec$partial, collapse = ", "), "]")
  t_fit <- system.time({
    caught <- aot_capture(lavaan::cfa(
      model = aot_model_syntax(aot_items),
      data = ctx$data,
      group = "Country",
      group.label = ctx$labels,
      estimator = "MLR", std.lv = TRUE, missing = "listwise",
      group.equal = spec$equal, group.partial = spec$partial,
      # A bounded limit helps diagnose non-convergence rather than
      # allowing individual problematic models to run for hours.
      control = list(iter.max = as.integer(max_iterations),
                     eval.max = as.integer(2L * max_iterations))
    ))
  })[["elapsed"]]
  message(format(Sys.time(), "%H:%M:%S"), " OPTIMISER FINISHED ",
          spec$name, " (", round(t_fit, 1), " s)")

  fit <- caught$value
  t_summary <- system.time({
    row <- aot_fit_row(fit, spec$name, caught$warnings)
  })[["elapsed"]]
  row$optimiser_seconds <- t_fit
  row$diagnostics_seconds <- t_summary
  row$model_path <- fit_path
  row$equal <- paste(spec$equal, collapse = ", ")
  row$partial <- paste(spec$partial, collapse = ", ")

  if (inherits(fit, "lavaan")) saveRDS(fit, fit_path)
  saveRDS(row, summary_path)
  message(format(Sys.time(), "%H:%M:%S"), " DONE ", spec$name,
    " CFI=", round(row$CFI, 4),
    " converged=", row$converged,
    " admissible=", row$admissible,
    " diagnostics=", round(t_summary, 1), " s")
  row
}

run_aot_partial_fast <- function(loading_items = character(),
                                 intercept_items = NULL,
                                 max_iterations = 1000L,
                                 cache_dir = "AOT_partial_fast_cache") {
  stopifnot(requireNamespace("lavaan", quietly = TRUE),
            requireNamespace("digest", quietly = TRUE))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  message(format(Sys.time(), "%H:%M:%S"), " Checking cached input data and legacy fits")
  ctx <- aot_context()
  message(format(Sys.time(), "%H:%M:%S"), " Input validation complete")

  # Descriptive first-pass screen using the *previously fitted* scalar model.
  # Large country-specific mean residuals suggest possible intercept DIF, but
  # do not constitute a formal constraint test.
  message(format(Sys.time(), "%H:%M:%S"), " Computing legacy item diagnostics")
  diagnostics <- aot_item_diagnostics(ctx)
  diagnostics <- diagnostics[order(-diagnostics$scalar_mean_residual_RMS), ]
  write.csv(diagnostics,
    file.path(cache_dir, "scalar-residual-screen.csv"), row.names = FALSE)
  message("Items ranked by legacy scalar-model mean residual RMS:")
  print(diagnostics[, c("item", "scalar_mean_residual_RMS")],
        row.names = FALSE)

  # Default is a SMALL exploratory analysis, not a 13-item deletion search.
  if (is.null(intercept_items)) {
    intercept_items <- head(diagnostics$item, 2L)
    message("Exploratory intercept candidates (please verify): ",
            paste(intercept_items, collapse = ", "))
  }
  loading_items <- unique(as.character(loading_items))
  intercept_items <- unique(as.character(intercept_items))
  if (!all(c(loading_items, intercept_items) %in% aot_items)) {
    stop("Unknown item name in loading_items or intercept_items")
  }
  # No need to hold the legacy fitted models in RAM during the new fits.
  ctx$legacy <- NULL
  gc(verbose = FALSE)

  specs <- aot_fast_specs(loading_items, intercept_items)
  rows <- do.call(rbind, lapply(specs, function(s) {
    aot_fast_fit(ctx, s, max_iterations, cache_dir)
  }))
  rownames(rows) <- NULL

  # Compare each partial model against the right less-constrained baseline:
  # partial scalar vs partial metric if substantive loadings were released.
  comparison <- data.frame(
    comparison = c("Metric - Configural", "Scalar - Metric",
      if (length(loading_items)) "Partial metric - Configural",
      if ("Partial scalar" %in% rows$stage) {
        if (length(loading_items)) "Partial scalar - Partial metric"
        else "Partial scalar - Metric"
      }),
    delta_CFI = NA_real_, delta_RMSEA = NA_real_, delta_SRMR = NA_real_
  )
  pairs <- list(c("Configural", "Metric"), c("Metric", "Scalar"))
  if (length(loading_items)) pairs[[length(pairs) + 1L]] <-
    c("Configural", "Partial metric")
  if ("Partial scalar" %in% rows$stage) pairs[[length(pairs) + 1L]] <-
    c(if (length(loading_items)) "Partial metric" else "Metric", "Partial scalar")
  for (i in seq_along(pairs)) {
    ix <- match(pairs[[i]], rows$stage)
    comparison$delta_CFI[i] <- diff(rows$CFI[ix])
    comparison$delta_RMSEA[i] <- diff(rows$RMSEA[ix])
    comparison$delta_SRMR[i] <- diff(rows$SRMR[ix])
  }

  write.csv(rows, file.path(cache_dir, "model-fit-summary.csv"), row.names = FALSE)
  write.csv(comparison, file.path(cache_dir, "model-comparisons.csv"), row.names = FALSE)
  message("Completed ", nrow(rows), " CFA models (instead of 54).")
  list(diagnostics = diagnostics, fits = rows, comparison = comparison,
       loading_items = loading_items, intercept_items = intercept_items)
}
