# Exploratory AOT item-reduction and partial-invariance analysis.
# All writes are confined to AOT_partial_invariance_cache/.
# Run: Rscript -e 'source("AOT_partial_invariance.R"); run_aot_invariance()'

aot_cache_dir <- "AOT_partial_invariance_cache"
aot_items <- paste0("AOT_", 1:13)
aot_negative <- paste0("AOT_", c(2, 4, 5, 6, 7, 8, 10, 13))

aot_capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = unique(warnings))
}

aot_fit_row <- function(fit, stage, warnings = character(), error = "") {
  measures <- setNames(rep(NA_real_, 7),
    c("cfi.robust", "tli.robust", "rmsea.robust", "srmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))
  converged <- admissible <- finite_se <- negative_scale_fixed <- FALSE
  N <- groups <- NA_integer_
  if (inherits(fit, "lavaan")) {
    converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
    check <- aot_capture(lavaan::lavInspect(fit, "post.check"))
    admissible <- isTRUE(check$value)
    warnings <- unique(c(warnings, check$warnings))
    N <- sum(lavaan::lavInspect(fit, "nobs"))
    groups <- lavaan::lavInspect(fit, "ngroups")
    if (converged) {
      pe <- lavaan::parameterTable(fit)
      nv <- pe[pe$lhs == "NEG" & pe$op == "~~" & pe$rhs == "NEG", ]
      negative_scale_fixed <- nrow(nv) == groups && all(nv$free == 0L) && all(abs(nv$est - 1) < 1e-10)
      finite_se <- all(is.finite(pe$se[pe$free > 0L]))
      fm <- aot_capture(lavaan::fitMeasures(fit, names(measures)))
      warnings <- unique(c(warnings, fm$warnings))
      if (!inherits(fm$value, "error")) measures[] <- fm$value
      else error <- conditionMessage(fm$value)
    }
  }
  if (inherits(fit, "error")) error <- conditionMessage(fit)
  data.frame(stage = stage, N = N, groups = groups,
    converged = converged, admissible = admissible, finite_se = finite_se,
    negative_scale_fixed = negative_scale_fixed,
    CFI = measures[1], TLI = measures[2], RMSEA = measures[3], SRMR = measures[4],
    chisq_scaled = measures[5], df_scaled = measures[6], p_scaled = measures[7],
    warning = paste(warnings, collapse = " | "), error = error,
    row.names = NULL, stringsAsFactors = FALSE)
}

aot_context <- function() {
  dir.create(file.path(aot_cache_dir, "models"), recursive = TRUE, showWarnings = FALSE)
  inputs <- c("results_v2_cache/derived/prepared-data.rds", "AOT_measurement_invariance.rds")
  provenance <- data.frame(path = inputs, md5 = unname(tools::md5sum(inputs)))
  d <- readRDS(inputs[1])$d
  countries <- names(which(table(d$Country) >= 100L))
  d <- d[d$Country %in% countries & complete.cases(d[, aot_items]),
         c("participant", "Country", aot_items)]
  legacy <- readRDS(inputs[2])
  base <- legacy$fit_AOT_config
  labels <- lavaan::lavInspect(base, "group.label")
  if (!setequal(unique(d$Country), labels)) stop("Baseline country selection differs from current data.")
  # Compare every scored response, not just sample sizes or covariance matrices.
  identical_data <- vapply(seq_along(labels), function(g) {
    isTRUE(all.equal(unname(as.matrix(d[d$Country == labels[g], aot_items])),
                    unname(base@Data@X[[g]]), check.attributes = FALSE))
  }, logical(1))
  if (!all(identical_data)) stop("Saved baseline does not match current scored AOT responses.")
  d$Country <- factor(d$Country, levels = labels)
  data_signature <- digest::digest(d)
  saveRDS(d, file.path(aot_cache_dir, "analysis-data.rds"))
  write.csv(provenance, file.path(aot_cache_dir, "input-provenance.csv"), row.names = FALSE)
  list(data = d, signature = data_signature, labels = labels, legacy = legacy,
       provenance = provenance, identical_data = identical_data)
}

aot_legacy_scaling_audit <- function(ctx) {
  do.call(rbind, lapply(c("fit_AOT_metric", "fit_AOT_scalar"), function(name) {
    m <- ctx$legacy[[name]]
    p <- lavaan::parameterTable(m)
    loads <- p[p$lhs == "NEG" & p$op == "=~", ]
    vars <- p[p$lhs == "NEG" & p$op == "~~" & p$rhs == "NEG", ]
    e <- lavaan::lavInspect(m, "est")[[2]]
    L <- e$lambda; P <- e$psi; A <- e$alpha
    covariance <- L %*% P %*% t(L)
    mean <- L %*% A
    L[, "NEG"] <- L[, "NEG"] / 2
    P["NEG", "NEG"] <- P["NEG", "NEG"] * 4
    A["NEG", 1] <- A["NEG", 1] * 2
    data.frame(model = name, free_wording_loadings = sum(loads$free > 0L),
      wording_loading_equalities = sum(nzchar(loads$label)),
      free_wording_variances = sum(vars$free > 0L),
      max_covariance_change = max(abs(covariance - L %*% P %*% t(L))),
      max_mean_change = max(abs(mean - L %*% A)))
  }))
}

aot_model_syntax <- function(items) {
  negative <- intersect(items, aot_negative)
  if (length(negative) < 3L || length(setdiff(items, aot_negative)) < 2L)
    stop("Too few indicators to retain the planned substantive/method-factor structure.")
  paste0("AOT =~ ", paste(items, collapse = " + "), "\nNEG =~ ",
         paste(negative, collapse = " + "), "\nAOT ~~ 0*NEG\nNEG ~~ 1*NEG")
}

aot_cached_fit <- function(ctx, items, stage, free_items = character()) {
  syntax <- aot_model_syntax(items)
  equal <- switch(stage, Configural = character(), Metric = "loadings",
                  Scalar = c("loadings", "intercepts"))
  partial <- paste0("NEG=~", intersect(items, aot_negative))
  if (stage != "Configural") partial <- c(partial, paste0("AOT=~", free_items))
  if (stage == "Scalar") partial <- c(partial, paste0(free_items, "~1"))
  # paste0(character(0), suffix) is not length zero in R; omit empty releases.
  if (!length(free_items)) partial <- paste0("NEG=~", intersect(items, aot_negative))
  spec <- list(syntax = syntax, equal = equal, partial = partial,
    estimator = "MLR", std.lv = TRUE, missing = "listwise", group.order = ctx$labels,
    control = list(iter.max = 20000L), analysis = "aot-invariance-identified-v2")
  key <- digest::digest(list(ctx$signature, spec, R.version.string,
                             as.character(utils::packageVersion("lavaan"))))
  path <- file.path(aot_cache_dir, "models", paste0(stage, "-", key, ".rds"))
  meta_path <- sub("\\.rds$", "-summary.rds", path)
  if (file.exists(path) && file.exists(meta_path)) return(readRDS(meta_path))
  cat(format(Sys.time(), "%H:%M:%S"), "fit", stage, "retain", length(items),
      "drop", paste(setdiff(aot_items, items), collapse = ","),
      "free", paste(free_items, collapse = ","), "\n")
  flush.console()
  start <- proc.time()[3]
  out <- aot_capture(lavaan::cfa(syntax, data = ctx$data, group = "Country",
    group.label = ctx$labels, estimator = "MLR", std.lv = TRUE, missing = "listwise",
    group.equal = equal, group.partial = partial, control = spec$control))
  fit <- out$value
  row <- aot_fit_row(fit, stage, out$warnings)
  row$seconds <- unname(proc.time()[3] - start)
  row$path <- path
  row$retained <- paste(items, collapse = ", ")
  row$released <- paste(free_items, collapse = ", ")
  # Store failures too, with their errors in the manifest rather than hiding them.
  saveRDS(list(fit = if (inherits(fit, "lavaan")) fit else NULL,
               specification = spec, summary = row,
               numerical_backend = extSoftVersion()["BLAS"]), path)
  saveRDS(row, meta_path)
  cat(format(Sys.time(), "%H:%M:%S"), "done", stage, "CFI", round(row$CFI, 4),
      "admissible", row$admissible, "seconds", round(row$seconds), "\n")
  flush.console()
  row
}

aot_baseline <- function(ctx) {
  if (!is.null(ctx$baseline_rows)) return(ctx$baseline_rows)
  stages <- c("Configural", "Metric", "Scalar")
  rows <- do.call(rbind, lapply(stages, function(stage) aot_cached_fit(ctx, aot_items, stage)))
  rows$candidate <- "Full 13 items"
  rows$kind <- "Baseline"
  rows$dropped <- ""
  rows$n_items <- 13L
  rows
}

aot_candidate <- function(ctx, dropped, label, kind = "Single deletion", free_items = character()) {
  items <- setdiff(aot_items, dropped)
  stages <- if (kind == "Partial constraints") c("Metric", "Scalar") else c("Configural", "Metric", "Scalar")
  rows <- lapply(stages, function(stage) aot_cached_fit(ctx, items, stage, free_items))
  rows <- do.call(rbind, rows)
  if (kind == "Partial constraints") {
    config <- aot_baseline(ctx)
    config <- config[config$stage == "Configural", ]
    rows <- rbind(config[, names(rows)], rows)
  }
  rows$released <- paste(free_items, collapse = ", ")
  rows$candidate <- label
  rows$kind <- kind
  rows$dropped <- paste(dropped, collapse = ", ")
  rows$n_items <- length(items)
  rows
}

aot_compare <- function(rows) {
  rows <- rows[match(c("Configural", "Metric", "Scalar"), rows$stage), ]
  rows$dCFI <- c(NA, diff(rows$CFI))
  rows$dRMSEA <- c(NA, diff(rows$RMSEA))
  rows$dSRMR <- c(NA, diff(rows$SRMR))
  rows
}

aot_rank <- function(rows) {
  candidates <- split(rows, rows$candidate)
  do.call(rbind, lapply(candidates, function(r) {
    r <- aot_compare(r)
    valid <- all(r$converged & r$admissible & r$finite_se & r$negative_scale_fixed) &&
      all(is.finite(r$CFI)) && all(r$error == "")
    data.frame(candidate = r$candidate[1], kind = r$kind[1], dropped = r$dropped[1],
      released = r$released[1], n_items = r$n_items[1], valid = valid,
      config_CFI = r$CFI[1], metric_CFI = r$CFI[2], scalar_CFI = r$CFI[3],
      metric_dCFI = r$dCFI[2], scalar_dCFI = r$dCFI[3],
      config_RMSEA = r$RMSEA[1], metric_dRMSEA = r$dRMSEA[2], scalar_dRMSEA = r$dRMSEA[3],
      metric_dSRMR = r$dSRMR[2], scalar_dSRMR = r$dSRMR[3],
      CFI_loss = if (valid) sum(pmax(0, -r$dCFI[2:3])) else Inf,
      within_01 = valid && all(r$dCFI[2:3] >= -.01),
      stringsAsFactors = FALSE, row.names = NULL)
  }))
}

aot_item_diagnostics <- function(ctx) {
  p <- lavaan::parameterEstimates(ctx$legacy$fit_AOT_config, standardized = TRUE)
  p <- p[p$op == "=~" & p$lhs == "AOT", ]
  residuals <- lavaan::lavResiduals(ctx$legacy$fit_AOT_scalar, zstat = FALSE, summary = FALSE)
  means <- do.call(rbind, lapply(residuals, function(r) r$mean[aot_items]))
  do.call(rbind, lapply(aot_items, function(item) {
    q <- p[p$rhs == item, "std.all"]
    data.frame(item = item, reversed = item %in% aot_negative,
      mean_loading = mean(q), sd_loading = sd(q), min_loading = min(q),
      max_loading = max(q), countries_negative_loading = sum(q < 0),
      scalar_mean_residual_RMS = sqrt(mean(means[, item]^2)))
  }))
}

run_aot_invariance <- function(workers = 2L) {
  if (!requireNamespace("lavaan", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE))
    stop("Install lavaan and digest before running this analysis.")
  ctx <- aot_context()
  signature <- digest::digest(list(ctx$signature, ctx$provenance,
    tools::md5sum("AOT_partial_invariance.R"), as.character(utils::packageVersion("lavaan")),
    R.version.string))
  result_path <- file.path(aot_cache_dir, "analysis-results.rds")
  if (file.exists(result_path)) {
    result <- readRDS(result_path)
    if (identical(result$signature, signature)) return(result)
  }
  scaling_audit <- aot_legacy_scaling_audit(ctx)
  legacy_rows <- do.call(rbind, lapply(seq_along(c("Configural", "Metric", "Scalar")), function(i) {
    m <- ctx$legacy[[c("fit_AOT_config", "fit_AOT_metric", "fit_AOT_scalar")[i]]]
    aot_fit_row(m, c("Configural", "Metric", "Scalar")[i])
  }))
  base <- aot_baseline(ctx)
  if (!all(base$converged & base$admissible & base$finite_se & base$negative_scale_fixed)) stop("Baseline is not admissible.")
  diagnostics <- aot_item_diagnostics(ctx)
  baseline_version <- ctx$legacy$fit_AOT_config@version
  ctx$baseline_rows <- base
  ctx$legacy <- NULL
  gc()
  write.csv(diagnostics, file.path(aot_cache_dir, "item-diagnostics.csv"), row.names = FALSE)
  # Exhaustive single-deletion screen. Small worker count limits peak memory.
  singles <- parallel::mclapply(aot_items, function(item) {
    aot_candidate(ctx, item, paste("Drop", item))
  }, mc.cores = if (.Platform$OS.type == "unix") workers else 1L, mc.preschedule = FALSE)
  if (any(vapply(singles, inherits, logical(1), "try-error"))) stop("A screening worker failed; rerun to resume caches.")
  singles <- do.call(rbind, singles)
  ranking <- aot_rank(singles)
  ranking <- ranking[order(!ranking$valid, ranking$CFI_loss, -ranking$config_CFI, ranking$dropped), ]
  write.csv(ranking, file.path(aot_cache_dir, "single-deletion-ranking.csv"), row.names = FALSE)
  ranked_items <- ranking$dropped[ranking$valid]
  if (length(ranked_items) < 3L) stop("Fewer than three admissible single-deletion candidates.")
  top <- head(ranked_items, 3L)
  cat("Top single-deletion candidates:", paste(top, collapse = ", "), "\n")
  # A bounded follow-up, not an exhaustive search over pairs or triples.
  followup_jobs <- c(lapply(2:3, function(n) list(type = "drop", n = n)),
                     lapply(1:3, function(n) list(type = "free", n = n)))
  followup <- parallel::mclapply(followup_jobs, function(job) {
    selected <- head(top, job$n)
    if (job$type == "drop") {
      aot_candidate(ctx, selected, paste0("Drop top ", job$n), kind = "Combined deletion")
    } else {
      aot_candidate(ctx, character(), paste0("Free top ", job$n),
        kind = "Partial constraints", free_items = selected)
    }
  }, mc.cores = if (.Platform$OS.type == "unix") workers else 1L, mc.preschedule = FALSE)
  if (any(vapply(followup, inherits, logical(1), "try-error"))) stop("A follow-up worker failed; rerun to resume caches.")
  all_rows <- rbind(base, singles, do.call(rbind, followup))
  all_rows <- do.call(rbind, lapply(split(all_rows, all_rows$candidate), aot_compare))
  summary <- aot_rank(all_rows)
  # Reliability and correspondence of the cumulative shortened raw mean scores.
  scores <- lapply(0:3, function(n) {
    dropped <- if (n == 0L) character() else head(top, n)
    items <- setdiff(aot_items, dropped)
    full <- rowMeans(ctx$data[, aot_items])
    reduced <- rowMeans(ctx$data[, items])
    alpha <- if (requireNamespace("psych", quietly = TRUE))
      aot_capture(psych::alpha(ctx$data[, items], check.keys = FALSE, warnings = FALSE))$value else NULL
    data.frame(removed = paste(dropped, collapse = ", "), n_items = length(items),
      alpha_pooled = if (is.list(alpha)) alpha$total$raw_alpha else NA_real_,
      correlation_with_full_score = cor(full, reduced), mean_score = mean(reduced), sd_score = sd(reduced))
  })
  result <- list(signature = signature, provenance = ctx$provenance, data_signature = ctx$signature,
    N = nrow(ctx$data), countries = length(ctx$labels),
    country_counts = as.data.frame(table(ctx$data$Country)),
    baseline_data_verified = all(ctx$identical_data), baseline_version = baseline_version,
    legacy_fits = legacy_rows, scaling_audit = scaling_audit,
    current_lavaan_version = as.character(utils::packageVersion("lavaan")),
    diagnostics = diagnostics, single_ranking = ranking, selected_items = top,
    fits = all_rows, comparison = summary, reliability = do.call(rbind, scores),
    session = utils::sessionInfo(), created = Sys.time())
  write.csv(all_rows, file.path(aot_cache_dir, "all-model-fit-statistics.csv"), row.names = FALSE)
  write.csv(summary, file.path(aot_cache_dir, "candidate-comparison.csv"), row.names = FALSE)
  saveRDS(result, result_path)
  result
}
