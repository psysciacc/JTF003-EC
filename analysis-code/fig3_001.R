# Alternative main-text Figure 3 drafts for JTF003-EC.
# Run from the analysis-code directory with: Rscript fig3_001.R
# Reads the cached models/data created by results.qmd and writes PNG and PDF versions of both options:
#   fig3_001_option_A.png  (H1, H2 and H3 estimates)
#   fig3_001_option_B.png  (H1 and correction probabilities by condition)
#
# Country intervals use the fitted country's conditional random effects and
# their conditional variances, added to fixed-effect uncertainty. These are
# approximate Wald intervals (the fixed/random covariance is omitted). Group
# and overall intervals use the fixed-effect covariance matrix.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(lme4)
  library(countrycode)
  library(patchwork)
})

cache_dir <- "results_cache"
html_cache <- file.path(cache_dir, "knitr-html")
output_dir <- Sys.getenv("FIG3_OUTPUT_DIR", unset = ".")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(cache_dir)) {
  stop("Run results.qmd once first; results_cache/ was not found.")
}

load_cached_chunk <- function(label) {
  files <- list.files(html_cache,
                      pattern = paste0("^", label, "_[[:xdigit:]]+\\.RData$"),
                      full.names = TRUE)
  if (length(files) != 1L) {
    stop("Expected one cached chunk for '", label, "' in ", html_cache,
         "; found ", length(files), ". Compile results.qmd to refresh its cache.")
  }
  env <- new.env(parent = emptyenv())
  load(files, envir = env)
  env
}

models_h1 <- load_cached_chunk("models_H1")
models_h23 <- load_cached_chunk("v2-chunk-016")
prepared <- load_cached_chunk("data_preparation_2")

required_objects <- c("mH1_all", "mH1_i", "mH1_c")
if (!all(vapply(required_objects, exists, logical(1), envir = models_h1,
                inherits = FALSE))) {
  stop("The cached H1 model chunk does not contain all expected models.")
}
required_objects <- c("m_all", "m_i", "m_c")
if (!all(vapply(required_objects, exists, logical(1), envir = models_h23,
                inherits = FALSE))) {
  stop("The cached H2/H3 model chunk does not contain all expected models.")
}
required_objects <- c("d", "dat_H1", "dat_Hother")
if (!all(vapply(required_objects, exists, logical(1), envir = prepared,
                inherits = FALSE))) {
  stop("The cached data-preparation chunk is incomplete.")
}

mH1_all <- models_h1$mH1_all
mH1_i <- models_h1$mH1_i
mH1_c <- models_h1$mH1_c
m_all <- models_h23$m_all
m_i <- models_h23$m_i
m_c <- models_h23$m_c
d <- prepared$d
dat_H1 <- prepared$dat_H1
dat_Hother <- prepared$dat_Hother

COLORS <- c(
  "Individualistic" = "#56B4E9",
  "Overall" = "#555555",
  "Collectivistic" = "#E69F00"
)
GROUP_LABELS <- c(
  "Individualistic" = "More individualistic countries",
  "Overall" = "Overall sample",
  "Collectivistic" = "More collectivistic countries"
)
GROUP_X <- c("Individualistic" = 1, "Overall" = 2, "Collectivistic" = 3)
CONDITION_LABELS <- c("Control", "Feedback", "Justification", "Both")

# Study country names are translated to ISO2 codes for FlagCDN. countrycode's
# built-in aliases cover the study countries; custom mappings make these
# historically variable names explicit.
country_aliases <- c(
  "Kosovo" = "XK", "Hong Kong" = "HK", "North Macedonia" = "MK",
  "Taiwan" = "TW"
)
country_lookup <- d %>%
  distinct(Country, culture) %>%
  filter(culture %in% c("Individualistic", "Collectivistic")) %>%
  mutate(
    ISO2 = countrycode::countrycode(
      Country, origin = "country.name", destination = "iso2c",
      custom_match = country_aliases, warn = FALSE
    ),
    Group = as.character(culture),
    flag_url = if_else(
      !is.na(ISO2),
      paste0("https://flagcdn.com/w40/", tolower(ISO2), ".png"),
      NA_character_
    )
  ) %>%
  select(Country, ISO2, Group, flag_url)

# Optional flag-image layer. If ggimage is unavailable, the plots remain
# self-contained and show small coloured country points instead.
use_flag_images <- !identical(Sys.getenv("FIG3_FLAGS", "true"), "false") &&
  requireNamespace("ggimage", quietly = TRUE)
if (use_flag_images) {
  flag_dir <- Sys.getenv("FIG3_FLAG_DIR", file.path(tempdir(), "jtf003-country-flags"))
  dir.create(flag_dir, showWarnings = FALSE, recursive = TRUE)
  flag_lookup <- country_lookup %>%
    filter(!is.na(ISO2)) %>%
    distinct(ISO2, flag_url) %>%
    mutate(flag_file = file.path(flag_dir, paste0(tolower(ISO2), ".png")))
  for (i in seq_len(nrow(flag_lookup))) {
    if (!file.exists(flag_lookup$flag_file[i])) {
      try(
        utils::download.file(flag_lookup$flag_url[i], flag_lookup$flag_file[i],
                             mode = "wb", quiet = TRUE),
        silent = TRUE
      )
    }
  }
  valid_flag <- vapply(flag_lookup$flag_file, function(path) {
    file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
  }, logical(1))
  flag_lookup$flag_file[!valid_flag] <- NA_character_
  country_lookup <- country_lookup %>%
    left_join(flag_lookup %>% select(ISO2, flag_file), by = "ISO2")
  if (!any(!is.na(country_lookup$flag_file))) {
    use_flag_images <- FALSE
  }
}
if (!use_flag_images) {
  message("Country marks will use small coloured points. Install ggimage and allow FlagCDN downloads to use flag images.")
}

# A model-matrix row in the coefficient order used by the four fixed terms.
design_row <- function(model, feedback = 0, justification = 0) {
  beta <- lme4::fixef(model)
  row <- c(
    "(Intercept)" = 1,
    "feedback" = feedback,
    "justification" = justification,
    "feedback:justification" = feedback * justification
  )
  row[names(beta)]
}

fixed_linear_estimate <- function(model, x) {
  beta <- lme4::fixef(model)
  V <- as.matrix(vcov(model))
  x <- x[names(beta)]
  estimate <- sum(x * beta)
  se <- sqrt(drop(t(x) %*% V %*% x))
  c(estimate = estimate, se = se)
}

# Country-specific conditional prediction from a random-intercept/random-slope
# model. The || specification in results.qmd estimates each country term
# independently in the prior. The interval approximation below omits
# posterior cross-term and fixed/random covariances; it is descriptive only.
country_linear_estimate <- function(model, x, country) {
  beta <- lme4::fixef(model)
  random <- lme4::ranef(model, condVar = TRUE)$Country
  post_var <- attr(random, "postVar")
  beta_names <- names(beta)
  x <- x[beta_names]
  ii <- match(country, rownames(random))
  if (is.na(ii)) return(c(estimate = NA_real_, se = NA_real_))

  country_beta <- beta + unlist(random[ii, beta_names, drop = TRUE])
  fixed_var <- drop(t(x) %*% as.matrix(vcov(model)) %*% x)
  random_var <- sum(vapply(beta_names, function(term) {
    term_post <- post_var[[term]]
    if (is.null(term_post)) return(0)
    x[[term]]^2 * term_post[1, 1, ii]
  }, numeric(1)))
  c(estimate = sum(x * country_beta),
    se = sqrt(max(0, fixed_var + random_var)))
}

# Keep the link-scale centre distinct from the output column name: tibble
# evaluates columns sequentially, so reusing 'estimate' would transform the
# centre twice when computing the limits.
to_probability <- function(estimate, se) {
  eta <- estimate
  tibble(estimate = plogis(eta), lower = plogis(eta - qnorm(.975) * se),
         upper = plogis(eta + qnorm(.975) * se))
}
to_odds_ratio <- function(estimate, se) {
  eta <- estimate
  tibble(estimate = exp(eta), lower = exp(eta - qnorm(.975) * se),
         upper = exp(eta + qnorm(.975) * se))
}

pooled_row <- function(model, group, feedback = 0, justification = 0,
                       measure = c("probability", "feedback_or", "justification_or")) {
  measure <- match.arg(measure)
  x <- switch(
    measure,
    probability = design_row(model, feedback, justification),
    feedback_or = design_row(model, feedback = 1, justification = 0) -
      design_row(model, feedback = 0, justification = 0),
    justification_or = design_row(model, feedback = 0, justification = 1) -
      design_row(model, feedback = 0, justification = 0)
  )
  est <- fixed_linear_estimate(model, x)
  vals <- if (measure == "probability") {
    to_probability(est[["estimate"]], est[["se"]])
  } else {
    to_odds_ratio(est[["estimate"]], est[["se"]])
  }
  bind_cols(tibble(Group = group), vals)
}

country_row <- function(model, country, group, feedback = 0, justification = 0,
                        measure = c("probability", "feedback_or", "justification_or")) {
  measure <- match.arg(measure)
  x <- switch(
    measure,
    probability = design_row(model, feedback, justification),
    feedback_or = design_row(model, feedback = 1, justification = 0) -
      design_row(model, feedback = 0, justification = 0),
    justification_or = design_row(model, feedback = 0, justification = 1) -
      design_row(model, feedback = 0, justification = 0)
  )
  est <- country_linear_estimate(model, x, country)
  vals <- if (measure == "probability") {
    to_probability(est[["estimate"]], est[["se"]])
  } else {
    to_odds_ratio(est[["estimate"]], est[["se"]])
  }
  bind_cols(tibble(Country = country, Group = group), vals)
}

# Country positions are repeatable and spread within the two cultural-group
# slots. The pooled estimates remain exactly centred on their group positions.
position_country_marks <- function(x, base_positions, jitter_width = 0.32) {
  x %>%
    group_by(Group) %>%
    arrange(Country, .by_group = TRUE) %>%
    mutate(
      x = unname(base_positions[Group]) +
        if (n() == 1L) 0 else seq(-jitter_width, jitter_width, length.out = n())
    ) %>%
    ungroup()
}

# Country-level rows for all panels, with the culture and ISO2/flag metadata.
country_models <- country_lookup %>%
  filter(Country %in% rownames(lme4::ranef(m_all)$Country))

# -------------------------------------------------------------------------
# Option A: H1 control probability; H2 feedback OR; H3 justification OR.
# -------------------------------------------------------------------------
h1_country <- map_dfr(seq_len(nrow(country_models)), function(i) {
  row <- country_models[i, ]
  country_row(mH1_all, row$Country, row$Group,
              feedback = 0, justification = 0, measure = "probability")
}) %>%
  left_join(country_models, by = c("Country", "Group")) %>%
  position_country_marks(GROUP_X) %>%
  mutate(Hypothesis = "H1")

make_or_country <- function(term) {
  map_dfr(seq_len(nrow(country_models)), function(i) {
    row <- country_models[i, ]
    country_row(m_all, row$Country, row$Group, measure = term)
  }) %>%
    left_join(country_models, by = c("Country", "Group")) %>%
    position_country_marks(GROUP_X) %>%
    mutate(Hypothesis = if_else(term == "feedback_or", "H2", "H3"))
}

country_A <- bind_rows(
  h1_country,
  make_or_country("feedback_or"),
  make_or_country("justification_or")
)

pooled_A <- bind_rows(
  pooled_row(mH1_i, "Individualistic", measure = "probability") %>% mutate(Hypothesis = "H1"),
  pooled_row(mH1_all, "Overall", measure = "probability") %>% mutate(Hypothesis = "H1"),
  pooled_row(mH1_c, "Collectivistic", measure = "probability") %>% mutate(Hypothesis = "H1"),
  pooled_row(m_i, "Individualistic", measure = "feedback_or") %>% mutate(Hypothesis = "H2"),
  pooled_row(m_all, "Overall", measure = "feedback_or") %>% mutate(Hypothesis = "H2"),
  pooled_row(m_c, "Collectivistic", measure = "feedback_or") %>% mutate(Hypothesis = "H2"),
  pooled_row(m_i, "Individualistic", measure = "justification_or") %>% mutate(Hypothesis = "H3"),
  pooled_row(m_all, "Overall", measure = "justification_or") %>% mutate(Hypothesis = "H3"),
  pooled_row(m_c, "Collectivistic", measure = "justification_or") %>% mutate(Hypothesis = "H3")
) %>%
  mutate(x = unname(GROUP_X[Group]))

country_A <- country_A %>%
  mutate(Group = factor(Group, levels = names(COLORS)))

condition_grid <- tibble(
  Condition = factor(CONDITION_LABELS, levels = CONDITION_LABELS),
  feedback = c(0, 1, 0, 1),
  justification = c(0, 0, 1, 1),
  condition_x = 1:4
)

make_model_condition_predictions <- function(model, group) {
  map_dfr(seq_len(nrow(condition_grid)), function(i) {
    cond <- condition_grid[i, ]
    x <- design_row(model, cond$feedback, cond$justification)
    est <- fixed_linear_estimate(model, x)
    bind_cols(
      tibble(Group = group, Condition = cond$Condition,
             condition_x = cond$condition_x),
      to_probability(est[["estimate"]], est[["se"]])
    )
  })
}

pooled_B <- bind_rows(
  make_model_condition_predictions(m_i, "Individualistic"),
  make_model_condition_predictions(m_all, "Overall"),
  make_model_condition_predictions(m_c, "Collectivistic")
) %>%
  mutate(
    Group = factor(Group, levels = names(COLORS)),
    offset = c(Individualistic = -0.23, Overall = 0,
               Collectivistic = 0.23)[as.character(Group)],
    x = condition_x + offset
  )

country_B <- map_dfr(seq_len(nrow(country_models)), function(i) {
  country <- country_models$Country[i]
  group <- country_models$Group[i]
  map_dfr(seq_len(nrow(condition_grid)), function(j) {
    cond <- condition_grid[j, ]
    country_row(m_all, country, group,
                feedback = cond$feedback, justification = cond$justification,
                measure = "probability") %>%
      mutate(Condition = cond$Condition, condition_x = cond$condition_x)
  })
}) %>%
  left_join(country_lookup, by = c("Country", "Group")) %>%
  group_by(Group, Condition) %>%
  arrange(Country, .by_group = TRUE) %>%
  mutate(
    group_offset = c(Individualistic = -0.23, Overall = 0,
                     Collectivistic = 0.23)[as.character(Group)],
    country_jitter = if (n() == 1L) 0 else seq(-0.12, 0.12, length.out = n()),
    x = condition_x + group_offset + country_jitter
  ) %>%
  ungroup()

# Validate EVERY displayed interval, including those for country overlays.
validate_intervals <- function(dat, name) {
  if (any(!is.finite(as.matrix(dat[c("estimate", "lower", "upper")]))))
    stop(name, ": missing or non-finite estimates/limits.")
  if (any(dat$lower > dat$estimate | dat$estimate > dat$upper))
    stop(name, ": a point estimate lies outside its confidence interval.")
}
for (nm in c("pooled_A", "country_A", "pooled_B", "country_B"))
  validate_intervals(get(nm), nm)
# Independent check against the original cached model's link-scale interval.
z <- qnorm(.975)
b <- fixef(m_all)["feedback"]
se <- sqrt(vcov(m_all)["feedback", "feedback"])
stopifnot(isTRUE(all.equal(
  as.numeric(unlist(filter(pooled_A, Group == "Overall", Hypothesis == "H2")[
    c("estimate", "lower", "upper")])),
  unname(exp(c(b, b-z*se, b+z*se))), tolerance = 1e-12)))
b <- fixef(mH1_all)["(Intercept)"]
se <- sqrt(vcov(mH1_all)["(Intercept)", "(Intercept)"])
stopifnot(isTRUE(all.equal(
  as.numeric(unlist(filter(pooled_A, Group == "Overall", Hypothesis == "H1")[
    c("estimate", "lower", "upper")])),
  unname(plogis(c(b, b-z*se, b+z*se))), tolerance = 1e-12)))

# Exact theme definition from results.qmd; no report execution is needed.
ml_theme <- function(base_size = 9, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(text = element_text(family = base_family, size = base_size),
          panel.border = element_blank(),
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = rel(0.8)),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 0.4),
          axis.text.x = element_text(size = 7, colour = "black"),
          axis.text.y = element_text(size = 7, colour = "black"),
          axis.ticks = element_line(colour = "black"))
}
figure_theme <- function(report_style = FALSE) {
  base <- if (report_style) ml_theme() else theme_minimal(base_size = 10)
  base + theme(legend.position = "bottom", legend.title = element_blank(),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey35"),
    plot.tag = element_text(face = "bold"),
    plot.margin = margin(8, 10, 8, 8), axis.title.x = element_blank())
}
figure_scales <- function() {
  list(scale_colour_manual(values = COLORS, limits = names(COLORS), guide = "none"),
       scale_fill_manual(values = COLORS, limits = names(COLORS),
                         labels = GROUP_LABELS, name = NULL),
       scale_shape_manual(values = c(Individualistic = 21, Overall = 23,
                                      Collectivistic = 21), guide = "none"),
       guides(fill = guide_legend(override.aes = list(shape = c(21,23,21),
                         colour = "white", size = 3.5, alpha = 1))))
}
country_layers <- function(dat, flags = TRUE, image_size = .028) {
  layers <- list(geom_linerange(data = dat,
    aes(x = x, ymin = lower, ymax = upper, colour = Group),
    inherit.aes = FALSE, linewidth = .25, alpha = .20, show.legend = FALSE))
  flag_rows <- rep(FALSE, nrow(dat))
  if (flags && use_flag_images && "flag_file" %in% names(dat)) {
    flag_rows <- !is.na(dat$flag_file)
    if (any(flag_rows)) layers <- c(layers, list(ggimage::geom_image(
      data = dat[flag_rows, ], aes(x = x, y = estimate, image = flag_file),
      inherit.aes = FALSE, size = image_size,
      # ggimage only applies its alpha argument when a colour tint is given.
      # Blend flags with white to keep their colours subdued on these white panels.
      image_fun = function(img) magick::image_colorize(img, 45, "white"))))
  }
  if (any(!flag_rows)) layers <- c(layers, list(geom_point(
    data = dat[!flag_rows, ], aes(x = x, y = estimate, colour = Group),
    inherit.aes = FALSE, size = 1.3, alpha = .45, show.legend = FALSE)))
  layers
}
summary_layers <- function(dat) {
  list(geom_errorbar(data = dat,
         aes(x = x, ymin = lower, ymax = upper, colour = Group),
         inherit.aes = FALSE, width = .07, linewidth = .8, show.legend = FALSE),
       geom_point(data = dat, aes(x = x, y = estimate, fill = Group, shape = Group),
         inherit.aes = FALSE, size = 3.4, colour = "white", stroke = .65))
}
# Both OR panels share limits that include ALL country and pooled intervals.
OR_LIMITS <- range(c(1, unlist(filter(country_A, Hypothesis != "H1")[c("lower","upper")]),
                          unlist(filter(pooled_A, Hypothesis != "H1")[c("lower","upper")])))
OR_LIMITS <- exp(log(OR_LIMITS) + c(-.10, .10))
H1_LABEL <- "Pr. first correct in reflective phase\n(among trials correct in either phase)"

plot_A_panel <- function(hypothesis, flags = TRUE, report_style = FALSE) {
  dat <- filter(country_A, Hypothesis == hypothesis)
  pooled <- filter(pooled_A, Hypothesis == hypothesis)
  p <- ggplot() +
    geom_hline(yintercept = if (hypothesis == "H1") .5 else 1,
               linetype = "dashed", linewidth = .4, colour = "grey55") +
    country_layers(dat, flags) + summary_layers(pooled) + figure_scales() +
    scale_x_continuous(breaks = 1:3,
      labels = c("Individualistic", "Overall", "Collectivistic"), limits = c(.45, 3.55)) +
    labs(title = switch(hypothesis, H1 = "H1 · Control condition", H2 = "H2 · Feedback",
                        H3 = "H3 · Justification"),
         subtitle = switch(hypothesis, H1 = "No feedback or justification",
                           H2 = "Compared with control; no justification",
                           H3 = "Compared with control; no feedback"),
         y = if (hypothesis == "H1") H1_LABEL else "Odds ratio (log scale)") +
    figure_theme(report_style)
  if (hypothesis == "H1") p + scale_y_continuous(limits = c(0,1),
      breaks = seq(0,1,.2), labels = scales::label_percent(accuracy = 1)) else
    p + scale_y_log10(limits = OR_LIMITS, breaks = c(.5,.75,1,1.25,1.5,2,2.5,3,4),
                      labels = scales::label_number())
}
plot_correction <- function(flags = TRUE, report_style = FALSE, facets = FALSE) {
  cc <- country_B
  pp <- pooled_B
  if (facets) {
    cc$x <- unname(GROUP_X[as.character(cc$Group)]) + cc$country_jitter * 2.3
    pp$x <- unname(GROUP_X[as.character(pp$Group)])
  }
  p <- ggplot() + country_layers(cc, flags, image_size = if (facets) .035 else .021) +
    summary_layers(pp) + figure_scales() +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2),
                       labels = scales::label_percent(accuracy = 1)) +
    labs(title = "H2/H3 · Error correction", subtitle = "Initially incorrect intuitive responses",
         y = "Probability of correcting an error") + figure_theme(report_style)
  if (facets) p + facet_wrap(~Condition, nrow = 1) +
    scale_x_continuous(breaks = 1:3, labels = c("Indiv.", "Overall", "Collect."),
                       limits = c(.45,3.55)) else
    p + scale_x_continuous(breaks = 1:4, labels = CONDITION_LABELS, limits = c(.5,4.5))
}
make_figures <- function(flags = TRUE, report_style = FALSE, facets = FALSE) {
  a <- wrap_plots(lapply(c("H1","H2","H3"), plot_A_panel,
                         flags = flags, report_style = report_style),
                 nrow = 1, guides = "collect") + plot_annotation(tag_levels = "a")
  b <- wrap_plots(plot_A_panel("H1", flags, report_style),
                 plot_correction(flags, report_style, facets), nrow = 1,
                 widths = if (facets) c(1, 2.8) else c(1, 1.8), guides = "collect") +
    plot_annotation(tag_levels = "a")
  list(A = a & theme(legend.position = "bottom"),
       B = b & theme(legend.position = "bottom"))
}
save_figures <- function(figures, prefix, report_style = FALSE) {
  for (version in names(figures)) {
    stem <- file.path(output_dir, paste0(prefix, "_option_", version))
    width <- if (report_style && version == "B") 14 else 12
    ggsave(paste0(stem, ".png"), figures[[version]], width = width, height = 4.8,
           dpi = 400, bg = "white")
    ggsave(paste0(stem, ".pdf"), figures[[version]], width = width, height = 4.8,
           device = grDevices::cairo_pdf, bg = "white")
    message("Saved ", stem, ".png and .pdf")
  }
}
# The second script sources this file in data-only mode, reusing the calculations
# and plotting helpers without writing the first script's figures.
if (!isTRUE(getOption("fig3.data_only", FALSE)))
  save_figures(make_figures(), "fig3_001")
