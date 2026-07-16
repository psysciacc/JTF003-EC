# Keep existing Country values.
# If Country is missing and fetch_id == 5, set Country = "NO" because all Norway data were collected in Norway.
# Else, if Country_Res is available, use country of residence.
# Else, use Q_Language only when that exact Q_Language maps to one unique country elsewhere in the dataset.


library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(countrycode)

# d <- d_no_test

# ---------------------------------------------------------------------------
# Conservative Country imputation from fetch_id, Country_Res, and Q_Language

norway_country_code <- "NO"  # change to "NOR" if that is your preferred code

clean_country_code <- function(x) {
  x <- as.character(x)
  x <- str_squish(str_to_upper(x))
  x <- na_if(x, "")
  x <- na_if(x, "NA")
  
  out <- map_chr(x, function(z) {
    if (is.na(z)) return(NA_character_)
    
    parts <- str_split(z, ",", simplify = FALSE)[[1]] %>%
      str_squish() %>%
      discard(~ .x == "")
    
    if (length(parts) == 0) {
      NA_character_
    } else if (length(unique(parts)) == 1) {
      unique(parts)
    } else {
      z
    }
  })
  
  out <- recode(
    out,
    "GB"  = "UK",
    "KU"  = "KW",  # Kuwait typo/non-standard code
    "ZB"  = "ZW",  # Zimbabwe typo/non-standard code
    "NOR" = "NO",
    .default = out
  )
  
  out
}

clean_country_name <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x <- na_if(x, "")
  x <- na_if(x, "NA")
  x <- na_if(x, "Prefer not to respond")
  x
}

country_name_to_code <- function(x) {
  x_clean <- clean_country_name(x)
  
  manual_lookup <- c(
    "United Kingdom" = "UK",
    "United States of America" = "US",
    "Hong Kong (SAR)" = "HK",
    "Russian Federation" = "RU",
    "North Macedonia" = "MK",
    "Türkiye" = "TR"
  )
  
  manual <- unname(manual_lookup[x_clean])
  
  cc <- countrycode(
    x_clean,
    origin = "country.name",
    destination = "iso2c",
    warn = FALSE
  )
  
  out <- if_else(!is.na(manual), manual, cc)
  clean_country_code(out)
}

is_finished_response <- function(x) {
  x %in% c(TRUE, "TRUE", "True", "true", "1", 1)
}

# ------------------------------------------------------------------------------

# store original Country value
if (!"Country_original" %in% names(d)) {
  d <- d %>%
    mutate(Country_original = Country)
}

d_country_base <- d %>%
  mutate(
    Country_clean = clean_country_code(Country),
    Country_Res_code = country_name_to_code(Country_Res),
    Country_Ori_code = country_name_to_code(Country_Ori)
  )

# Build a conservative Q_Language -> Country lookup:
# only use Q_Language values that map to exactly one observed country.
qlang_country_lookup <- d_country_base %>%
  filter(!is.na(Q_Language), !is.na(Country_clean)) %>%
  filter(str_detect(Country_clean, "^[A-Z]{2}$")) %>%
  distinct(Q_Language, Country_clean) %>%
  group_by(Q_Language) %>%
  summarise(
    n_countries_seen = n_distinct(Country_clean),
    countries_seen = str_c(sort(unique(Country_clean)), collapse = ", "),
    Country_from_QLanguage = if_else(
      n_countries_seen == 1,
      first(Country_clean),
      NA_character_
    ),
    .groups = "drop"
  )

# view language lookup:
# view(qlang_country_lookup)

### FINALLY generated proposed imputations

d_country_proposed <- d_country_base %>%
  left_join(
    qlang_country_lookup %>%
      select(
        Q_Language,
        Country_from_QLanguage,
        n_countries_seen,
        countries_seen
      ),
    by = "Q_Language"
  ) %>%
  mutate(
    Country_was_missing = is.na(Country_clean),
    
    Country_proposed = case_when(
      !is.na(Country_clean) ~ Country_clean,
      fetch_id == 5 ~ norway_country_code,
      !is.na(Country_from_QLanguage) ~ Country_from_QLanguage,
      TRUE ~ NA_character_
    ),
    
    Country_imputation_rule = case_when(
      !is.na(Country_clean) ~ "already present",
      fetch_id == 5 ~ "filled from fetch_id == 5: Norway offline data",
      !is.na(Country_from_QLanguage) ~ "filled from unique Q_Language mapping",
      TRUE ~ "still missing after conservative checks"
    ),
    
    Country = Country_proposed
  )

# ------------------------------------------------------------------------------
### Table to show original missing country cases
country_missing_review <- d_country_proposed %>%
  filter(Country_was_missing, is_finished_response(Finished)) %>%
  group_by(
    Q_Language,
    Country_Res_YN,
    Country_Res,
    fetch_id,
    Country_Res_code,
    Country_from_QLanguage,
    countries_seen,
    Country_proposed,
    Country_imputation_rule
  ) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(
    Country_imputation_rule,
    fetch_id,
    Q_Language,
    desc(N)
  )

print(country_missing_review, n = 200)
# view(country_missing_review)

### Export a smaller table showing only actual imputations
country_substitution_review <- country_missing_review %>%
  filter(!is.na(Country_proposed)) %>%
  arrange(fetch_id, Q_Language, Country_imputation_rule, desc(N))

print(country_substitution_review, n = 200)

write_csv(
  country_substitution_review,
  "country_substitutions_for_review.csv"
)

### also a table with unresolved ones
country_still_missing_review <- country_missing_review %>%
  filter(is.na(Country_proposed)) %>%
  arrange(fetch_id, Q_Language, desc(N))

print(country_still_missing_review, n = 200)

write_csv(
  country_still_missing_review,
  "country_still_missing_after_fill.csv"
)


### Save in final dataset
d <- d_country_proposed %>%
  select(
    -Country_clean,
    -Country_proposed,
    -Country_Res_code,
    -Country_Ori_code,
    -Country_from_QLanguage,
    -n_countries_seen,
    -countries_seen
  )

# export
write_csv(d, "all_withNorway_noTest_countryFill_raw.csv")

# ------------------------------------------------------------------------------
##### Summaries on how many countries were filled in, in total and broken down by rule
country_missing_summary <- d %>%
  summarise(
    total_responses = n(),
    country_missing_originally = sum(Country_was_missing, na.rm = TRUE),
    country_filled = sum(Country_was_missing & !is.na(Country), na.rm = TRUE),
    country_still_missing = sum(Country_was_missing & is.na(Country), na.rm = TRUE),
    pct_missing_filled = 100 * country_filled / country_missing_originally
  )

t(country_missing_summary)


