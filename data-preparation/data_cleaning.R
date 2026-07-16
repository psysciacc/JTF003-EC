# data cleaning

rm(list=ls())
library(tidyverse)

# source("fetch_data_v3.R")
d <- read_csv("all_withNorway_raw.csv")

unique(d$Test)
mean(is.na(d$Test))

# check if we find flagged ID
which(d$ResponseId=="R_2MiehJBgYVcN7Xj")

which(d$LabID=="1589") # 1589 AR-EG-1589
range(d$StartDate[which(d$LabID=="1589")])

# ---------------------------------------------------------------------------
# Correct known cases where real participants used the test link
# These should be treated as genuine data, so set Test to NA.

d <- d %>%
  mutate(
    # ensure date-only version
    # If RecordedDate is already POSIXct, as_date() works.
    # If it is character, parse_date_time() handles common Qualtrics formats.
    RecordedDate_date = case_when(
      inherits(RecordedDate, "POSIXt") ~ as_date(RecordedDate),
      TRUE ~ as_date(parse_date_time(
        RecordedDate,
        orders = c("ymd HMS", "ymd HM", "dmy HMS", "dmy HM", "mdy HMS", "mdy HM")
      ))
    ),
    
    # store original Test value for audit trail
    Test_original = Test,
    
    # Flag known real-data-via-test-link cases
    real_data_test_link = case_when(
      ResponseId == "R_2MiehJBgYVcN7Xj" ~ TRUE,
      LabID == "1589" & RecordedDate_date >= dmy("16/09/2025") ~ TRUE,
      LabID == "2143" ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Correct Test coding
    Test = if_else(real_data_test_link, NA_character_, as.character(Test))
  )

# this is what has been corrected
d %>%
  filter(real_data_test_link) %>%
  select(
    ResponseId,
    LabID,
    RecordedDate,
    RecordedDate_date,
    Test_original,
    Test,
    Country,
    Sample
  ) %>%
  arrange(LabID, RecordedDate_date)


# ---------------------------------------------------------------------------
# Mark known test responses that were submitted through the real link
# Applies only to the final online dataset: fetch_id == 4
# These should have Test = "Test" if recorded before real data collection started.

as_qualtrics_datetime <- function(x, tz = "Europe/London") {
  if (inherits(x, "POSIXt")) return(x)
  if (inherits(x, "Date")) return(as.POSIXct(x, tz = tz))

  parse_date_time(
    as.character(x),
    orders = c("ymd HMS", "ymd HM", "ymd",
               "dmy HMS", "dmy HM", "dmy",
               "mdy HMS", "mdy HM", "mdy"),
    tz = tz
  )
}

# # Use the timezone of RecordedDate if already set; otherwise choose the timezone
# # used by the Qualtrics export/account.
# qualtrics_tz <- attr(d$RecordedDate, "tzone")
# if (is.null(qualtrics_tz) || length(qualtrics_tz) == 0 || is.na(qualtrics_tz[1]) || qualtrics_tz[1] == "") {
#   qualtrics_tz <- "Europe/London"
# } else {
#   qualtrics_tz <- qualtrics_tz[1]
# }

real_collection_starts <- tribble(
  ~Q_Language,    ~collection_start_chr,
  "FR-CM-B",      "2026-05-15 09:00:00", # Cameroon
  "AR-EG-B",      "2026-06-25 00:00:00", # Egypt
  "AR-KW-B",      "2026-06-25 00:00:00", # Kuwait
  "EN-GB-HK-B",   "2026-05-18 00:00:00", # Hong Kong
  "EN-GB-ZB-B",   "2026-05-18 00:00:00", # Zimbabwe
  "EN-GB-ZA-B",   "2026-05-18 00:00:00"  # South Africa
) %>%
  mutate(
    collection_start = ymd_hms(collection_start_chr, tz = "UTC") # qualtrics_tz
  ) %>%
  select(Q_Language, collection_start)

# store original Test values once, for audit trail
if (!"Test_original" %in% names(d)) {
  d <- d %>%
    mutate(Test_original = Test)
}

d <- d %>%
  mutate(
    RecordedDate_dt = as_qualtrics_datetime(RecordedDate, tz = qualtrics_tz)
  ) %>%
  left_join(real_collection_starts, by = "Q_Language") %>%
  mutate(
    real_link_test_before_collection = fetch_id == 4 &
      !is.na(collection_start) &
      RecordedDate_dt < collection_start,
    
    Test = if_else(
      real_link_test_before_collection,
      "Test",
      as.character(Test)
    )
  ) %>%
  select(-collection_start)

d %>%
  filter(real_link_test_before_collection) %>%
  select(
    ResponseId,
    LabID,
    Q_Language,
    RecordedDate,
    RecordedDate_dt,
    Test_original,
    Test
  ) %>%
  arrange(Q_Language, RecordedDate_dt)

# apply filter
d_no_test <- d %>%
  filter(is.na(Test))
nrow(d_no_test)

write_csv(d_no_test, "all_withNorway_noTest.csv")


