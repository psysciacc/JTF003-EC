# THis script collate all steps for producing the dataset in the final form

hablar::set_wd_to_script_path()
rm(list=ls())

# load and combina all data (online + offline surve in Norway)
# source("fetch_data_v4.R")

# data cleaning — remove `Test` responses, taking into account exceptions in
# "Test Link and Other Deviations"
source("data_cleaning.R")

# Finally fill in missing Country responses were possible (Norway + unique langauges)
d <- d_no_test
source("fill_in_country.R")

# ---------------------------------------------------------------
# any additional steps here:

# remove people who failed to give consent
excluded_consent <- nrow(d) - nrow(d %>% filter(Consent == "I consent")) # 164
remaining_after_consent <- nrow(d %>% filter(Consent == "I consent")) # 24087
d <- d %>% 
  filter(Consent == "I consent")

# fully anonimise
d$Gender_4_TEXT[which(d$Gender_4_TEXT == "my name is Zerin fatima. i am a students of IUBAT. I want to be a teacher")] <- "my name is ##### ######. i am a students of IUBAT. I want to be a teacher"
d$PROLIFIC_PID <- NULL

# ---------------------------------------------------------------
# save final dataset
write_csv(d, "all_withNorway_noTest_countryFill_raw.csv")

# Copy it to analysis folder
file.copy("all_withNorway_noTest_countryFill_raw.csv", 
          "../analysis-code/all_withNorway_noTest_countryFill_raw.csv",
          overwrite = TRUE)

