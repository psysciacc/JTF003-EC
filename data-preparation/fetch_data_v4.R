# Data extraction script v001
# Matteo Lisi 2026, based on tracker code by Erin Buchanan

hablar::set_wd_to_script_path()
rm(list=ls())

library(tidyverse)
library(purrr)
library(qualtRics)

# helpers
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

clean_qtext <- function(x) {
  x %>%
    str_replace_all("<[^>]+>", " ") %>%
    str_squish()
}

# load API keys 
# NOTE: API key not included in github repository, so this will throw an error
source("../tracker_code/api_key.R")

# copied from ../study_data.R
survey_id <- "SV_9mNnD0W9qxsEgrY"
survey_id2 <- "SV_2013PbPEco6rGhE"
survey_id3 <- "SV_d4N9JjGJyLzedim"
survey_id4 <- "SV_9tceQThx7Pm9JVc"

# ---------------------------------------------------------------------------
# make lookup tables of QID

get_qid_lookup <- function(survey_id, fetch_id) {
  
  desc <- fetch_description(
    surveyID = survey_id,
    elements = "questions"
  )
  
  qs <- desc$questions %||% desc$result$questions
  
  # Works for the usual list-of-questions structure returned by Qualtrics
  lookup <- imap_dfr(qs, function(q, nm) {
    tibble(
      fetch_id = fetch_id,
      survey_id = survey_id,
      qid = q$QuestionID %||% q$questionId %||% nm,
      question_text = q$QuestionText %||% q$questionText %||% NA_character_
    )
  }) %>%
    mutate(
      question_text = clean_qtext(question_text)
    ) %>%
    filter(!is.na(qid), !is.na(question_text))
  
  lookup
}

survey_lookup <- tibble(
  fetch_id = c(1, 2, 3, 4),
  survey_id = c(survey_id, survey_id2, survey_id3, survey_id4)
)

qid_lookup <- survey_lookup %>%
  mutate(lookup = map2(survey_id, fetch_id, get_qid_lookup)) %>%
  select(lookup) %>%
  unnest(lookup)

# # sanity check
# qid_lookup %>%
#   filter(qid %in% c("QID1124", "QID1127", "QID1119", "QID1122", "QID1120", "QID1118"))

# -----------------------------
# Here we define a resolver/replacer function to replaced piped text across all columns
resolve_piped_question_text <- function(data, qid_lookup, id_col = "fetch_id") {
  
  lookup_vec <- setNames(
    qid_lookup$question_text,
    paste(qid_lookup[[id_col]], qid_lookup$qid, sep = "__")
  )
  
  replace_one_value <- function(value, fetch_id) {
    
    if (is.na(value)) return(NA_character_)
    
    if (!str_detect(value, "\\$\\{q://QID\\d+/QuestionText\\}")) {
      return(value)
    }
    
    qids <- str_extract_all(value, "QID\\d+")[[1]]
    out <- value
    
    for (qid in qids) {
      key <- paste(fetch_id, qid, sep = "__")
      replacement <- lookup_vec[[key]]
      
      if (!is.null(replacement) && !is.na(replacement)) {
        token <- paste0("${q://", qid, "/QuestionText}")
        out <- gsub(token, replacement, out, fixed = TRUE)
      }
    }
    
    out
  }
  
  replace_vector <- function(x, fetch_id_vec) {
    map2_chr(x, fetch_id_vec, replace_one_value)
  }
  
  char_cols <- names(data)[map_lgl(data, is.character)]
  
  data %>%
    mutate(across(
      all_of(char_cols),
      ~ replace_vector(.x, .data[[id_col]])
    ))
}


# --------------------------------------------------------------------------
# other misc/helpers

# Identify the Qualtrics response-ID column
response_id_col <- function(x) {
  if ("ResponseId" %in% names(x)) return("ResponseId")
  if ("ResponseID" %in% names(x)) return("ResponseID")
  
  stop("Could not find ResponseId or ResponseID.")
}


# Replace selected labelled columns with their Qualtrics recode values
replace_with_numeric_codes <- function(
    labelled,
    coded,
    pattern = "^crt[1-6]_[ir]$"
) {
  
  labelled_id <- response_id_col(labelled)
  coded_id   <- response_id_col(coded)
  
  code_cols <- intersect(
    names(labelled)[str_detect(names(labelled), pattern)],
    names(coded)
  )
  
  if (anyDuplicated(coded[[coded_id]])) {
    stop("Response IDs are not unique in the coded export.")
  }
  
  row_match <- match(
    labelled[[labelled_id]],
    coded[[coded_id]]
  )
  
  if (anyNA(row_match)) {
    stop("Some responses in the labelled export were not found in the coded export.")
  }
  
  for (nm in code_cols) {
    
    # Preserve the variable/question label attribute
    variable_label <- attr(labelled[[nm]], "label")
    
    labelled[[nm]] <- readr::parse_double(
      as.character(coded[[nm]][row_match]),
      na = c("", "NA")
    )
    
    attr(labelled[[nm]], "label") <- variable_label
  }
  
  labelled
}


# Fetch a mostly-labelled survey, with selected variables kept as codes
fetch_survey_mixed <- function(
    survey_id,
    code_pattern = "^crt[1-6]_[ir]$",
    ...
) {
  
  labelled <- fetch_survey(
    surveyID = survey_id,
    ...,
    convert = FALSE,
    label = TRUE
  )
  
  coded <- fetch_survey(
    surveyID = survey_id,
    ...,
    convert = FALSE,
    label = FALSE
  )
  
  replace_with_numeric_codes(
    labelled = labelled,
    coded = coded,
    pattern = code_pattern
  )
}
# ---------------------------------------------------------------------------
# also from study_data
the_study <- fetch_survey_mixed(survey_id) %>%
  mutate(LabID = if_else(Q_Language == "EN-NL-381", "381", LabID)) %>%  #need to change in data
  mutate(Test = if_else(Q_Language == "EN-NG-187C", NA_character_, Test)) %>% #need to change in data except first session
  mutate(LabID = as.character(LabID),
         id = as.character(id),
         fetch_id = 1)

the_study_2 <- fetch_survey_mixed(survey_id2) %>%
  mutate(LabID = as.character(LabID),
         fetch_id = 2) 

the_study_3 <- fetch_survey_mixed(survey_id3) %>%
  mutate(LabID = as.character(LabID))%>%
  mutate(id=NA,
         fetch_id = 3)

the_study_4 <- fetch_survey_mixed(survey_id4) %>%
  mutate(LabID = as.character(LabID))%>%
  mutate(id=NA,
         fetch_id = 4)

# ------------------------------------------------------------------------------
# ADD NORWAY DATA (collected offline)
dN_labelled <- read_survey(
  "norway_data/PSA_Errorcorrection_April+14,+2026_16.02_labels.csv",
  strip_html = TRUE,
  import_id = FALSE,
  add_column_map = TRUE,
  add_var_labels = TRUE
)

dN_coded <- read_survey(
  "norway_data/PSA_Errorcorrection_April+14,+2026_16.03_values.csv",
  strip_html = TRUE,
  import_id = FALSE,
  add_column_map = TRUE,
  add_var_labels = TRUE
)

dN <- replace_with_numeric_codes(dN_labelled, dN_coded) %>%
  mutate(
    LabID = as.character(LabID),
    id = NA_character_,
    fetch_id = 5,
    Country = "NOR"
  )

# Norway should uses the same survey/QID structure as the online surveys.
# Re-use fetch_id == 1 lookup, but relabel it as fetch_id == 5.
qid_lookup_norway <- qid_lookup %>%
  filter(fetch_id == 1) %>%
  mutate(
    fetch_id = 5,
    survey_id = "Norway_local_export"
  )
qid_lookup_all <- bind_rows(qid_lookup, qid_lookup_norway)

# ------------------------------------------------------------------------------
# Now compile and save ALL data 

# alternative
dN$Age <- as.character(dN$Age)
dN$Defocus_main <- as.character( dN$Defocus_main)
dN$Defocus_crt_r <- as.character(dN$Defocus_crt_r)
dN$Defocus_crt_i <- as.character(dN$Defocus_crt_i)

all_data <- bind_rows(the_study, the_study_2, the_study_3, the_study_4, dN) %>% 
  mutate(Country = if_else(Country == "CS", "CZ",
                           if_else(Country == "ZH", "CN",
                                   if_else(LabID == "266", "UK",
                                           if_else(LabID == "381", "NL",
                                                   if_else(Country == "ET,ET", "ET", 
                                                           if_else(Country == "NG,NG", "NG",
                                                                   if_else(Country == "MY,MY", "MY",
                                                                           Country)))))))) %>%
  mutate(Country = if_else(LabID == "381", "NL", Country)) %>% 
  mutate(Country = if_else(LabID == "397", "IT", Country)) %>% 
  mutate(Country = if_else(LabID == "2796", "GR", Country)) %>% 
  mutate(LabID = case_when(
    LabID == "1257,1257" ~ "1257",
    LabID == "187,187" ~ "187",
    LabID == "64,64" ~ "64", 
    TRUE ~ LabID
  )) %>% 
  mutate(Sample = str_replace(trimws(Sample), ".*S.*", "S")) %>%
  mutate(Sample = recode(Sample, "C" = "Community", "S" = "Student")) %>%
  select(-matches("_DO_")) %>%
  resolve_piped_question_text(qid_lookup_all)

str(all_data)

write_csv(all_data, "all_withNorway_raw.csv")

# these returns TRUE
# all(all_data$fetch_id[which(!is.na(all_data$ResponseID))]==2)
# all(the_study_2$ResponseId == the_study_2$ResponseID)
