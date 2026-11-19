

tab_N_condition <- d %>%
  mutate(
    condition = case_when(
      feedback == "F"  & justification == "NJ" ~ "feedback",
      feedback == "NF" & justification == "J"  ~ "justification",
      feedback == "F"  & justification == "J"  ~ "feedback + justification",
      .default = "baseline"
    )
  ) %>%
  count(condition, name = "N") %>%
  knitr::kable(
    booktabs = TRUE,
    caption = "Sample size by experimental condition."
  )


d %>%
  mutate(
    condition = case_when(
      feedback == "F"  & justification == "NJ" ~ "feedback",
      feedback == "NF" & justification == "J"  ~ "justification",
      feedback == "F"  & justification == "J"  ~ "feedback + justification",
      .default = "baseline"
    )
  ) %>%
  count(condition, culture, name = "N") %>%
  pivot_wider(
    names_from = culture,
    values_from = N,
    values_fill = 0
  ) %>%
  knitr::kable(
    booktabs = TRUE,
    caption = "Sample size by experimental condition and cultural group."
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("scale_down", "striped")
  )

###explicit condition

d_desc <- d_all_H1 %>%
  mutate(
    Condition = case_when(
      feedback == "NF" & justification == "NJ" ~ "Baseline",
      feedback == "F"  & justification == "NJ" ~ "Feedback",
      feedback == "NF" & justification == "J"  ~ "Justification",
      feedback == "F"  & justification == "J"  ~ "Feedback + justification"
    )
  )



# timeouts
timeouts <- d_desc %>%
  group_by(Condition) %>%
  summarise(
    N = n(),
    n_timeouts = sum(is.na(response_intuitive)),
    pct_timeouts = 100 * n_timeouts / N,
    .groups = "drop"
  ) %>%
  bind_rows(
    d_desc %>%
      summarise(
        Condition = "Overall",
        N = n(),
        n_timeouts = sum(is.na(response_intuitive)),
        pct_timeouts = 100 * n_timeouts / N
      )
  )

knitr::kable(timeouts, digits = 2)

# correct intuitive changed to error
reversals <- d_desc %>%
  filter(accuracy_intuitive == 1,
         !is.na(accuracy_reflective)) %>%
  group_by(Condition) %>%
  summarise(
    N = n(),
    n_reversals = sum(accuracy_reflective == 0),
    pct_reversals = 100 * n_reversals / N,
    .groups = "drop"
  ) %>%
  bind_rows(
    d_desc %>%
      filter(accuracy_intuitive == 1,
             !is.na(accuracy_reflective)) %>%
      summarise(
        Condition = "Overall",
        N = n(),
        n_reversals = sum(accuracy_reflective == 0),
        pct_reversals = 100 * n_reversals / N
      )
  )

knitr::kable(reversals, digits = 2)

## correct responses arising through correction

#d_desc <- d_all_H1 %>%
d_desc <- d_all %>%
  filter(#accuracy_intuitive == 0,        # keep only intuitive errors
         ext_resources_use=="No")  %>%
  mutate(
    Condition = case_when(
      feedback == "NF" & justification == "NJ" ~ "Control",
      feedback == "F"  & justification == "NJ" ~ "Feedback",
      feedback == "NF" & justification == "J"  ~ "Justification",
      feedback == "F"  & justification == "J"  ~ "Feedback + justification"
    )
  )

summarise_patterns <- function(x) {
  x %>%
    summarise(
      N = n(),
      `Correct → Correct` = sum(accuracy_intuitive == 1 & accuracy_reflective == 1, na.rm = TRUE),
      `Correct → Incorrect` = sum(accuracy_intuitive == 1 & accuracy_reflective == 0, na.rm = TRUE),
      `Incorrect → Correct` = sum(accuracy_intuitive == 0 & accuracy_reflective == 1, na.rm = TRUE),
      `Incorrect → Incorrect` = sum(accuracy_intuitive == 0 & accuracy_reflective == 0, na.rm = TRUE),
      Missing = sum(is.na(accuracy_intuitive) | is.na(accuracy_reflective)),
      .groups = "drop"
    ) %>%
    mutate(
      `% correct via correction` = round(
        100 * `Incorrect → Correct` /
          (`Correct → Correct` + `Incorrect → Correct`), 1
      )
    )
}

h1_descriptives <- bind_rows(
  d_desc %>%
    group_by(Condition) %>%
    summarise_patterns(),
  
  d_desc %>%
    summarise_patterns() %>%
    mutate(Condition = "Overall")
)

knitr::kable(h1_descriptives, digits = 1)



# Participant-level descriptives
participant_desc <- d_all_H1 %>%
  mutate(
    Condition = case_when(
      feedback == "NF" & justification == "NJ" ~ "Control",
      feedback == "F"  & justification == "NJ" ~ "Feedback",
      feedback == "NF" & justification == "J"  ~ "Justification",
      feedback == "F"  & justification == "J"  ~ "Feedback + justification"
    )
  ) %>%
  group_by(Condition, participant) %>%
  summarise(
    CC = sum(accuracy_intuitive == 1 & accuracy_reflective == 1, na.rm = TRUE),
    CI = sum(accuracy_intuitive == 1 & accuracy_reflective == 0, na.rm = TRUE),
    IC = sum(accuracy_intuitive == 0 & accuracy_reflective == 1, na.rm = TRUE),
    II = sum(accuracy_intuitive == 0 & accuracy_reflective == 0, na.rm = TRUE),
    Missing = sum(is.na(accuracy_intuitive) | is.na(accuracy_reflective)),
    .groups = "drop"
  ) %>%
  mutate(
    pct_corrected = if_else(
      CC + IC > 0,
      100 * IC / (CC + IC),
      NA_real_
    )
  )

# Format as Mean (SD)
fmt <- function(x) {
  if (all(is.na(x))) return("NA")
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

# Summary by condition
h1_descriptives <- participant_desc %>%
  group_by(Condition) %>%
  summarise(
    N_participants = n(),
    `Correct → Correct` = fmt(CC),
    `Correct → Incorrect` = fmt(CI),
    `Incorrect → Correct` = fmt(IC),
    `Incorrect → Incorrect` = fmt(II),
    Missing = fmt(Missing),
    `% correct via correction` = fmt(pct_corrected),
    N_pct = sum(!is.na(pct_corrected)),
    .groups = "drop"
  ) %>%
  arrange(factor(Condition, levels = c(
    "Control", "Feedback", "Justification",
    "Feedback + justification"
  )))

knitr::kable(
  h1_descriptives,
  caption = "H1 response patterns: participant-level mean (SD)"
)
