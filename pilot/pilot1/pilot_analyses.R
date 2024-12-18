rm(list=ls())
hablar::set_wd_to_script_path()

library(tidyverse)
library(kableExtra)

d <- read_csv("PSA+Error+Correction+PILOT.csv")
str(d)


# ------------------------------------------------------------------------------
### data checks

# ------------------------------------------------------------------------------
# dropouts
dropout_rates <- d %>%
  filter(Status == 0) %>%
  summarise(Completed = max(Finished == 1)) %>%
  summarise(
    Total_Participants = n(),
    Dropouts = sum(Completed == 0),
    Dropout_Rate = Dropouts / Total_Participants
  )
print(dropout_rates)

# ------------------------------------------------------------------------------
# exclusions
total_initial <- nrow(d)

# Status == 0
excluded_status <- total_initial - nrow(d %>% filter(Status == 0))
remaining_after_status <- nrow(d %>% filter(Status == 0))


# Finished == 1
excluded_finished <- remaining_after_status - nrow(d %>% filter(Status == 0, 
                                                                Finished == 1))
#excluded_finished <- total_initial - nrow(d %>% filter(Finished == 1))
remaining_after_finished <- nrow(d %>% filter(Status == 0, 
                                              Finished == 1))

# Q_RecaptchaScore > 0.5
excluded_recaptcha <- remaining_after_finished - nrow(d %>% filter(Status == 0, 
                                                                   Finished == 1, 
                                                                   Progress >= 99, 
                                                                   Q_RecaptchaScore >= 0.5))
remaining_after_recaptcha <- nrow(d %>% filter(Status == 0, 
                                               Finished == 1, 
                                               Progress >= 99, 
                                               Q_RecaptchaScore >= 0.5))

# Total number of participants excluded (excluded because they fail at least 1 criterion)
total_excluded <- total_initial - remaining_after_recaptcha

# Create summary table
summary_table <- data.frame(
  Criterion = c( "Status == 0",
                 "Finished == 1", 
                 "Q_RecaptchaScore >= 0.5", 
                 "Total"),
  Excluded = c(excluded_status, 
               excluded_finished, 
               excluded_recaptcha, 
               total_excluded),
  Fraction = c(excluded_status/total_initial, 
               excluded_finished/total_initial, 
               excluded_recaptcha/total_initial, 
               total_excluded/total_initial),
  
  Percentage = c(excluded_status/total_initial * 100, 
                 excluded_finished/total_initial * 100, 
                 excluded_recaptcha/total_initial * 100, 
                 total_excluded/total_initial * 100)
)

# Print the summary table as an HTML table
kable(summary_table, format = "html", booktabs = TRUE, 
      col.names = c("Exclusion Criterion", "Number Excluded", "Fraction Excluded", "Percentage Excluded"), 
      caption="Each line indicate how many participants failed each criterion. Exclusion criteria are applied sequentially in order; the last line indicate how many participants failed at least one criterion in total.") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

# ------------------------------------------------------------------------------
demographic_summary <- d %>%
  summarise(
    `Avearge age (Std.)` = sprintf("%.2f (%.2f)", mean(Age, na.rm = TRUE), sd(Age, na.rm = TRUE)),
    `Age Range` = paste(min(Age, na.rm = TRUE), "to", max(Age, na.rm = TRUE)),
    `N. males (Gender = 1)` = sprintf("%i", sum(Gender == 1, na.rm = TRUE)),
    `N. female (Gender = 2)` = sprintf("%i", sum(Gender == 2, na.rm = TRUE)),
    `N. non binary (Gender = 3)` = sprintf("%i", sum(Gender == 3, na.rm = TRUE))
  ) %>%
  pivot_longer(everything(), names_to = "Demographic", values_to = "Value") 

knitr::kable(demographic_summary, format = "html", booktabs = TRUE, caption = "Demographic Information of Participants", col.names = c("Demographic", "Value"))%>%
  kable_styling(latex_options = c("scale_down","striped"))


# ------------------------------------------------------------------------------
# Reshape data for analyzing errors

d_long <- d %>%
  # Add participant ID
  mutate(participant = row_number()) %>%
  # Select participant ID and columns matching 'crt#_i' or 'crt#_r'
  select(participant, matches('^crt[1-6]_[ir]$')) %>%
  # Pivot longer to reshape data into long format
  pivot_longer(
    cols = matches('^crt[1-6]_[ir]$'),
    names_to = 'item',
    values_to = 'response'
  ) %>%
  # Extract 'problem' number and 'condition' from 'item'
  extract(
    item,
    into = c('crt', 'problem', 'condition'),
    regex = '(crt)([1-6])_([ir])'
  ) %>%
  # Convert 'problem' to integer and recode 'condition'
  mutate(
    problem = as.integer(problem),
    condition = recode(condition, 'i' = 'intuitive', 'r' = 'reflexive')
  ) %>%
  # Select and arrange the final columns
  select(participant, problem, condition, response)

# View the transformed data
print(d_long)

d_long$accuracy <- ifelse(!is.na(d_long$response), 
                       ifelse(d_long$response==1,1,0),
                       d_long$response)

mean_accuracy_table <- d_long %>%
  # Group by 'problem' and 'condition'
  group_by(problem, condition) %>%
  # Calculate mean accuracy for each group
  summarize(mean_accuracy = mean(accuracy, na.rm = TRUE)) %>%
  ungroup() %>%
  # Reshape the data to have conditions as columns
  pivot_wider(
    names_from = condition,
    values_from = mean_accuracy
  )

# View the resulting table
print(mean_accuracy_table)

knitr::kable(mean_accuracy_table, format = "html", booktabs = TRUE, caption = "Mean accuracy by CRT problem and phase", digits=2)%>%
  kable_styling(latex_options = c("scale_down","striped"))
