# participant_conf <- dat_conf %>%
#   filter(
#     !is.na(confidence),
#     !is.na(correctness),
#     !is.na(IH_score)
#   ) %>%
#   group_by(
#     participant,
#     condition,
#     correctness
#   ) %>%
#   summarise(
#     IH_score = first(IH_score),
#     mean_confidence = mean(confidence, na.rm = TRUE),
#     n_trials = n(),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Phase = recode(
#       condition,
#       "intuitive"  = "Intuitive phase",
#       "reflective" = "Reflective phase"
#     )
#   )
# 

binned_conf_trial <- dat_conf %>%
  mutate(
    Phase = recode(
      condition,
      "intuitive"  = "Intuitive phase",
      "reflective" = "Reflective phase"
    )
  ) %>%
  group_by(Phase, correctness) %>%
  mutate(
    IH_bin = ntile(IH_score, 20)
  ) %>%
  group_by(Phase, correctness, IH_bin) %>%
  summarise(
    IH_score = mean(IH_score, na.rm = TRUE),
    mean_confidence = mean(confidence, na.rm = TRUE),
    n_trials = n(),
    se_confidence =
      sd(confidence, na.rm = TRUE) / sqrt(n_trials),
    .groups = "drop"
  )

adjust_for_problem <- function(model, data) {
  
  # Prediction using the participant's actual problem
  pred_actual <- predict(
    model$gam,
    newdata = data,
    type = "response"
  )
  
  # For every observation, predict it at every problem
  problem_levels <- levels(data$problem)
  
  pred_average_problem <- sapply(
    problem_levels,
    function(p) {
      
      nd <- data
      nd$problem <- factor(
        p,
        levels = problem_levels
      )
      
      predict(
        model$gam,
        newdata = nd,
        type = "response"
      )
    }
  ) %>%
    rowMeans()
  
  # Difference caused specifically by the actual problem
  problem_effect <- pred_actual - pred_average_problem
  
  data %>%
    mutate(
      confidence_adjusted =
        confidence - problem_effect
    )
}

dat_int_adj <- adjust_for_problem(
  m_conf_int,
  dat_int
)

dat_ref_adj <- adjust_for_problem(
  m_conf_ref,
  dat_ref
)

dat_conf_adj <- bind_rows(
  dat_int_adj %>%
    mutate(Phase = "Intuitive phase"),
  
  dat_ref_adj %>%
    mutate(Phase = "Reflective phase")
)

participant_conf_adj <- dat_conf_adj %>%
  group_by(
    participant,
    Phase,
    correctness
  ) %>%
  summarise(
    IH_score = first(IH_score),
    
    mean_confidence =
      mean(confidence_adjusted, na.rm = TRUE),
    
    .groups = "drop"
  )

binned_conf_adj <- participant_conf_adj %>%
  group_by(Phase, correctness) %>%
  mutate(
    IH_bin = ntile(IH_score, 30)
  ) %>%
  group_by(Phase, correctness, IH_bin) %>%
  summarise(
    IH_score = mean(IH_score, na.rm = TRUE),
    
    n_participants = n(),
    
    se_confidence =
      sd(mean_confidence, na.rm = TRUE) /
      sqrt(n_participants),
    
    mean_confidence =
      mean(mean_confidence, na.rm = TRUE),
    
    .groups = "drop"
  )

mean_accuracy <- dat_conf %>%
  group_by(condition, participant) %>%
  summarise(
    accuracy = mean(accuracy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(condition) %>%
  summarise(
    accuracy_pct = 100 * mean(accuracy),
    .groups = "drop"
  ) %>%
  mutate(
    Phase = recode(
      condition,
      "intuitive" = "Intuitive phase",
      "reflective" = "Reflective phase"
    )
  )

# binned_conf_trial <- dat_conf %>%
#   mutate(
#     Phase = recode(
#       condition,
#       "intuitive"  = "Intuitive phase",
#       "reflective" = "Reflective phase"
#     )
#   ) %>%
#   group_by(Phase, correctness) %>%
#   mutate(
#     IH_bin = ntile(IH_score, 10)
#   ) %>%
#   group_by(Phase, correctness, IH_bin) %>%
#   summarise(
#     IH_score = mean(IH_score, na.rm = TRUE),
#     mean_confidence = mean(confidence, na.rm = TRUE),
#     n_trials = n(),
#     se_confidence =
#       sd(confidence, na.rm = TRUE) / sqrt(n_trials),
#     .groups = "drop"
#   )

ggplot(
  pred_conf,
  aes(
    x = IH_score,
    y = fit,
    colour = correctness,
    fill = correctness
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  
  # SE around binned observed means
  geom_errorbar(
    data = binned_conf_adj,
    aes(
      x = IH_score,
      ymin = mean_confidence - se_confidence,
      ymax = mean_confidence + se_confidence,
      colour = correctness
    ),
    inherit.aes = FALSE,
    width = 0,
    linewidth = 0.5
  ) +
  
  geom_point(
    data = binned_conf_adj,
    aes(
      x = IH_score,
      y = mean_confidence,
      colour = correctness
    ),
    pch=21,
    inherit.aes = FALSE,
    size = 2
  ) +
  
  # geom_line(
  #   data = binned_conf_adj,
  #   aes(
  #     x = IH_score,
  #     y = mean_confidence,
  #     colour = correctness
  #   ),
  #   inherit.aes = FALSE,
  #   linewidth = 0.7,
  #   alpha = 0.7
  # ) +
  
  facet_wrap(~Phase, nrow = 1) +
  
  labs(
    x = "Intellectual humility",
    y = "Confidence",
    colour = "Response accuracy",
    fill = "Response accuracy"
  ) +
  
  geom_hline(
    data = mean_accuracy,
    aes(yintercept = accuracy_pct),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.8
  ) +
  
  geom_text(
    data = mean_accuracy,
    aes(
      x = Inf,
      y = accuracy_pct,
      label = sprintf("Accuracy: %.1f%%", accuracy_pct)
    ),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = -0.5,
    size = 3.5,
    colour = "grey30"
  ) +
  
  ml_theme()

# # Participant-level observed means
# participant_conf <- dat_conf %>%
#   filter(
#     !is.na(confidence),
#     !is.na(correctness),
#     !is.na(IH_score)
#   ) %>%
#   group_by(
#     participant,
#     condition,
#     correctness
#   ) %>%
#   summarise(
#     IH_score = first(IH_score),
#     mean_confidence = mean(confidence, na.rm = TRUE),
#     n_trials = n(),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Phase = recode(
#       condition,
#       "intuitive" = "Intuitive phase",
#       "reflective" = "Reflective phase"
#     )
#   )
# 
# 
# 
# 
# ggplot(
#   pred_conf,
#   aes(
#     x = IH_score,
#     y = fit,
#     colour = correctness,
#     fill = correctness
#   )
# ) +
#   
#   # Participant-level observed means
#   geom_point(
#     data = participant_conf,
#     aes(
#       x = IH_score,
#       y = mean_confidence,
#       colour = correctness
#     ),
#     inherit.aes = FALSE,
#     alpha = 0.08,
#     size = 0.7,
#     position = position_jitter(
#       width = 0.025,
#       height = 0
#     )
#   ) +
#   
#   # Model confidence intervals
#   geom_ribbon(
#     aes(
#       ymin = lower,
#       ymax = upper
#     ),
#     alpha = 0.15,
#     colour = NA
#   ) +
#   
#   # Model-predicted functions
#   geom_line(
#     linewidth = 1.1
#   ) +
#   
#   facet_wrap(
#     ~Phase,
#     nrow = 1
#   ) +
#   
#   labs(
#     x = "Intellectual humility",
#     y = "Confidence",
#     colour = "Response accuracy",
#     fill = "Response accuracy"
#   ) +
#   
#   ml_theme()
