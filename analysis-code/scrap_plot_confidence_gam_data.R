# Participant-level observed means
participant_conf <- dat_conf %>%
  filter(
    !is.na(confidence),
    !is.na(correctness),
    !is.na(IH_score)
  ) %>%
  group_by(
    participant,
    condition,
    correctness
  ) %>%
  summarise(
    IH_score = first(IH_score),
    mean_confidence = mean(confidence, na.rm = TRUE),
    n_trials = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Phase = recode(
      condition,
      "intuitive" = "Intuitive phase",
      "reflective" = "Reflective phase"
    )
  )


ggplot(
  pred_conf,
  aes(
    x = IH_score,
    y = fit,
    colour = correctness,
    fill = correctness
  )
) +
  
  # Participant-level observed means
  geom_point(
    data = participant_conf,
    aes(
      x = IH_score,
      y = mean_confidence,
      colour = correctness
    ),
    inherit.aes = FALSE,
    alpha = 0.08,
    size = 0.7,
    position = position_jitter(
      width = 0.025,
      height = 0
    )
  ) +
  
  # Model confidence intervals
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # Model-predicted functions
  geom_line(
    linewidth = 1.1
  ) +
  
  facet_wrap(
    ~Phase,
    nrow = 1
  ) +
  
  labs(
    x = "Intellectual humility",
    y = "Confidence",
    colour = "Response accuracy",
    fill = "Response accuracy"
  ) +
  
  ml_theme()
