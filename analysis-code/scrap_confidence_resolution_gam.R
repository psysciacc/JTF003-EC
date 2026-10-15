library(dplyr)
library(tidyr)
library(ggplot2)

# Mean/SD used to convert IH_z back to the original IH scale
IH_par <- dat_conf %>%
  distinct(participant, IH_score, IH_z)

IH_mean <- mean(IH_par$IH_score, na.rm = TRUE)
IH_sd   <- sd(IH_par$IH_score, na.rm = TRUE)


# ------------------------------------------------------------
# Function to obtain predicted confidence curves
# ------------------------------------------------------------

make_conf_predictions <- function(model, data, phase, n = 200) {
  
  # IH sequence
  IH_seq <- seq(
    min(data$IH_z, na.rm = TRUE),
    max(data$IH_z, na.rm = TRUE),
    length.out = n
  )
  
  # Prediction grid:
  # every IH × correctness combination × all 6 problems
  newdat <- expand_grid(
    IH_z = IH_seq,
    correctness = factor(
      c("Incorrect", "Correct"),
      levels = c("Incorrect", "Correct")
    ),
    problem = levels(data$problem)
  ) %>%
    mutate(
      problem = factor(
        problem,
        levels = levels(data$problem)
      )
    )
  
  # Preserve sum contrasts
  contrasts(newdat$problem) <- contr.sum(nlevels(data$problem))
  
  # GAM design matrix
  Xp <- predict(
    model$gam,
    newdata = newdat,
    type = "lpmatrix"
  )
  
  # IDs identifying rows to average together
  newdat$pred_id <- interaction(
    newdat$IH_z,
    newdat$correctness,
    drop = TRUE
  )
  
  ids <- unique(newdat$pred_id)
  
  # Average design matrix over the six problems
  Xbar <- do.call(
    rbind,
    lapply(ids, function(id) {
      colMeans(Xp[newdat$pred_id == id, , drop = FALSE])
    })
  )
  
  # Corresponding IH/correctness values
  pred_dat <- newdat %>%
    distinct(pred_id, IH_z, correctness) %>%
    arrange(match(pred_id, ids))
  
  # Fixed/GAM coefficients and covariance matrix
  beta <- coef(model$gam)
  V <- vcov(model$gam)
  
  pred_dat %>%
    mutate(
      fit = as.numeric(Xbar %*% beta),
      
      se = sqrt(
        rowSums((Xbar %*% V) * Xbar)
      ),
      
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      
      Phase = phase
    )
}




pred_conf_int <- make_conf_predictions(
  m_conf_int,
  dat_int,
  "Intuitive phase"
)

pred_conf_ref <- make_conf_predictions(
  m_conf_ref,
  dat_ref,
  "Reflective phase"
)

pred_conf <- bind_rows(
  pred_conf_int,
  pred_conf_ref
)

IH_par <- dat_conf %>%
  distinct(participant, IH_score, IH_z)

IH_mean <- mean(IH_par$IH_score, na.rm = TRUE)
IH_sd   <- sd(IH_par$IH_score, na.rm = TRUE)

pred_conf <- pred_conf %>%
  mutate(
    IH_score = IH_z * IH_sd + IH_mean
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
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~Phase, nrow = 1) +
  labs(
    x = "Intellectual humility",
    y = "Predicted confidence",
    colour = "Response accuracy",
    fill = "Response accuracy"
  ) +
  ml_theme()


