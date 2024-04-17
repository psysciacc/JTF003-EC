#!/usr/bin/Rscript
# expect command line arguments in order

cmd_args=commandArgs(TRUE)
N = as.numeric(cmd_args[1])
N_c = as.numeric(cmd_args[2])
N_rep = as.numeric(cmd_args[3])
sim_label = cmd_args[4]


# library
library(tidyverse)

# parameter settings ----

# these are plausible parameter values based on Sirota 2023 and similar 
# article in the literature

# p(error @ stage 1)
err_mu <- 1.22
err_sd <- 0.92

# p(correct @ stage 2 | error @ stage 1) 
corr_mu <- -2.53  
corr_sd <- 1.12

# set effect sizes ----

# determine the country-level random effect variance
rfx_sd_country  <- 0.2 

# Set effect values in standardized as Cohen's d and rescale 
# to an approximately similar value for the logistic distribution; 
# specifically we multiply the d values by the ratio of the standard deviation
# of the logistic distribution to that of the standard normal distribution, 
# which is (pi/sqrt(3)).
d_fdbk <- 0.2 * (pi/sqrt(3)) # feedback minus baseline
d_just <- 0.2 * (pi/sqrt(3)) # justification minus baseline
d_both <- (d_fdbk + d_just)  # (feedback + justification) minus baseline 
# Note: we expects the effects of feedback to be additive

# vector of fixed-effects coefficients
beta_lo_p <- c(corr_mu, NA, NA, NA)
beta_lo_p[2] <- d_fdbk
beta_lo_p[3] <- d_just
beta_lo_p[4] <- d_both

# power simulation parameters ----
n_trial <- 6

# wrapper function to simulate data ----

# first custom function to sample countries from clusters with some constraints
# with how little and how many countries there will be for each cluster
# (based on the results of data collection in Bago et al 2022, https://doi.org/10.1038/s41562-022-01319-5)
generate_labels <- function(N_c) {
  
  # min/max occurrences for each cluster
  min_occurrences <- round(0.15 * N_c)
  max_occurrences <- round(0.65 * N_c)
  
  # sample with constraints
  N_southern <- sample(min_occurrences:max_occurrences, 1)
  N_eastern <- sample(min_occurrences:max_occurrences, 1)
  while((N_southern + N_eastern) > (N_c - min_occurrences) || (N_southern + N_eastern) < (N_c - max_occurrences)) {
    N_southern <- sample(min_occurrences:max_occurrences, 1)
    N_eastern <- sample(min_occurrences:max_occurrences, 1)
  }
  N_western <- N_c - (N_southern + N_eastern)
  labels <- c(rep("Southern", N_southern), rep("Eastern", N_eastern), rep("Western", N_western))
  labels <- sample(labels)
  
  return(labels)
}

# data simulation function
gen_data <- function(N, N_c){
  
  N_tot <- N * N_c * 4 # 4 is the number of between subjects conditions
  
  ## simulate dataset
  
  # simulate number of errors in phase 1
  error_lo <- rnorm(N_tot, mean=err_mu, sd=err_sd)
  
  # simulate random effects and between-individual variability 
  # (to introduce heterogeneity and make simulations more realistic)
  rfx_id <- rnorm(N_tot, mean=0, sd=corr_sd)
  
  # sample countries
  rfx_country <- rnorm(N_c, mean=0, sd=rfx_sd_country)
  cluster_country <- generate_labels(N_c)   
  # create dataset frame
  dat <- expand_grid(country = 1:N_c,
                     condition = 1:4,
                     id = 1:N,
                     trial=1:6) %>%
    mutate(id = str_c(id,country,condition,sep="_"))
  dat$cluster <- cluster_country[dat$country]
  
  # simulate responses
  dat$corrected <- NA
  dat$error <- NA
  dat$id <- as.integer(factor(dat$id, labels=1:length(unique(dat$id))))
  dat$lo <- NA
  dat$error_lo <- error_lo[dat$id]
  dat$rfx_country <- rfx_country[dat$country]
  dat$rfx_id <- rfx_id[dat$id]
  
  # set log-odds of p(correction) for each datapoint
  dat <- dat %>%
    mutate(lo = beta_lo_p[1]+rfx_country+rfx_id) %>%
    mutate(lo = case_when(
      condition==1 ~ lo,
      condition==2 ~ lo + beta_lo_p[2],
      condition==3 ~ lo + beta_lo_p[3],
      condition==4 ~ lo + beta_lo_p[4]
    ))
  
  # simulate number of errors
  dat$error <- rbinom(nrow(dat), size=1, p=exp(dat$error_lo)/(1+exp(dat$error_lo)))
  
  # simulate number of errors that are corrected in phase 2
  dat$corrected <- ifelse(dat$error==1,
                          rbinom(nrow(dat), size=1, p=exp(dat$lo)/(1+exp(dat$lo))),
                          0)
  
  # total correct responses
  dat$correct <- 1-dat$error + dat$corrected
  
  # keep only trials with errors in phase 1
  dat <- dat %>%
    # filter(error==1) %>%
    mutate(c_F = ifelse(condition==2,1,0),
           c_FR = ifelse(condition==3,1,0)) %>%
    mutate(feedback=ifelse(condition==2 | condition==4, 1, 0),
           justification=ifelse(condition==3 | condition==4, 1, 0),
           condition = case_when(
             condition==1 ~ "baseline",
             condition==2 ~ "feedback",
             condition==3 ~ "justification",
             condition==4 ~ "feedback+justification"
           ))
  
  dat$cluster <- factor(dat$cluster)
  contrasts(dat$cluster) <- contr.treatment(levels(dat$cluster), base=3)
  
  # Set attributes
  attr(dat, "N") <- N
  attr(dat, "N_c") <- N_c
  
  # return simulated data
  return(dat)
}

# wrapper function to fit models and test hypotheses ----

test_hypotheses <- function(dat){
  
  require(lme4)
  # dat <- gen_data(100, 30)
  
  dat_H1 <- dat %>%
    filter(correct==1)
  
  dat <- dat %>%
    filter(error==1)
  
  ## FIT MODELS
  
  # fit overall model
  m_all <- glmer(corrected ~ feedback * justification + (1|country),
                 family=binomial("logit"), data=dat,
                 control = glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=1e5)))
  
  # fit cluster-level models
  m_west <- glmer(corrected ~  feedback * justification + (1|country),
                  family=binomial("logit"), data=dat[dat$cluster=="Western",],
                  control = glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=1e5)))
  
  m_south <- glmer(corrected ~  feedback * justification + (1|country),
                   family=binomial("logit"), data=dat[dat$cluster=="Southern",],
                   control = glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=1e5)))
  
  m_east <- glmer(corrected ~  feedback * justification + (1|country),
                  family=binomial("logit"), data=dat[dat$cluster=="Eastern",],
                  control = glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=1e5)))
  
  # now fit models for H1 specifically: note that here the dependent variable is different
  # as the hypothesis is about the proportion of correct response (out of all correct responses)
  # that are corrected in phase 2 ('data' arguments is set to 'dat_H1' instead of just 'dat')
  # 
  mH1_all <- glmer(corrected ~ feedback * justification + (1|country),
                 family=binomial("logit"), data=dat_H1,
                 control = glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=1e5)))

  # fit cluster-level models
  mH1_west <- glmer(corrected ~ feedback * justification + (1|country),
                  family=binomial("logit"), data=dat_H1[dat_H1$cluster=="Western",],
                  control = glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=1e5)))

  mH1_south <- glmer(corrected ~ feedback * justification + (1|country),
                   family=binomial("logit"), data=dat_H1[dat_H1$cluster=="Southern",],
                   control = glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=1e5)))

  mH1_east <- glmer(corrected ~ feedback * justification + (1|country),
                  family=binomial("logit"), data=dat_H1[dat_H1$cluster=="Eastern",],
                  control = glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=1e5)))
  
  ## TEST HYPOTHESES
  
  alpha <- 0.05
  # note: the tests on separate clusters are independent one another 
  # (different countries, participants, labs, etc.), so the expected false positive
  # rate of the composite test, run on the 3 cluster separately and using the 
  # conventional alpha=0.05 for each test, would be 0.05^3 = 0.000125 (the 
  # hypothesis supported if and only if all 3 independent tests are significant)
  # We don't also apply correction for multiple comparisons otherwise it would become 
  # too conservative.
  
  # -----------
  # H1 
  # # p_correction < 0.5 
  
  # One-tailed test for H1a
  # Testing if the intercept is significantly less than 0
  b_H1a <- fixef(mH1_all)["(Intercept)"]
  SE_H1a <- sqrt(diag(vcov(mH1_all))["(Intercept)"])
  z_H1a <- b_H1a / SE_H1a
  p_value_H1a <- pnorm(z_H1a) # left-tailed test
  H1a <- ifelse(p_value_H1a < alpha, 1, 0)
  
  # One-tailed tests for H1b - for 'west', 'south', 'east'
  models <- list(mH1_west, mH1_south, mH1_east)
  names <- c("west", "south", "east")
  H1b_pvals <- numeric(length(models))
  H1b_results <- numeric(length(models))
  
  for (i in 1:length(models)) {
    b <- fixef(models[[i]])["(Intercept)"]
    SE <- sqrt(diag(vcov(models[[i]]))["(Intercept)"])
    z <- b / SE
    p_value <- pnorm(z) # left-tailed test
    H1b_pvals[i] <- p_value
    H1b_results[i] <- ifelse(p_value < alpha, 1, 0)
  }
  
  names(H1b_pvals) <- names
  names(H1b_results) <- names
  
  # complete test
  H1b <- ifelse(all(H1b_results == 1), 1, 0)

  # -----------
  # H2
  
  # feedback increases correction rate

  # One-tailed test for H2a
  # Testing if the coefficient for 'feedback' is significantly greater than 0
  b_H2a <- fixef(m_all)["feedback"]
  SE_H2a <- sqrt(diag(vcov(m_all))["feedback"])
  z_H2a <- b_H2a / SE_H2a
  p_value_H2a <- 1 - pnorm(z_H2a) # right-tailed test
  H2a <- ifelse(p_value_H2a < alpha, 1, 0)
  
  # One-tailed tests for H2b - for 'west', 'south', 'east'
  models <- list(m_west, m_south, m_east)
  names <- c("west", "south", "east")
  H2b_pvals <- numeric(length(models))
  H2b_results <- numeric(length(models))
  
  for (i in 1:length(models)) {
    b <- fixef(models[[i]])["feedback"]
    SE <- sqrt(diag(vcov(models[[i]]))["feedback"])
    z <- b / SE
    p_value <- 1 - pnorm(z) # right-tailed test
    H2b_pvals[i] <- p_value
    H2b_results[i] <- ifelse(p_value < alpha, 1, 0)
  }
  
  names(H2b_pvals) <- names
  names(H2b_results) <- names
  
  # complete test
  H2b <- ifelse(all(H2b_results == 1), 1, 0)
  
  # -----------
  # H3
  
  # justification increases correction rate
  
  # One-tailed test for H3a
  # Testing if the coefficient for 'justification' is significantly greater than 0
  b_H3a <- fixef(m_all)["justification"]
  SE_H3a <- sqrt(diag(vcov(m_all))["justification"])
  z_H3a <- b_H3a / SE_H3a
  p_value_H3a <- 1 - pnorm(z_H3a) # right-tailed test
  H3a <- ifelse(p_value_H3a < alpha, 1, 0)
  
  # One-tailed tests for H3b - for 'west', 'south', 'east'
  models <- list(m_west, m_south, m_east)
  names <- c("west", "south", "east")
  H3b_pvals <- numeric(length(models))
  H3b_results <- numeric(length(models))
  
  for (i in 1:length(models)) {
    b <- fixef(models[[i]])["justification"]
    SE <- sqrt(diag(vcov(models[[i]]))["justification"])
    z <- b / SE
    p_value <- 1 - pnorm(z) # right-tailed test
    H3b_pvals[i] <- p_value
    H3b_results[i] <- ifelse(p_value < alpha, 1, 0)
  }
  
  names(H3b_pvals) <- names
  names(H3b_results) <- names
  
  # complete test
  H3b <- ifelse(all(H3b_results == 1), 1, 0)
  
  
  # save results -----------
  
  # retrieve settings
  N <- attr(dat, "N")
  N_c <- attr(dat, "N_c")
  N_tot <- N * N_c * 4
  
  d_line <- data.frame(N_per_condition = N,
                       N_tot = N_tot,
                       N_countries = N_c,
                       H1a, H1b, 
                       H2a, H2b, 
                       H3a, H3b)
  
  return(d_line)
  
}

# ---------------------------------------------------------------------------
# run simulations

filename <- str_c("PSAsim_wnCllvsm_N",N,"_Nc",N_c,"_Nrep",N_rep,"_",sim_label,".csv")
d_res <- {}
for(it in  1:N_rep ){
  d_line <- test_hypotheses(gen_data(N, N_c))
  d_res <- rbind(d_res, d_line)
  write_csv(d_res,filename)
}
