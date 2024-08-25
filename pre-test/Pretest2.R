# ******************************************************************************************
# This file contains the code for the analyses reported in the manuscript:
#
# Sirota, M, et al (2024) Mapping and increasing error correction behaviour in a culturally diverse sample. 
# Manuscript submitted for publication.
# 
# last update: 03/May/2024
#
# This script file is licensed under a CC-BY 4.0 license. 
# see http://creativecommons.org/licenses/by/4.0/
# 
# written by Miroslav Sirota (email: msirota@essex.ac.uk)
# Please email me if you see any errors or have any questions.
# ******************************************************************************************


# ------------------------------------------------------------------------------
# Install and load required libraries
# ------------------------------------------------------------------------------

# Define the list of packages
packages <- c("psych", "dplyr", "MASS", "effsize")

# Check if each package is already installed and if not, install it
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  } else {
    library(package, character.only = TRUE)
  }
}

# ------------------------------------------------------------------------------
# Clear a work space and load data 
# ------------------------------------------------------------------------------

# Clear work space
rm(list = ls())

# Set working directory
setwd(dirname(file.path(rstudioapi::getActiveDocumentContext()$path)))

# Load data
Pre2 <- read.csv("Pre2.csv",header=T,fileEncoding = 'UTF-8-BOM')

# Check the structure of the data
str(Pre2)



#------------------------------------------------------------------------------
# Codebook for "Pretest 2 (Pretest2.csv)
# Form: Variable name / Variable label / Variable values 
#------------------------------------------------------------------------------
#

#Id/ Identification Number/ NA	
#Age/ "What is your age?"/ NA
#Sex/ "What is your sex? If you are considering how to answer, use the sex recorded on your birth certificate or Gender Recognition Certificate."/ 0 = "male", 1 = "female"
#Duration/ The time taken to complete the whole questionnaire (in seconds) / NA
#Pp1_t TO Pp2_t/The time taken to read the practice problem (problems 1 to 2; in seconds)/ NA  
#Pp1_a TO Pp2_a/Answers to the practice problem (problems 1 to 2)/ 1 = "Option A", 2 = "Option B", 3 = "Option C", 4 = "Option D"
#Pp1_ch TO Pp2_ch/Comprehension checks of the practice problem (problems 1 to 2)/ 1 = "Option A", 2 = "Option B", 3 = "Option C", 4 = "Option D"
#CRT1_t TO CRT6_t/The time taken to read the cognitive reflection problem (problems 1 to 6; in seconds)/ NA  
#CRT1_a TO CRT6_a/Answers to the cognitive reflection problem (problems 1 to 6)/ 1 = "Option A", 2 = "Option B", 3 = "Option C", 4 = "Option D"
#CRT1_ch TO CRT6_ch/Comprehension checks of the cognitive reflection problem (problems 1 to 6)/ 1 = "Option A", 2 = "Option B", 3 = "Option C", 4 = "Option D"

# ------------------------------------------------------------------------------
# Sample description 
# ------------------------------------------------------------------------------

##Description

#Age
describe(Pre2$Age)

#Sex
table(Pre2$Sex)
prop.table(table(Pre2$Sex))*100

# ------------------------------------------------------------------------------
# Analysis 
# ------------------------------------------------------------------------------

# Reading times

#Extract relevant variables
crt_columns <- grep("^CRT\\d+_t$", names(Pre2), value = TRUE)
a_columns <- paste0(sub("_t$", "_a", crt_columns))
ch_columns <- paste0(sub("_t$", "_ch", crt_columns))

#Define functions

# Define a function to apply describe with conditions
describe_with_conditions <- function(column, a, ch) {
  describe(column[a == 1 & ch == 1])
}

# Define a function to apply describe with conditions after logging (Bago & De Neys 2019 method)
describe_with_conditions_logged <- function(column, a, ch) {
  logged_data <- log(column) #logged
  desc <- describe(logged_data[a == 1 & ch == 1]) #described
  desc$EXPmean <- exp(desc$mean) #exp mean
  desc$EXPsd <- ((exp(desc$sd) - 1) * exp(desc$mean)) #exp sd
  return (desc)
}

#Reading times all responses (raw)
results1 <- lapply(Pre2[crt_columns], describe)
results1

# Reading times for those who answered with comprehension (raw)
results2 <- mapply(describe_with_conditions, Pre2[crt_columns], Pre2[a_columns], Pre2[ch_columns])
results2

# Reading times for those who answered correctly (logged + exponentiated)
results3 <- mapply(describe_with_conditions_logged, Pre2[crt_columns], Pre2[a_columns], Pre2[ch_columns])
results3


# wrapper function to fit gamma distribution
fit_gamma_and_plot <- function(x, mean_time=NULL, observed_value=NULL, cohen_d_value=NULL, title="", x_limit=50, y_limit=NULL) {
  # Fit the gamma distribution to the data
  fit <- MASS::fitdistr(x, "gamma")
  
  shape <- fit$estimate[1]
  rate <- fit$estimate[2]
  
  # If mean_time is provided, calculate the corresponding percentile
  quantile_value <- NULL
  if (!is.null(mean_time)) {
    quantile_value <- mean_time
  }
  
  # Determine y-axis limit if not provided
  if (is.null(y_limit)) {
    y_limit <- max(hist(x, breaks = 20, plot = FALSE)$density, na.rm = TRUE) * 1.1
  }
  
  # Plot histogram of the data with x-axis and y-axis limits
  hist(x, breaks = 20, freq = FALSE, col="darkgrey", border = "white", xlab = " ", main = title, xlim=c(0, x_limit), ylim=c(0, y_limit))
  curve(dgamma(x, shape = shape, rate = rate), add = TRUE, col = "black", lwd=2)
  
  if (!is.null(quantile_value)) {
    abline(v = quantile_value, col = "blue", lty = 1, lwd=2)
  }
  
  if (!is.null(observed_value)) {
    abline(v = observed_value, col = "red", lty = 2, lwd=2)
  }
  
  # If both mean_time and observed_value are provided, draw a line with arrows and label it with Cohen's d value
  if (!is.null(mean_time) && !is.null(observed_value) && !is.null(cohen_d_value)) {
    mid_point <- (mean_time + observed_value) / 2
    arrows(x0 = mean_time, y0 = y_limit * 0.9, x1 = observed_value, y1 = y_limit * 0.9, code = 3, angle = 90, length = 0.1, col = "purple", lwd = 2)
    text(mid_point, y_limit * 0.95, labels = paste("d =", round(cohen_d_value, 2)), col = "black")
  }
  
  return(list(shape = shape, rate = rate, mean_time = quantile_value, percentile = if (!is.null(observed_value)) pgamma(observed_value, shape, rate) else NULL))
}

# Calculate effect sizes when compared with solution times in the control condition of Manipulation Pretest
data <- data.frame(
  measure = c("crt1", "crt2", "crt3", "crt4", "crt5", "crt6"),
  mean1 = c(17.37, 14.23, 32.17, 29.58, 29.88, 31.80),
  sd1 = c(15.50, 12.21, 32.60, 29.44, 32.95, 30.18),
  mean2 = c(10.8, 7.2, 8.6, 10.9, 11.7, 8.7),
  sd2 = c(4.8, 3.3, 3.7, 7.7, 5.2, 3.7)
)

# Function to calculate Cohen's d using two different methods (pooled sd and reading group sd)
calculate_cohens_d <- function(m1, sd1, m2, sd2, method = "pooled") {
  if (method == "pooled") {
    sd_pooled <- sqrt((sd1^2 + sd2^2) / 2)
    d <- (m1 - m2) / sd_pooled
  } else if (method == "reading_sd") {
    d <- (m1 - m2) / sd2
  }
  return(d)
}

# Calculate Cohen's d for each measure using both methods
data$cohens_d_pooled <- mapply(calculate_cohens_d, data$mean1, data$sd1, data$mean2, data$sd2, MoreArgs = list(method = "pooled"))
data$cohens_d_reading_sd <- mapply(calculate_cohens_d, data$mean1, data$sd1, data$mean2, data$sd2, MoreArgs = list(method = "reading_sd"))

# Print results
print(data)


#Create Figure S1
png(filename = "Figure S1.png", width = 2400, height = 1600, res = 300)

par(mfrow=c(2,3),mar = c(4, 4, 4, 1))
fit_gamma_and_plot(Pre2$CRT1_t, 17.37, observed_value=(results2[,1]$mean), cohen_d_value = 1.4, title="CRT 1")
fit_gamma_and_plot(Pre2$CRT2_t, 14.23, observed_value=(results2[,2]$mean), cohen_d_value = 2.1, title="CRT 2")
fit_gamma_and_plot(Pre2$CRT3_t, 32.17, observed_value=(results2[,3]$mean), cohen_d_value = 6.4, title="CRT 3")
fit_gamma_and_plot(Pre2$CRT4_t, 29.58, observed_value=(results2[,4]$mean), cohen_d_value = 2.4, title="CRT 4")
fit_gamma_and_plot(Pre2$CRT5_t, 29.88, observed_value=(results2[,5]$mean), cohen_d_value = 3.5, title="CRT 5")
fit_gamma_and_plot(Pre2$CRT6_t, 31.80, observed_value=(results2[,6]$mean), cohen_d_value = 6.2, title="CRT 6")

dev.off()






