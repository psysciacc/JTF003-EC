rm(list=ls())

library(tidyverse)
library(qualtRics)
library(sjPlot)
library(ordinal)
library(lme4)
library(patchwork)

# the data in `manip_test_noOverlap_anonim.csv` have been already anonimized (`PROLIFIC_PID` has been replaced with an anonymous code) and do not contain the participants that also did the reading pre-test.
d <- read_csv("./manip_test_noOverlap_anonim.csv")