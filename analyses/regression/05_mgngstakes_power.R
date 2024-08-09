#!/usr/bin/env Rscript
# ============================================================================ #
## 05_mgngstakes_power.R
## MGNGStakes estimate post-hoc power for responses and RTs.
# Johannes Algermissen, 2023.

rm(list = ls())

# ============================================================================ #
#### 00) Set directories, load packages and custom functions: ####

## Set codeDir:
codeDir   <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/regression/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs.R")) # Load packages and options settings
dirs <- set_dirs(rootDir)

## Load packages:
library(lme4)
library(pwr)

## Load custom functions:
source(paste0(dirs$codeDir, "functions/00_mgngstakes_functions_regression.R")) # Load functions

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### 01) Read behavior: ####

## Find all raw data files:
data <- read_behavior()

## Preprocessing:
data <- wrapper_preprocessing(data)
cat("Finished pre-processing :-)\n"); # beep()

## Select data, standardize variables:
modData <- subset(data, !(subject_n %in% c(52)))

# ============================================================================ #
#### 02) Estimate ICC: ####

# ----------------------------------------------------- #
### responses:
formula <- "response_n ~ 1 + (1|subject_f)"
mod <- glmer(formula = formula, data = modData, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

## Extract variances and compute ICC:
varData <- as.data.frame(VarCorr(mod))
ICC <- varData$vcov[1] / (varData$vcov[1] + (pi^2)/3); ICC
# Residual variance is fixed to pi^2/3, see 
# https://stats.stackexchange.com/questions/62770/calculating-icc-for-random-effects-logistic-regression
# https://stats.stackexchange.com/questions/128750/residual-variance-for-glmer

## response: 0.04287225

# ----------------------------------------------------- #
### RTs:

formula <- "RTcleaned_n ~ 1 + (1|subject_f)"
mod <- lmer(formula = formula, data = modData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

## Extract variances and compute ICC:
varData <- as.data.frame(VarCorr(mod))
ICC <- varData$vcov[1] / (varData$vcov[1] + varData$vcov[2]); ICC

# RTs: 0.09446111

# ============================================================================ #
#### 03) Compute r achieved with 80% given FIXED number of subjects: ####

## Fixed:
# ICC <- 0.04287225
# ICC <- 0.09446111

## Compute nTrials:
nCue <- 4
nRep <- 20
nTrialBlock <- nCue * nRep
nBlock <- 4
nTrial <- nBlock * nTrialBlock; nTrial

## Compute effective sample size:
nSub <- 54
nData <- nSub*nTrial; nData
Neff <- nData/ (1 + (nSub - 1) * ICC); Neff

## Estimate power:
power <- 0.80
pwr.r.test(n = Neff, sig.level = 0.05, power = power)

## responses: 0.03852042
## RTs: 0.05223095

# END OF FILE.
