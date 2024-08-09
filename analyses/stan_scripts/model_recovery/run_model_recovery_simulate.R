#!/usr/bin/env Rscript
# ======================================================================================================== #
# run_model_recovery_simulate.R
# Simulate data sets for model recovery.
# Johannes Algermissen, 2023.
# ======================================================================================================== #
#### 00a) Retrieve inputs from bash: ####

### Clear workspace and console.
cat("\014")
rm(list = ls()); 

### Retrieve inputs from bash:
args <- commandArgs(trailingOnly = TRUE)

### Extract info:
## Model to fit:
model2gen <- as.numeric(args[1])
if (is.na(model2gen)){
  model2gen <- 1
  cat("Errorneous input argument, set model2gen to", model2gen, "\n")
}

suffix <- args[2]
if (is.na(suffix)){
  suffix <- "_sepParam" # _standard, _deltaIntSlope, _sepParam
  cat("Not suffix provided, set to",suffix, "\n")
}

nIter <- as.numeric(args[3])
if (is.na(nIter)){
  nIter <- 10
  cat("Errorneous input argument, set nIter to", nIter, "\n")
}

## Manual:
# model2gen <- 10
# suffix <- "_sepParam"
# nIter <- 1000

# ======================================================================================================== #
#### 00b) Load libraries: ####

library(boot) # for inv.logit
library(stringr) # for str_pad
library(mvtnorm) # for rmvnorm

## Load libraries:
# install.packages("rstan", dependencies = TRUE)
library(rstan) # for reading posterior
# if (!require("rstan", character.only = TRUE)){
#   require("rstan"); 
#   cat("Load stan library\n")
# }
library(loo) # for extract

library(RWiener) # for simulating with RWiener

library(lattice)
library(pastecs)

# ======================================================================================================== #
#### 00c) Convert model ID to string: ####

# model2gen <- 12
model2gen_str <- str_pad(model2gen, width = 2, side = "left", pad = "0")

## Check inputs:
cat(paste0("Model to generate data from is ", model2gen_str), "\n")
cat(paste0("Suffix is ", suffix, "\n"))
cat(paste0("Simulate data over ", nIter, " iterations\n"))

# ======================================================================================================== #
#### 00d) Set directories: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)

targetDir <- paste0(dirs$modelRecov, "modelRecov_nIter", nIter, "/")
dir.create(targetDir, showWarnings = FALSE)

# ======================================================================================================== #
#### 00e) Initialize functions for parameter transformation: ####

## Transformation functions:
log1p_exp <- function(x){return(log(1 + exp(x)))}
log_exp_m1 <- function(x){return(log(exp(x) - 1))}
logit <- function(x){return(log(x / (1 - x)))}
inv_logit <- function(x){return(1/(1+exp(-x)))}

forwardTransform <- function(mat){
  ## Loop over parameters, transform:
  nCol <- ncol(mat)
  for (iCol in 1:nCol){
    if (grepl("alpha", subParamNames[iCol])){
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log1p_exp\n"))
      mat[, iCol] <- log1p_exp(mat[, iCol])
    }
    if (grepl("tau", subParamNames[iCol])){
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log1p_exp\n"))
      mat[, iCol] <- log1p_exp(mat[, iCol])
    }
    if (grepl("beta", subParamNames[iCol])){
      cat(paste0("Transform parameter ", subParamNames[iCol], " using inv.logit\n"))
      mat[, iCol] <- inv_logit(mat[, iCol])
    }
    if (grepl("deltaSlope", subParamNames[iCol])){
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log1p_exp\n"))
      mat[, iCol] <- log1p_exp(mat[, iCol])
    }
    if (grepl("epsilon", subParamNames[iCol])){
      cat(paste0("Transform parameter ", subParamNames[iCol], " using inv.logit\n"))
      mat[, iCol] <- inv_logit(mat[, iCol])
    }
    if (grepl("pi", subParamNames[iCol]) & model2gen %in% c(5, 6, 9, 10, 11, 12)){
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log1p_exp\n"))
      mat[, iCol] <- log1p_exp(mat[, iCol])
    }
    if (grepl("pi", subParamNames[iCol]) & model2gen %in% c(7)){
      cat(paste0("Transform parameter ", subParamNames[iCol], " using inv.logit\n"))
      mat[, iCol] <- inv_logit(mat[, iCol])
    }
    if (grepl("theta", subParamNames[iCol]) & model2gen %in% c(9, 12)){
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log1p_exp\n"))
      mat[, iCol] <- log1p_exp(mat[, iCol])
    }
  }
  return(mat)
}

inverseTransform <- function(mat){
  ## Loop over parameters, transform:
  nCol <- ncol(mat)
  for (iCol in 1:nCol){
    if (grepl("alpha", subParamNames[iCol])){
      mat[which(mat[, iCol] < 0), iCol] <- 0 # set to minimum 0
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log_exp_m1\n"))
      mat[, iCol] <- log_exp_m1(mat[, iCol])
    }
    if (grepl("tau", subParamNames[iCol])){
      mat[which(mat[, iCol] < 0), iCol] <- 0 # set to minimum 0
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log_exp_m1\n"))
      mat[, iCol] <- log_exp_m1(mat[, iCol])
    }
    if (grepl("beta", subParamNames[iCol])){
      mat[which(mat[, iCol] < 0), iCol] <- 0 # set to minimum 0
      cat(paste0("Transform parameter ", subParamNames[iCol], " using logit\n"))
      mat[, iCol] <- logit(mat[, iCol])
    }
    if (grepl("deltaSlope", subParamNames[iCol])){
      mat[which(mat[, iCol] < 0), iCol] <- 0 # set to minimum 0
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log_exp_m1\n"))
      mat[, iCol] <- log_exp_m1(mat[, iCol])
    }
    if (grepl("epsilon", subParamNames[iCol])){
      mat[which(mat[, iCol] < 0), iCol] <- 0 # set to minimum 0
      cat(paste0("Transform parameter ", subParamNames[iCol], " using logit\n"))
      mat[, iCol] <- logit(mat[, iCol])
    }
    if (grepl("pi", subParamNames[iCol]) & model2gen %in% c(5, 6, 9, 10, 11, 12)){
      mat[which(mat[, iCol] < 0), iCol] <- 0 # set to minimum 0
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log_exp_m1\n"))
      mat[, iCol] <- log_exp_m1(mat[, iCol])
    }
    if (grepl("pi", subParamNames[iCol]) & model2gen %in% c(7)){
      mat[which(mat[, iCol] < 0), iCol] <- 0 # set to minimum 0
      cat(paste0("Transform parameter ", subParamNames[iCol], " using logit\n"))
      mat[, iCol] <- logit(mat[, iCol])
    }
    if (grepl("theta", subParamNames[iCol]) & model2gen %in% c(9, 12)){
      mat[which(mat[, iCol] < 0), iCol] <- 0 # set to minimum 0
      cat(paste0("Transform parameter ", subParamNames[iCol], " using log_exp_m1\n"))
      mat[, iCol] <- log_exp_m1(mat[, iCol])
    }
  }
  return(mat)
}

# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
#### 01) Load and prepare data of single subject: ####

cat("Load and prepare data\n")

## Load data:
load(paste0(dirs$results, "data4stan.Rdata"))
str(dList)

cat("Data successfully loaded\n")

## Extract data format:
nSub <- dList$nSub
nTrial <- dList$nTrial
nData <- dList$nData

nStim <- length(unique(dList$stimuli))
nResp <- length(unique(dList$resp))

## Select single subject:
iSub <- 1
cat(paste0("Select design structure for subject ", iSub, "\n"))
subIdx <- which(dList$subject == iSub)

## Extract relevant variables:
data <- as.data.frame(cbind(dList$subject[subIdx], dList$trialnr[subIdx],
                            dList$stimuli[subIdx], dList$reqAction[subIdx], dList$valence[subIdx], dList$congruency[subIdx], 
                            dList$stakes[subIdx], dList$validity[subIdx], 
                            dList$resp[subIdx], dList$rt[subIdx], dList$outcome[subIdx]))

## Set variable names:
names(data) <- c("subject", "trialnr", "stimuli", "reqAction", "valence", "congruency", "stakes", "validity",
                 "resp", "rt", "outcome")

# ======================================================================================================== #
#### 02) Load PREVIOUSLY FITTED model, extract subject-level parameters: ####

cat(paste0("Load parameters from model M", model2gen_str, "\n"))

## Set model manually ÄGAIN:
# model2gen <- 12
# model2gen_str <- str_pad(model2gen, width = 2, side = "left", pad = "0")

## Load parameter file:
load(paste0(dirs$results, "M", model2gen_str, suffix, "_fitted.Rdata"))
# load(paste0(dirs$results, "fit_2023_08_15/M", model2gen_str, suffix, "_fitted.Rdata"))

cat("Parameters successfully loaded\n")

## Extract relevant data dimensions: 
nSample <- dim(f2)[1]; nChain <- dim(f2)[2]; nParamAll <- dim(f2)[3]
cat(paste0("Model has ", nSample, " samples per chain, ", nChain, " chains, ", nParamAll, " parameters\n")) # nSamples, nChains, nParameters

## Extract posterior:
cat("Collapse across chains: Create list per parameter with nSample elements \n")
posterior <- rstan::extract(f2, inc_warmup = FALSE, permuted = TRUE) # 

## Identify subject-level parameter names:
allparamNames <- names(posterior) # names of all parameters
subParamNames <- allparamNames[grepl("sub_", allparamNames, fixed = T)] # names of subject level parameters
paramNames <- gsub("sub_", "", subParamNames) # remove sub
nParamSub <- length(paramNames) # number subject level parameters

cat(paste0("Extract ", nParamSub, " subject parameters: ", paste(subParamNames, collapse = ", ")), "\n")

# -------------------------------------------------------------------------- #
### Display results of group-level parameters:

## Plot means:
summary(f2)$summary[1:nParamSub, ]
round(summary(f2)$summary[1:nParamSub, 1], 2) # means
round(summary(f2)$summary[1:nParamSub, 3], 2) # sd

# -------------------------------------------------------------------------- #
## Cast in one big array of dimensions subject x parameters x samples:

cat("Cast subject parameters (per subject, parameter, sample) into matrix \n")
subParSampleMat <- array(NA, c(nSub, nParamSub, nSample*nChain)) # subjects, parameters, iterations
for (iParam in 1:nParamSub){ # iParam <- 1 # for each parameter
  subParSampleMat[, iParam, ] <- t(posterior[[subParamNames[iParam]]]) # extract lists of all subjects for this parameter, transpose, save in matrix
}
cat(paste0("Extracted ", nParamSub, " parameters for ", nSub, " subjects, ", nSample*nChain, " iterations\n"))

## Inspect selected iteration for selected subject:
# subParamNames
# iSub <- 1; iIter <- 1
# subParSampleMat[iSub, , iIter]

# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
#### 03) Determine ground-truth parameter values:

# ======================================================================================================== #
#### 03a) Option A: Capture parameter correlations across subjects, sample from multivariate normal distribution: ####

## 1) Mean per subject across samples:
cat("Compute mean parameter value per subject\n")
subParMat <- apply(subParSampleMat, c(1, 2), mean) # average over samples

## 2) Transform to sampling space (unbounded):
cat("Transform mean parameter values per subject into sampling space\n")
subParMatTransf <- inverseTransform(subParMat)

## 3) Compute means and covariance across subjects:
cat("Compute mean and covariance of parameter values across subjects\n")
meanVec <- apply(subParMatTransf, c(2), mean); meanVec
sdVec <- apply(subParMatTransf, c(2), sd); sdVec
covMat <- cov(subParMatTransf); covMat

## 4) Sample new parameters from multivariate normal:
nIterSample <- nIter * 1000 # number samples to draw
cat(paste0("Sample ", nIterSample, " parameter combinations for model M", model2gen_str, "\n"))
set.seed(70) # set RNG seed
sampledIterParMatTransf <- rmvnorm(n = nIterSample, mean = meanVec, sigma = covMat)

## 5) Now undo transform so everything is in bounded space:
cat("Transform parameters back to normal space\n")
sampledIterParMat <- forwardTransform(sampledIterParMatTransf)

# ======================================================================================================== #
#### 03b) Enforce constraints: ####

### Initialize:
validIdx <- rep(1, nIterSample)

# --------------------------------------------------------- #
### a) Check that learning rate sufficiently high:
#         0%        10%        20%        30%        40%        50%        60%        70%        80%        90%       100% 
# 0.01351338 0.05372806 0.07250291 0.08769763 0.10526289 0.12442260 0.14506835 0.17109247 0.20692901 0.26498964 0.55015058 
if(any(grepl("epsilon", paramNames))){
  threshVal <- 0.05
  cat(paste0("Enforce that epsilon > ", threshVal, "\n"))
  epsIdx <- which(grepl("epsilon", paramNames))
  # densityplot(sampledIterParMat[, epsIdx])
  # quantile(sampledIterParMat[, epsIdx], seq(0, 1, 0.10), na.rm = T)
  validIdx <- validIdx * as.numeric(sampledIterParMat[, epsIdx] > threshVal) # keep only sufficiently high learning rates
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

# --------------------------------------------------------- #
## b) Check that betas sufficiently far apart:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.072  0.054  0.084  0.107  0.129  0.148  0.166  0.188  0.215  0.259  0.447 
if(any(grepl("betaWin", paramNames))){
  threshVal <- 0.20 # 0.05
  cat(paste0("Enforce that betaWin - betaAvoid (starting point bias) > ", threshVal, "\n"))
  betaWinIdx <- which(grepl("betaWin", paramNames))
  betaAvoidIdx <- which(grepl("betaAvoid", paramNames))
  # densityplot(sampledIterParMat[, betaWinIdx] - sampledIterParMat[, betaAvoidIdx])
  # round(quantile(sampledIterParMat[, betaWinIdx] - sampledIterParMat[, betaAvoidIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric((sampledIterParMat[, betaWinIdx] - sampledIterParMat[, betaAvoidIdx]) > threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

# --------------------------------------------------------- #
### c) that deltas sufficiently far apart:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.958  0.497  0.776  0.982  1.146  1.316  1.493  1.675  1.880  2.172  3.150 
if(any(grepl("deltaWin", paramNames))){
  threshVal <- 1.00 # 0.50
  cat(paste0("Enforce that deltaWin - deltaAvoid (drift rate intercepts) > ", threshVal, "\n"))
  deltaWinIdx <- which(grepl("deltaWin", paramNames))
  deltaAvoidIdx <- which(grepl("deltaAvoid", paramNames))
  # densityplot(sampledIterParMat[, deltaWinIdx] - sampledIterParMat[, deltaAvoidIdx])
  # round(quantile(sampledIterParMat[, deltaWinIdx] - sampledIterParMat[, deltaAvoidIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric((sampledIterParMat[, deltaWinIdx] - sampledIterParMat[, deltaAvoidIdx]) > threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

# --------------------------------------------------------- #
### d) that alphas sufficiently far apart:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.225 -0.026  0.009  0.034  0.055  0.074  0.092  0.114  0.139  0.174  0.346 
if(model2gen %in% c(5, 9, 10)){
  threshVal <- 0.14 # 0.05 0.10
  cat(paste0("Enforce that pi - alpha (threshold) > ", threshVal, "\n"))
  alphaIdx <- which(grepl("alpha", paramNames))
  piIdx <- which(grepl("pi", paramNames))
  # densityplot(sampledIterParMat[, piIdx] - sampledIterParMat[, alphaIdx])
  # round(quantile(sampledIterParMat[, piIdx] - sampledIterParMat[, alphaIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric((sampledIterParMat[, piIdx] - sampledIterParMat[, alphaIdx]) > threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

# --------------------------------------------------------- #
### e) that taus sufficiently far apart:

## Pi minus tau:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.220 -0.045 -0.020 -0.003  0.011  0.023  0.035  0.047  0.061  0.084  0.214 
if(model2gen %in% c(6, 11, 12)){
  threshVal <- 0.07 # 0.023 0.05
  cat(paste0("Enforce that pi - tau > (non-decision time) ", threshVal, "\n"))
  tauIdx <- which(grepl("tau", paramNames))
  piIdx <- which(grepl("pi", paramNames))
  # densityplot(sampledIterParMat[, piIdx] - sampledIterParMat[, tauIdx])
  # round(quantile(sampledIterParMat[, piIdx] - sampledIterParMat[, tauIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric((sampledIterParMat[, piIdx] - sampledIterParMat[, tauIdx]) > threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

## Theta minus tau:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.184 -0.064 -0.040 -0.024 -0.011  0.002  0.015  0.028  0.045  0.071  0.232 
if(model2gen %in% c(9)){
  threshVal <- 0.07 # 0.03 0.06
  cat(paste0("Enforce that abs(theta - tau) (non-decision time) > ", threshVal, "\n"))
  tauIdx <- which(grepl("tau", paramNames))
  thetaIdx <- which(grepl("theta", paramNames))
  # densityplot(sampledIterParMat[, thetaIdx] - sampledIterParMat[, tauIdx])
  # round(quantile(sampledIterParMat[, thetaIdx] - sampledIterParMat[, tauIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric(abs(sampledIterParMat[, thetaIdx] - sampledIterParMat[, tauIdx]) > threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

## Theta minus pi:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.215 -0.080 -0.054 -0.035 -0.019 -0.004  0.012  0.029  0.051  0.081  0.291 
if(model2gen %in% c(12)){
  threshVal <- 0.07 # 0.03 0.06
  cat(paste0("Enforce that abs(theta - pi) (non-decision time) > ", threshVal, "\n"))
  piIdx <- which(grepl("pi", paramNames))
  thetaIdx <- which(grepl("theta", paramNames))
  # densityplot(sampledIterParMat[, thetaIdx] - sampledIterParMat[, piIdx])
  # round(quantile(sampledIterParMat[, thetaIdx] - sampledIterParMat[, piIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric(abs(sampledIterParMat[, thetaIdx] - sampledIterParMat[, piIdx]) > threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

# --------------------------------------------------------- #
### f) that betas sufficiently far apart:

## Theta minus pi:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.088 -0.019 -0.006  0.004  0.012  0.021  0.029  0.037  0.047  0.062  0.143 
if(model2gen %in% c(7)){
  threshVal <- 0.07 # 0.03 0.06
  cat(paste0("Enforce that beta - pi > ", threshVal, "\n"))
  betaIdx <- which(grepl("beta", paramNames))
  piIdx <- which(grepl("pi", paramNames))
  # densityplot(sampledIterParMat[, betaIdx] - sampledIterParMat[, piIdx])
  # round(quantile(sampledIterParMat[, betaIdx] - sampledIterParMat[, piIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric((sampledIterParMat[, betaIdx] - sampledIterParMat[, piIdx]) > threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

# --------------------------------------------------------- #
### g) that drift bonus positive:

## pi:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.164 -0.140 -0.136 -0.133 -0.131 -0.128 -0.126 -0.123 -0.120 -0.116 -0.089 
if(model2gen %in% c(8)){
  threshVal <- -0.14 # -0.10 -0.14
  cat(paste0("Enforce that pi (drift bonus) < ", threshVal, "\n"))
  piIdx <- which(grepl("pi", paramNames))
  # densityplot(sampledIterParMat[, piIdx])
  # round(quantile(sampledIterParMat[, piIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric(sampledIterParMat[, piIdx] < threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}

## theta:
#     0%    10%    20%    30%    40%    50%    60%    70%    80%    90%   100% 
# -0.078 -0.064 -0.061 -0.058 -0.056 -0.054 -0.053 -0.051 -0.048 -0.045 -0.029 
if(model2gen %in% c(10, 11)){
  threshVal <- -0.07 # -0.05 -0.06
  cat(paste0("Enforce that theta (drift bonus) < ", threshVal, "\n"))
  thetaIdx <- which(grepl("theta", paramNames))
  # densityplot(sampledIterParMat[, thetaIdx])
  # round(quantile(sampledIterParMat[, thetaIdx], seq(0, 1, 0.10), na.rm = T), 3)
  validIdx <- validIdx * as.numeric(sampledIterParMat[, thetaIdx] < threshVal) # keep only sufficiently high differences
  cat(paste0("which leaves ", sum(validIdx), " valid rows\n"))
}
cat(paste0("M", model2gen_str, ": after enforcing constraints: ", sum(validIdx), " parameter combinations left! :-)\n"))

# ======================================================================================================== #
#### 03d) Select and save parameters: ####

## Select data:
if (sum(validIdx) < nIter){stop(paste0("Only ", sum(validIdx), " rows left; not enough for ", nIter, " iterations\n"))}
cat(paste0("Final parameter matrix has ", sum(validIdx), " valid rows; select first ", nIter, " of them and save\n"))
cat(paste0("Original parameter matrix has ", nrow(sampledIterParMat), " rows\n"))
sampledIterParMat <- sampledIterParMat[which(validIdx == 1), ]
cat(paste0("Parameter matrix after constraints enforcement has ", nrow(sampledIterParMat), " rows\n"))
sampledIterParMat <- sampledIterParMat[1:nIter, ]
cat(paste0("Final parameter matrix has ", nrow(sampledIterParMat), " rows\n"))

## File name:
fileName <- paste0("M", model2gen_str, suffix, "_modelRecovery_trueParams.csv")

## Save file:
write.csv(sampledIterParMat, paste0(targetDir, fileName), row.names = F)

## Read back in:
# sampledIterParMat <- read.csv(paste0(targetDir, fileName))

# ======================================================================================================== #
#### 03e) Load code to simulate from model: ####

source(paste0(dirs$simscripts, "RLDDM_modSims", suffix, "_flat.R")) # model code for simulations
run_simulations <- eval(parse(text = paste0("modSim_M", model2gen_str)))

# ======================================================================================================== #
#### 04) Simulate data: ####

## Initialize:
dListIter <- list()

cat("Start simulations\n")
start.time.loop <- Sys.time()

for (iIter in 1:nIter){ # iIter <- 1
  
  # ------------------------------------------------------------------- #
  iIter_str <- str_pad(iIter, width = 4, side = "left", pad = "0")

  # ------------------------------------------------------------------- #
  ## Extract ground-truth parameters:
  
  par <- sampledIterParMat[iIter, ]
  
  # ------------------------------------------------------------------- #
  ## Simulate data for single subject:
  
  cat(paste0(">>> Model M", model2gen_str, ", iteration ", iIter_str, ": Simulate data\n"))
  dListSub <- run_simulations(data, par) # takes around 16 sec. for 54 people
  dListSub$rt[dListSub$resp == 0] <- 0 # se RTs for NoGo to 0
  dListIter[[iIter]] <- dListSub # save in long list
    
}
cat("* ----------------------------------------------------------------------------- *\n")
cat("Finished simulations!!! :-)\n")
end.time.loop <- Sys.time()
duration <- difftime(end.time.loop, start.time.loop)
cat(paste0("\nEntire loop of ", nIter, " iterations took ", capture.output(duration), "\n"))

# ======================================================================================================== #
#### 06a) SAVE subject list: ####

## File name:
modelRecovFile <- paste0(targetDir, "M", model2gen_str, suffix, "_modelRecovery_data.Rdata")
cat(paste0("Start saving ", modelRecovFile, "\n"))

## Save as R workspace:
save(sampledIterParMat, dListIter, subParamNames, nParamSub,
     file = modelRecovFile)
cat(paste0("Finished saving ", modelRecovFile, "! :-)\n"))

cat("END OF SCRIPT\n")

# END OF FILE.
