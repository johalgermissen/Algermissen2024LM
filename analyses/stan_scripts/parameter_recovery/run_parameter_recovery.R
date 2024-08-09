#!/usr/bin/env Rscript
# ======================================================================================================== #
# run_parameter_recovery.R
# Perform parameter recovery by sampling parameters, simulating data, fitting model to them using Stan.
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
# model2gen <- 4
# suffix <- "_sepParam"
# nIter <- 1000

# ======================================================================================================== #
#### 00b) Load libraries: ####

library(boot) # for inv.logit
library(stringr) # for str_pad
library(mvtnorm) # for rmvnorm

library(rstan) # for reading posterior
library(loo) # for extract

library(RWiener) # for simulating with RWiener

library(lattice)
library(pastecs)

# ======================================================================================================== #
#### 00c) Convert model ID to string: ####

model2gen_str <- str_pad(model2gen, width = 2, side = "left", pad = "0")

## Check inputs:
cat(paste0("Model to fit is ", model2gen_str), "\n")
cat(paste0("Suffix is ", suffix, "\n"))
cat(paste0("Perform parameter recovery over ", nIter, " iterations\n"))

# ======================================================================================================== #
#### 00d) Set directories: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)

targetDir <- paste0(dirs$paramRecov, "subParMvtnorm", nIter, "/")
dir.create(targetDir, showWarnings = FALSE)

# ======================================================================================================== #
#### 00e) Load helper files: ####

cat("Load helper scripts\n")

source(paste0(dirs$simscripts, "RLDDM_modSims", suffix, "_flat.R")) # model code for simulations
source(paste0(dirs$models, "RLDDM4stan", suffix, "_flat.R")) # model code for simulations
source(paste0(dirs$helpers, "initVals4stan", suffix, "_flat.R")) # model code for simulations

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
      cat(paste0("Transform parameter ", subParamNames[iCol], " using og1p_exp\n"))
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
# dList$subject <- rep(1:nSub, each = nTrial)
# dList$trialnr <- rep(1:nTrial, times = nSub)

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
#### 02) Load previously fitted model, extract subject-level parameters: ####

cat("Load and prepare parameters\n")

## Set model manually ÄGAIN:
# model2gen <- 4
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
allParNames <- names(posterior) # names of all parameters
subParamNames <- allParNames[grepl("sub_", allParNames, fixed = T)] # names of subject level parameters
nParamSub <- length(subParamNames) # number subject level parameters

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

# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
#### 03) Determine ground-truth parameter values:

# ======================================================================================================== #
#### 03a) Option A: Capture parameter correlations across subjects, sample from multivariate normal distribution: ####

## 1) Mean per subject across samples:
subParMat <- apply(subParSampleMat, c(1, 2), mean) # average over samples

## 2) Transform to sampling space (unbounded):
subParMatTransf <- inverseTransform(subParMat)

## 3) Compute means and covariance across subjects:
meanVec <- apply(subParMatTransf, c(2), mean); meanVec
sdVec <- apply(subParMatTransf, c(2), sd); sdVec
covMat <- cov(subParMatTransf); covMat

## 4) Sample new parameters from multivariate normal:
# # nIter <- 100 # number samples to draw
set.seed(70) # set RNG seed
sampledIterParMatTransf <- rmvnorm(n = nIter, mean = meanVec, sigma = covMat)

## 5) Now undo transform so everything is in bounded space:
sampledIterParMat <- forwardTransform(sampledIterParMatTransf)

# ======================================================================================================== #
#### 03b) Save parameters: ####

fileName <- paste0("M", model2gen_str, suffix, "_parameterRecovery_trueParams.csv")
write.csv(sampledIterParMat, paste0(targetDir, fileName), row.names = F)

# ======================================================================================================== #
#### 03c) Load code to simulate from model: ####

run_simulations <- eval(parse(text = paste0("modSim_M", model2gen_str)))

# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
#### 04) Prepare fitting model: ####

model2fit <- model2gen

## Model index as string:  
model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")

## Model code:
# source(paste0(dirs$models, "RLDDM4stan", suffix, "_flat.R")) # model code for simulations
mstr <- get(paste("mstr", model2fit_str, sep = "")) # string to be used as input for stan fitting

## Parameter initialization:
# source(paste0(dirs$helpers, "initVals4stan", suffix, "_flat.R")) # model code for simulations
init_fun <- eval(parse(text = paste0("init_fun_mod", model2fit_str)))

# ------------------------------------------------------------------- #
### Input and output parameters:

parType = c("alpha", "tau", "beta", "betaWin", "betaAvoid", "deltaInt", "deltaSlope", "deltaWin", "deltaAvoid", "epsilon", "pi", "theta") 
paramOut4Mod <- list(c(1,2,3,6), c(1,2,3,6,7,10), c(1,2,4,5,6,7,10), c(1,2,3,7,8,9,10),
                     c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11),
                     c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), 
                     c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12)
)

## Select parameters for selected model:
paramIdx <- paramOut4Mod[[model2fit]] # extract indices for selected model
paramInit <- parType[paramIdx]
paramOut <- c(paramInit, "log_lik") # add log likelihood
paramInit <- paste0("nd_", paramInit) # add nd_

cat(paste0("Input  parameters are ", paste0(paramInit, collapse = ", "), "\n"))
cat(paste0("Output parameters are ", paste0(paramOut, collapse = ", "), "\n"))

# ------------------------------------------------------------------- #
### Fitting settings:

nStanIter <- 1000
nWarmup <- floor(nStanIter/2)
nChain <- 4
cat(paste0("Fit each model with ", nStanIter, " iterations (", nWarmup, " warmup) over ", nChain, " chains\n"))

# ======================================================================================================== #
#### 05) Initialize output arrays for saving: ####

cat("Initialize output objects\n")

## Initialize objects with summary statistics of each model fit:
fittedIterParMat <- matrix(NA, nIter, nParamSub)
neffMat <- matrix(NA, nIter, nParamSub)
rhatMat <- matrix(NA, nIter, nParamSub)
waicVec <- rep(NA, nIter)
looVec <- rep(NA, nIter)
durationVec <- rep(NA, nIter)
divergentVec <- rep(NA, nIter)
maxTreeVec <- rep(NA, nIter)
energyMat <- matrix(NA, nIter, nChain)

# --------------------------------------------------------------------- #
## Start simulations:

cat("Start simulations\n")
start.time.loop <- Sys.time()

for (iIter in 1:nIter){ # iIter <- 1
  
  # ------------------------------------------------------------------- #
  cat("* ----------------------------------------------------------------------------- *\n")
  cat(paste0(">>> Iteration ", iIter, ": Start\n"))
  
  # ------------------------------------------------------------------- #
  ## Extract ground-truth parameters:
  
  par <- sampledIterParMat[iIter, ]
  
  # ------------------------------------------------------------------- #
  #### 05a) Simulate data for single subject: ####
  
  cat(paste0(">>> Iteration ", iIter, ": Simulate data\n"))
  dListSub <- run_simulations(data, par) # takes around 16 sec. for 54 people
    
  # ------------------------------------------------------------------- #
  #### 05b) Fit model: ####
  
  seed = 70 # arbitrary, but setting a fixed seed makes models reproducible
  
  ## Fit, track duration:
  cat(paste0(">>> Iteration ", iIter, ": Start fitting model\n"))
  start.time.model <- Sys.time(); cat("Start at", as.character(start.time.model), "\n")
  fit <- stan(model_code = mstr, data = dListSub, init = init_fun,
              warmup = nWarmup, iter = nStanIter, chains = nChain, cores = nChain,
              seed = seed+1, pars = paramOut)
  end.time.model <- Sys.time(); cat("Finish at", as.character(end.time.model), "\n")
  duration <- difftime(end.time.model, start.time.model); cat(paste0("\nFitting iteration ", iIter, " took ", capture.output(duration)), "\n")
  durationVec[iIter] <- duration
  ## Fit DDM with 4 chains, 1000 samples each, in 1.5 minutes
  
  # ------------------------------------------------------------------- #
  #### 05c) Extract and store output: ####
  
  ## Extract from posterior:
  cat(paste0(">>> Iteration ", iIter, ": Save model output\n"))
  posterior <- rstan::extract(fit, inc_warmup = FALSE, permuted = TRUE) # collapse across chains
  parAvgTable <- summary(fit)$summary # extract summary statistics
  
  ## Parameters:
  fittedIterParMat[iIter, ] <- parAvgTable[1:nParamSub, 1]
  
  ## Convergence criteria:
  neffMat[iIter, ] <- parAvgTable[1:nParamSub, 9]
  rhatMat[iIter, ] <- parAvgTable[1:nParamSub, 10]
  
  ## Model fit indices:
  waicVec[iIter] <- loo::waic(extract_log_lik(fit))[[1]][3, 1]
  looVec[iIter] <- loo::loo(fit)[[1]][3, 1]
  
  ## Divergent transitions:
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  divergentVec[iIter] <- sum(divergent)
  
  ## Maximum tree depth reached:
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  max_depth <- 10
  maxTreeVec[iIter] <- length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  
  ## Energy:
  for (iChain in 1:length(sampler_params)) { # loop over chains # n <- 1
    energies = sampler_params[iChain][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    energyMat[iIter, iChain] <- numer/denom
  } 
  
  # ======================================================================================================== #
  #### 06a) SAVE fitted parameters: ####
  
  # paramRecovFile <- paste0(targetDir, "M", model2gen_str, suffix, "_parameterRecovery.Rdata")
  paramRecovFile <- paste0(targetDir, "M", model2gen_str, suffix, "_parameterRecovery.Rdata")
  cat(paste0("Start saving ", paramRecovFile, " for iteration ", iIter, "\n"))
  
  ## Save as R workspace:
  save(sampledIterParMat, fittedIterParMat, 
       neffMat, rhatMat, waicVec, looVec, durationVec,
       divergentVec, maxTreeVec, energyMat,
       paramInit, paramOut, subParamNames, nParamSub,
       nStanIter, nWarmup, nChain,
       file = paramRecovFile)
  cat(paste0("Finished saving ", paramRecovFile, " for iteration ", iIter, "! :-)\n"))
  
}
cat("* ----------------------------------------------------------------------------- *\n")
cat("* ----------------------------------------------------------------------------- *\n")
cat("* ----------------------------------------------------------------------------- *\n")
cat("Finished simulations!!! :-)\n")
end.time.loop <- Sys.time()
duration <- difftime(end.time.loop, start.time.loop)
cat(paste0("\nEntire loop of ", nIter, " iterations took ", capture.output(duration), "\n"))

cat("END OF SCRIPT\n")

# END OF FILE.
