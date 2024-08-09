#!/usr/bin/env Rscript
# ======================================================================================================== #
# run_model_recovery_fit.R
# Load simulate data sets and fit model with rstan.
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

model2fit <- as.numeric(args[3])
if (is.na(model2fit)){
  model2fit <- 1
  cat("Errorneous input argument, set model2fit to", model2fit, "\n")
}

## Manual:
# model2gen <- 12
# suffix <- "_sepParam"
# model2fit <- 4

# ======================================================================================================== #
#### 00b) Load libraries: ####

library(boot) # for inv.logit
library(stringr) # for str_pad

library(rstan) # for reading posterior
library(loo) # for extract

library(lattice)
library(pastecs)

# ======================================================================================================== #
#### 00c) Convert model ID to string: ####

# model2gen <- 12
model2gen_str <- str_pad(model2gen, width = 2, side = "left", pad = "0")
# model2fit <- 12
model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")

## Check inputs:
cat(paste0("Model to use generated data from is M", model2gen_str), "\n")
cat(paste0("Suffix is ", suffix, "\n"))
cat(paste0("Model to fit is M", model2fit_str), "\n")

# ======================================================================================================== #
#### 00d) Set directories: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)

targetDir <- paste0(dirs$modelRecov, "modelRecov_nIter1000/")
dir.create(targetDir, showWarnings = FALSE)
cat(paste0("Read simulated data from ", targetDir, " \n"))

# ======================================================================================================== #
#### 00e) Load helper files: ####

cat("Load helper scripts\n")

source(paste0(dirs$simscripts, "RLDDM_modSims", suffix, "_flat.R")) # model code for simulations
source(paste0(dirs$models, "RLDDM4stan", suffix, "_flat.R")) # model code for simulations
source(paste0(dirs$helpers, "initVals4stan", suffix, "_flat.R")) # model code for simulations

# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
#### 01) Load data: ####

## File name:
dataFile <- paste0(targetDir, "M", model2gen_str, suffix, "_modelRecovery_data.Rdata")
cat(paste0("Load data from file ", dataFile, "\n"))

## Load as R workspace:
load(dataFile)

## Determine number iterations:
nIter <- length(dListIter)
nParamInput <- ncol(sampledIterParMat)
cat(paste0("Found data for ", nIter, " iterations\n"))

# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
#### 02) Prepare fitting model: ####

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
nParamOutput <- length(paramOut) - 1 # without log_lik

cat(paste0("Input  parameters are ", paste0(paramInit, collapse = ", "), "\n"))
cat(paste0("Output parameters are ", paste0(paramOut, collapse = ", "), "\n"))

# ------------------------------------------------------------------- #
### Fitting settings:

nStanIter <- 1000
nWarmup <- floor(nStanIter/2)
nChain <- 4
cat(paste0("Fit each model with ", nStanIter, " iterations (", nWarmup, " warmup) over ", nChain, " chains\n"))

# ------------------------------------------------------------------- #
### Output file Name:

# outputFile <- paste0(targetDir, "genM", model2gen_str, "_fitM", model2fit_str, suffix, "_modelRecovery_fit.Rdata"
outputFile <- paste0(targetDir, "genM", model2gen_str, "_fitM", model2fit_str, suffix, "_modelRecovery_fit.Rdata")

# ======================================================================================================== #
#### 03) Initialize output arrays for saving: ####

cat("Initialize output objects\n")

## Initialize objects with summary statistics of each model fit:
fittedIterParMat <- matrix(NA, nIter, nParamOutput)
neffMat <- matrix(NA, nIter, nParamOutput)
rhatMat <- matrix(NA, nIter, nParamOutput)
waicVec <- rep(NA, nIter)
looVec <- rep(NA, nIter)
durationVec <- rep(NA, nIter)
divergentVec <- rep(NA, nIter)
maxTreeVec <- rep(NA, nIter)
energyMat <- matrix(NA, nIter, nChain)

# --------------------------------------------------------------------- #
### Check if file exists, if yes, load and continue iterations:

startIter <- 1

# outputFile <- "/project/2420133.01/MGNGStakes/results/stan_results/model_recovery/modelRecov_nIter1000/genM03_fitM09_sepParam_modelRecovery_fit.Rdata"
if (file.exists(outputFile)){
  cat(paste0("File ", outputFile, " already exists; load and continue:\n"))
  load(outputFile)
  startIter <- sum(complete.cases(waicVec)) + 1 # number complete cases + 1
  cat(paste0("Continue from iteration ", startIter, "\n"))
  if(startIter > nIter){
    cat("startIter > nIter, already all simulations run!\n")
    stop("startIter > nIter, already all simulations run!")
    }
}

# --------------------------------------------------------------------- #
### Start fitting:

cat("Start fitting\n")
start.time.loop <- Sys.time()
for (iIter in startIter:nIter){ # iIter <- 1
  
  # ------------------------------------------------------------------- #
  cat("* ----------------------------------------------------------------------------- *\n")
  cat(paste0(">>> Iteration ", iIter, ": Start\n"))
  
  # ------------------------------------------------------------------- #
  #### 04) Extract data: ####
  
  cat(paste0(">>> Iteration ", iIter, ": Extract data\n"))
  dListSub <- dListIter[[iIter]] # takes around 16 sec. for 54 people
  
  ## Inspect:
  # str(dListSub)
  # table(dListSub$resp)
  # densityplot(dListSub$rt)
  
  # ------------------------------------------------------------------- #
  #### 05) Fit model: ####
  
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
  fittedIterParMat[iIter, ] <- parAvgTable[1:nParamOutput, 1]
  
  ## Convergence criteria:
  neffMat[iIter, ] <- parAvgTable[1:nParamOutput, 9]
  rhatMat[iIter, ] <- parAvgTable[1:nParamOutput, 10]
  
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
  
  cat(paste0("Start saving ", outputFile, " for iteration ", iIter, "\n"))
  
  ## Save as R workspace:
  save(sampledIterParMat, fittedIterParMat, dListIter,
       neffMat, rhatMat, waicVec, looVec, durationVec,
       divergentVec, maxTreeVec, energyMat,
       paramInit, paramOut, subParamNames, nParamSub, nParamInput, 
       paramInit, paramOut, nParamOutput,
       nStanIter, nWarmup, nChain,
       file = outputFile)
  cat(paste0("Finished saving ", outputFile, " for iteration ", iIter, "! :-)\n"))
  
}
cat("* ----------------------------------------------------------------------------- *\n")
cat("* ----------------------------------------------------------------------------- *\n")
cat("* ----------------------------------------------------------------------------- *\n")
cat("Finished simulations!!! :-)\n")
end.time.loop <- Sys.time()
duration <- difftime(end.time.loop, start.time.loop)
cat(paste0("\nEntire loop of ", nIter, " iterations took ", capture.output(duration), "\n"))
warning(paste0("Finished running and saving all ", nIter, " iterations, reached end of script\n"))
cat("END OF SCRIPT\n")

# END OF FILE.
