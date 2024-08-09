#!/usr/bin/env Rscript
# ======================================================================================================== #
# run_simulations.R
# Model fitting using rstan.
# Johannes Algermissen, 2023.
# ======================================================================================================== #
#### RETRIEVE INPUTS FROM BASH: ####

### Clear workspace and console.
cat("\014")
rm(list = ls()); 

### Retrieve inputs from bash:
args <- commandArgs(trailingOnly = TRUE)

### Extract info:
## Model to fit:
iMod <- as.numeric(args[1])
if (is.na(iMod)){
  iMod <- 1
  cat("Errorneous input argument, set iMod to", iMod, "\n")
}

suffix <- args[2]
if (is.na(suffix)){
  suffix <- "_sepParam" # _standard, _deltaIntSlope, _sepParam
  cat("Not suffix provided, set to",suffix, "\n")
}


simType <- args[3]
if (is.na(simType)){
  simType <- "modSim" # osap, modSim
  cat("Not simulation type provided, set to",simType, "\n")
}


## Manual:
# iMod <- 5
# suffix <- "_sepParam"
# simType <- "modSim"

# ======================================================================================================== #
#### SET SEED: ####

set.seed(70)

# ======================================================================================================== #
#### LOAD LIBRARIES: ####

library(boot) # for inv.logit
library(stringr) # for str_pad
library(rstan) # for reading posterior

library(RWiener) # for simulating with RWiener

library(lattice)
library(pastecs)

# ======================================================================================================== #
## CONVERT MODEL ID TO STRING:

iMod_str <- str_pad(iMod, width = 2, side = "left", pad = "0")

## Check inputs:
cat(paste0("Model to fit is ", iMod_str), "\n")
cat(paste0("Suffix is ", suffix, "\n"))
cat(paste0("Simulation type is ", simType, "\n"))

# ======================================================================================================== #
#### SET DIRECTORIES: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)

# ======================================================================================================== #
#### LOAD HELPER SCRIPTS: ####

cat("Load helper scripts\n")

source(paste0(dirs$simscripts, "RLDDM_modSims", suffix, ".R")) # model code for simulations

# ======================================================================================================== #
#### LOAD AND PREPARE DATA: ####

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

nRep <- nTrial/nStim

## Extract relevant variables:
data <- as.data.frame(cbind(dList$subject, dList$trialnr,
                            dList$stimuli, dList$reqAction, dList$valence, dList$congruency, dList$stakes, dList$validity, 
                            dList$resp, dList$rt, 
                            dList$outcome))

## Set variable names:
names(data) <- c("subject", "trialnr", "stimuli", "reqAction", "valence", "congruency", "stakes", "validity",
                 "resp", "rt", "outcome")

## Delete RTs for NoGo (resp == 0) (handled as 0 within stan):
data$rt[data$resp == 0] <- NA

# ======================================================================================================== #
#### LOAD MODEL AND FITTED PARAMETERS: ####

cat("Load and prepare parameters\n")

## Load parameter file:
load(paste0(dirs$results, "M", iMod_str, suffix, "_fitted.Rdata"))
fit <- f2

cat("Parameters successfully loaded\n")

## Extract relevant data dimensions: 
nSample <- dim(fit)[1]; nChain <- dim(fit)[2]; nParamAll <- dim(fit)[3]
cat(paste0("Model has ", nSample, " samples per chain, ", nChain, " chains, ", nParamAll, " parameters\n")) # nSamples, nChains, nParameters

## Extract posterior:
cat("Collapse across chains: Create list per parameter with nSample elements \n")
posterior <- rstan::extract(fit, inc_warmup = FALSE, permuted = TRUE) # 

## Extract subject-level parameter names:
allParNames <- names(posterior) # names of all parameters
subParNames <- allParNames[grepl("sub_", allParNames, fixed = T)] # names of subject level parameters
nParamSub <- length(subParNames) # number subject level parameters

cat(paste0("Use ", nParamSub, " subject parameters: ", paste(subParNames, collapse = ", ")), "\n")

# -------------------------------------------------------------------------- #
## Cast in one big array of dimensions subject x subject parameters x samples:
cat("Cast subject parameters (per subject, parameter, sample) into matrix \n")
allSubParMat <- array(NA, c(nSub, nParamSub, nSample*nChain)) # subjects, parameters, iterations
for (iParam in 1:nParamSub){ # iParam <- 1 # for each parameter
  allSubParMat[, iParam, ] <- t(posterior[[subParNames[iParam]]]) # extract lists of all subjects for this parameter, transpose, save in matrix
}

# ======================================================================================================== #
#### INITIALIZE ARRAYS FOR SAVING: ####

cat("Initialize output objects\n")

## Number of samples in chain:
nSampleAll <- nSample * nChain # based on number of samples in parameter object

## Subset to number of iterations:
nIter <- 1000
selIter <- sample(1:nSampleAll, nIter, replace = F) # draw posterior samples without replacement
selIter <- sort(selIter) # sort iterations

## Initialize objects to extract and store:
simQGo <- array(NaN, dim = c(nIter, nData))
simQNoGo <- array(NaN, dim = c(nIter, nData))
simResp <- array(NaN, dim = c(nIter, nData))
simRT <- array(NaN, dim = c(nIter, nData))

# ======================================================================================================== #
#### RUN SIMULATIONS: ####

## Initialize function for running simulations:
run_simulations <- eval(parse(text = paste0(simType, "_M", iMod_str)))

cat("Start simulations\n")

start.time <- Sys.time()
for (iIter in 1:nIter){ # iIter <- 1
  
  # ------------------------------------------------------------------- #
  cat(paste0("Start iteration ", iIter, "\n"))
  
  # ------------------------------------------------------------------- #
  ## Extract correct parameters:

  par <- allSubParMat[, , selIter[iIter]]
  colnames(par) <- subParNames # add column ways

  # ------------------------------------------------------------------- #
  ## Run simulations:

  tmp <- run_simulations(data, par) # takes around 16 sec. for 54 people
    
  # ------------------------------------------------------------------- #
  ## Extract simulated objects and store in large matrices of dimension nIter x nData:
  
  simQGo[iIter, ] <- tmp$QGo
  simQNoGo[iIter, ] <- tmp$QNoGo
  simResp[iIter, ] <- tmp$sim_resp
  simRT[iIter, ] <- tmp$sim_rt
  
  ## Delete NoGo RTs: 
  idx <- simResp[iIter, ] == 0 # localize NoGos
  simRT[iIter, idx] <- NA # delete NoGo RTs

}
cat("Finished simulations\n")
end.time <- Sys.time()
duration <- difftime(end.time, start.time)
duration

# ======================================================================================================== #
#### SAVE SIMULATED DATA: ####

cat("Set output file name\n")
simFile <- paste0(dirs$sims, "M", iMod_str, suffix, "_", simType, "_simulations.Rdata")
cat(paste0("Start saving ", simFile, "\n"))

## Save as R workspace:
save(data, selIter, simQGo, simQNoGo, simResp, simRT,
     file = simFile) # save as workspace

cat("END OF SCRIPT\n")

# END OF FILE.
