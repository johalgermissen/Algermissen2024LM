#!/usr/bin/env Rscript
# ======================================================================================================== #
# test_parameter_recovery.R
# Test simulating data and fitting model to it using Stan.
# Johannes Algermissen, 2023.
# ======================================================================================================== #

### Clear workspace and console.
cat("\014")
rm(list = ls()); 

# ======================================================================================================== #
#### 00) Load libraries: ####

require(stringr)

require(rstan)
require(bayesplot)
require(loo)

# ======================================================================================================== #
#### 00) Set directories: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)

## Load helper files:
suffix <- "_sepParam"
source(paste0(dirs$simscripts, "RLDDM_modSims", suffix, "_flat.R")) # model code for simulations
source(paste0(dirs$models, "RLDDM4stan", suffix, "_flat.R")) # model code for simulations
source(paste0(dirs$helpers, "initVals4stan", suffix, "_flat.R")) # model code for simulations

# ======================================================================================================== #
#### 01) General set up: ####

# ------------------------------------------------------------------- #
### Input and output parameters:

parType = c("alpha", "tau", "beta", "betaWin", "betaAvoid", "deltaInt", "deltaSlope", "deltaWin", "deltaAvoid", "epsilon", "pi", "theta") 
paramOut4Mod <- list(c(1,2,3,6), c(1,2,3,6,7,10), c(1,2,4,5,6,7,10), c(1,2,3,7,8,9,10),
                     c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11),
                     c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), 
                     c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12)
)

# ------------------------------------------------------------------- #
### Load data:

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
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
#### 02) Select model ID: ####

## Select model:
model2gen <- 1

# suffix <- "_sepParam"

nIter <- 1
iIter <- 1

# ======================================================================================================== #
#### 03) Generate data: ####

## Turn into model string:
model2gen_str <- str_pad(model2gen, width = 2, side = "left", pad = "0")

## Load previous model fit:
load(paste0(dirs$results, "fit_2023_08_15/M", model2gen_str, suffix, "_fitted.Rdata"))
nSample <- dim(f2)[1]; nChain <- dim(f2)[2]; nParamAll <- dim(f2)[3]

## Extract posterior:
posteriorGen <- rstan::extract(f2, inc_warmup = FALSE, permuted = TRUE)
meanIdx <- which(grepl("mu_", names(posteriorGen)) & !grepl("transf_", names(posteriorGen)))
sdIdx <- which(grepl("sd_", names(posteriorGen)))
nParam <- length(meanIdx)
cat(paste0("Summary for model ", model2gen_str, ":\n"))

# Inspect:
summary(f2)$summary[1:nParam, ]
nRound <- 2
round(summary(f2)$summary[meanIdx, 1], nRound) # mean
round(summary(f2)$summary[sdIdx, 1], nRound) # sd
cbind(round(summary(f2)$summary[meanIdx, 1], nRound), round(summary(f2)$summary[sdIdx, 1], nRound) )

cat("Adjust priors....\n")



source(paste0(dirs$models, "RLDDM4stan", suffix, "_flat.R")) # model code for simulations

## Identify subject-level parameter names:
allParNames <- names(posteriorGen) # names of all parameters
subParamNames <- allParNames[grepl("sub_", allParNames, fixed = T)] # names of subject level parameters
nParamSub <- length(subParamNames) # number subject level parameters

## Cast in one big array of dimensions subject x parameters x samples:
subParSampleMat <- array(NA, c(nSub, nParamSub, nSample*nChain)) # subjects, parameters, iterations
for (iParam in 1:nParamSub){ # iParam <- 1 # for each parameter
  subParSampleMat[, iParam, ] <- t(posteriorGen[[subParamNames[iParam]]]) # extract lists of all subjects for this parameter, transpose, save in matrix
}

## Mean per subject across samples:
subParMat <- apply(subParSampleMat, c(1, 2), mean) # average over samples

## Means and covariance across subjects:
subParamNames
meanVec <- apply(subParMat, c(2), mean); meanVec
sdVec <- apply(subParMat, c(2), sd); sdVec
covMat <- cov(subParMat); covMat

## Sample new parameters from multivariate normal:
# # nIter <- 100 # number samples to draw
set.seed(70) # set RNG seed
sampledIterParMat <- rmvnorm(n = nIter, mean = meanVec, sigma = covMat)
sampledIterParMat

## Data to simulate:
run_simulations <- eval(parse(text = paste0("modSim_M", model2gen_str)))

# ======================================================================================================== #
#### 04) Model to fit: ####

## Model index as string:
model2fit <- model2gen
model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")

## Model code:
source(paste0(dirs$models, "RLDDM4stan", suffix, "_flat.R")) # model code for simulations
mstr <- get(paste("mstr", model2fit_str, sep = "")) # string to be used as input for stan fitting

## Parameter initialization:
# source(paste0(dirs$helpers, "initVals4stan", suffix, "_flat.R")) # model code for simulations
init_fun <- eval(parse(text = paste0("init_fun_mod", model2fit_str)))

## Select parameters for selected model:
paramIdx <- paramOut4Mod[[model2fit]] # extract indices for selected model
paramInit <- parType[paramIdx]
paramOut <- c(paramInit, "log_lik") # add log likelihood
paramInit <- paste0("nd_", paramInit) # add nd_

### Select parameters:
# iIter <- 1
par <- sampledIterParMat[iIter, ]

# --------------------------------- #
### Simulate data:
# source(paste0(dirs$simscripts, "RLDDM_modSims", suffix, "_flat.R")) # model code for simulations
# run_simulations <- eval(parse(text = paste0("modSim_M", model2gen_str)))
dListSub <- run_simulations(data, par) # takes around 16 sec. for 54 people

# --------------------------------- #
### Fit, track duration:
nStanIter <- 200
nWarmup <- floor(nStanIter/2)
nChain <- 4
seed = 70 # arbitrary, but setting a fixed seed makes models reproducible
fit <- stan(model_code = mstr, data = dListSub, init = init_fun,
            warmup = nWarmup, iter = nStanIter, chains = nChain, cores = nChain,
            seed = seed+1, pars = paramOut)

# --------------------------------- #
### Evaluate output:
posterior <- rstan::extract(fit, inc_warmup = FALSE, permuted = TRUE) # collapse across chains
cat(paste0("Summary for model ", model2fit_str, ":\n"))
summary(fit)$summary[1:nParam, ] # extract summary statistics
round(par, 2)
round(summary(fit)$summary[1:nParam, 1], 2)

## END OF FILE.
