#!/usr/bin/env Rscript
# ======================================================================================================== #
# fit_stan_main.R
# Model fitting using rstan.
# Johannes Algermissen, 2021.
# ======================================================================================================== #
#### RETRIEVE INPUTS FROM BASH: ####

### Clear workspace and console.
cat("\014\n")
rm(list = ls()); 

### Retrieve inputs from bash:
args <- commandArgs(trailingOnly = TRUE)

### Extract info:
## Model to fit:
model2fit <- as.numeric(args[1])
if (is.na(model2fit)){
  model2fit <- 1
  cat("Errorneous input argument, set model2fit to", model2fit, "\n")
}

suffix <- args[2]
if (is.na(suffix)){
  suffix <- "_standard"
  cat("Not suffix provided, set to", suffix, "\n")
}

isTest <- FALSE

# ------------------------------------- #
# cat("USE HARDCODED DEFAULT SETTINGS\n")
# model2fit <- 12
# suffix <- "_sepParam" # _standard _deltaIntSlope _sepParam
# isTest <- TRUE

# model2fit <- 12
# suffix <- "_deltaIntSlope" # _standard _deltaIntSlope _sepParam
# isTest <- TRUE

# ======================================================================================================== #
## Convert to string:

require(stringr)
model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")

## suffix of model version to use:
# suffix <- as.character(args[2])

## Check inputs:
cat(paste0("Model to fit is ", model2fit_str), "\n")
cat(paste0("Suffix is ", suffix, "\n"))

# ======================================================================================================== #
#### LOAD LIBRARIES: ####

## Load libraries:
if (!require("rstan",character.only = TRUE)){
	require("rstan"); 
  cat("Load stan library\n")
}

## For multiple cores:
require(parallel);
options(mc.cores = parallel::detectCores()) # for multiple cores

# ======================================================================================================== #
#### SET DIRECTORIES: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)

obs4stanFile    = paste0(dirs$results, "data4stan", ".Rdata") # turn mat-file into vector for stan

## Where output file is saved:
if(!dir.exists(dirs$results)){dir.create(dirs$results)}
stanFitsFile    = paste0(dirs$results, "M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2
# if(exists(stanFitsFile)){rm(stanFitsFile)} # delete if already exists

# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
# ======================================================================================================== #
#### LOAD HELPER SCRIPTS: ####

source(paste0(dirs$helpers, "prepObs4stan.R")) # prepares singletrial matlab data for use in Rstan.
source(paste0(dirs$helpers, "initVals4stan", suffix, ".R")) # prepares singletrial matlab data for use in Rstan.
source(paste0(dirs$models, "RLDDM4stan", suffix, ".R")) # specification of model used for fitting.
source(paste0(dirs$helpers, "fitmodel.R")) # performs actual model fitting.

# set working directory.
# setwd(dirs$root)

# ======================================================================================================== #
#### PREPARE DATA: ####

if(!file.exists(obs4stanFile)){ # if fMRIEEGPav4stan.Rdata doesn't exist
  # prepObs4stan <- function(inputDir, targetFile, filePattern = ".csv", sub2exclude = c()) {
  # inputDir = dirs$raw; targetFile = obs4stanFile; filePattern = ".csv"; sub2exclude = c(52) # dry settings
  prepObs4stan(inputDir = dirs$raw, targetFile = obs4stanFile, filePattern = ".csv", sub2exclude = c(52)) # saves results
} else {
  cat("Input data already exists, do not generate again; delete file to change it\n")
}

# ---------------------------------------------------- #
## Inspect input data:

# load(obs4stanFile)
# str(dList)

# Check for NAs:
# for (iVar in 1:length(dList)){
#   cat(paste0("Variable ", names(dList)[iVar], ": ", sum(is.na(dList[iVar])), " NAs\n"))
# }

# Check single variables:
# names(dList)
# dList$nStim
# dList$nResp
# dList$nTrial
# dList$nSub
# dList$nData # dList$nSub * dList$nTrial
# table(dList$stimuli) # balanced
# table(dList$resp)
# library(lattice); densityplot(as.numeric(dList$rt)) # note many 0s for NoGo
# table(dList$outcome) # +1 or 0 or -1
# dList$Qi

# ======================================================================================================== #
#### FIT MODEL: ####

cat("Ready for fitting model\n")
if(!file.exists(stanFitsFile)){ # if fits_obsData.Rdata doesn?t exist
  cat("Fit model\n")
  fitmodel(model2fit, obs4stanFile = obs4stanFile,  
           stanFitsFile = stanFitsFile, 
           isTest = isTest)
} else {
  cat("Output file", stanFitsFile, "already exist; do not fit again; delete old output file to refit\n")
}

cat("END OF SCRIPT")
# ======================================================================================================== #
# END OF SCRIPT.