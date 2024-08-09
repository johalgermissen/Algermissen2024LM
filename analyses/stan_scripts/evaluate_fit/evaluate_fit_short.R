#### evaluate_fit_short.R ####

cat("\0140")
rm(list = ls())

# ==================================================================================== #
#### Load packages: ####

# Base packages:
require(lattice)
require(stringr)

# Rstan family:
require(rstan)
require(bayesplot)
require(loo)

options(scipen = 20)

# ==================================================================================== #
#### Initialize directories: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)

# ==================================================================================== #
#### Add utility files: ####

# waicfile <- paste0(scriptdir, "waic.R")
# source(paste0(scriptdir, "getWAIC.R"))
source(paste0(dirs$helpers, "stan_utility.R"))
# https://github.com/betanalpha/knitr_case_studies/blob/master/rstan_workflow/stan_utility.R

# Mode of vector:
Mode <- function(x) {
  ux <- unique(x) # all variables
  ux[which.max(tabulate(match(x, ux)))] # highest frequency count
}

# ================================================================================================================================= #
# ================================================================================================================================= #
# ================================================================================================================================= #
# ================================================================================================================================= #
#### SELECT MODEL: ####

# Delete in case of previously loaded model:
if (exists("f1")){rm(f1, f2, dList, posterior, parSummaryTable)}

model2fit <- 12
# suffix <- "_standard" # _fixThresh
# suffix <- "_deltaIntSlope" # _fixThresh
suffix <- "_sepParam" # _fixThresh
model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")

# ==================================================================================== #
#### Load output file: ####

## Define output file name:
stanfitsfile = paste0(dirs$results, "M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2

## Load output file:
cat("Start loading model", model2fit_str, "\n")
load(stanfitsfile)
cat("Finished loading model", model2fit_str, "\n")

## Reassign name of object:
fit <- f2 # 
# rm(f1,f2) # delete old fitting objects
nSample <- dim(fit)[1]; nChain <- dim(fit)[2]; nParam <- dim(fit)[3]
cat(paste0("Model has ", nSample, " samples per chain, ", nChain, " chains, ", nParam, " parameters\n")) # nSamples, nChains, nParameters

# ==================================================================================== #
#### Duration: ####

# duration
cat(paste0("Duration in hours: ", 
           round(as.numeric(duration, units = "hours"), 1), "\n"))

# ==================================================================================== #
#### Extract all samples of each parameter: ####
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html

posterior <- rstan::extract(fit, inc_warmup = FALSE, permuted = TRUE) # collapse across chains

# ==================================================================================== #
#### Names of all parameters: ####
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html#posterior-summary-statistics-and-convergence-diagnostics

## Extract parameter names for later:
parSummaryTable <- summary(fit)$summary
paramVec <- rownames(parSummaryTable) # extract names of all parameters (row names)
# dim(parSummaryTable) # nParameters x 10 summary statistics
# head(parSummaryTable)

# ==================================================================================== #
#### Compute WAIC: ####

# waicfile     = paste0(dirs$results, "M", model2fit_str, "_waic.Rdata")
# getWAIC(stanfitsfile,waicfile, model2fit_str)
## does not work because exp(-745) smallest number for which exp > 0
# https://stackoverflow.com/questions/21013309/calcute-the-expoenential-of-big-negative-value-in-r

# Try loo package:
loo::waic(extract_log_lik(fit))

# ==================================================================================== #
#### Compute LOO: ####
# https://avehtari.github.io/modelselection/CV-FAQ.html

loo(fit)
# help('pareto-k-diagnostic')

# kfold(fit)

# ==================================================================================== #
#### Assess log-likelihood: ####

log_lik <- posterior$log_lik; cat("Median log-likelihood:", median(log_lik), "\n")

# ==================================================================================== #
#### Evaluate traceplots: ####

parName <- "transf_mu_"

## Select all parameters with this character subset:
selParamVec <- paramVec[grep(parName, paramVec, fixed = T)] # sub-selection based on pattern

## Traceplot for given parameter name:
if (length(dev.list() != 0)){dev.off()}
for (iPar in 1:length(selParamVec)){ # loop
  cat(paste0("Trace plot for parameter ", selParamVec[iPar], "\n"))
  print(traceplot(fit,pars = selParamVec[iPar])) # orange and purple
  # print(mcmc_trace(fit,  pars = selParamVec[iPar])) # blue
  readline(prompt = "Press [enter] to continue")
  if(iPar == length(selParamVec)){cat("Finished :-)\n")}
}

# print(traceplot(fit,pars = "sub_beta[1]")) # orange and purple
# ==================================================================================== #
#### Evaluate divergent transitions, maximum tree depth, energy: ####
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html#sampler-diagnostics
# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html

## Check all together:
check_hmc_diagnostics(fit)

# -------------------------------------------------------- #
### 1) Extract divergent transitions manually:

# sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
# divergent <- do.call(rbind, sampler_params)[,'divergent__']
# n = sum(divergent); n

# -------------------------------------------------------- #
### 2) Extract maximum tree-depth reached manually:

# sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
# treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
# max_depth = 10
# n = length(treedepths[sapply(treedepths, function(x) x == max_depth)]); n

# -------------------------------------------------------- #
### 3) Extract energy Bayesian fraction of missing information (E-BFMI) manually:

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
for (n in 1:length(sampler_params)) { # loop over chains # n <- 1
  energies = sampler_params[n][[1]][,'energy__']
  numer = sum(diff(energies)**2) / length(energies)
  denom = var(energies)
  cat(paste0("Chain ", n, ": E-BFMI = ", numer/denom, "\n")) # shouldn't be < 0.2
}  

# ==================================================================================== #
#### Evaluate Rhats: ####

## Code by Michael Betancourt: check if any rhats > 1.1
check_rhat(fit) 

## Extract rhats from summary:
Rhat <- parSummaryTable[, 10]

# ------------------------------------------------------ #
## Maximum:
cat("Maximum Rhat =", max(Rhat), "\n")

## Sorted in descending order:
sort(Rhat, decreasing = T)[1:10]
# sort(Rhat, decreasing = T)[1:40]

# ------------------------------------------------------ #
## Above certain criterion:

crit <- 2
cat(paste0(sum(Rhat>crit), " parameters for which Rhat > ", crit, "\n"))
# which(Rhat>crit)

crit <- 1.1
cat(paste0(sum(Rhat>crit), " parameters for which Rhat > ", crit, "\n"))
# which(Rhat>crit)

crit <- 1.01
cat(paste0(sum(Rhat>crit), " parameters for which Rhat > ", crit, "\n"))
# which(Rhat>crit)

crit <- 1.001
cat(paste0(sum(Rhat > crit), " parameters for which Rhat > ", crit, "\n"))
which(Rhat > crit)

# ==================================================================================== #
#### Evaluate nEffs: ####

## Code by Michael Betancourt: check if any Neff ratio of < 0.001
check_n_eff(fit) # check whether 0.001 > of iterations

# Extract nEff from summary:
nEff <- parSummaryTable[, 9]

# ------------------------------------------------------ #
# Minimum:
cat("Minimum nEff =", min(nEff), "\n")

## Sorted in ascending order:
sort(nEff, decreasing = F)[1:10]
# sort(nEff, decreasing = F)[1:40]

# cor(Rhat, nEff)
# plot(Rhat, nEff)

# ------------------------------------------------------ #
# Below certain criterion:
# crit <- 3
# cat(paste0(sum(nEff<crit), " parameters for which nEff < ",crit, "\n"))
# which(nEff<crit)

crit <- 10
cat(paste0(sum(nEff<crit), " parameters for which nEff < ",crit, "\n"))
# which(nEff<crit)

crit <- 100
cat(paste0(sum(nEff<crit), " parameters for which nEff < ",crit, "\n"))
# which(nEff<crit)

crit <- 1000
cat(paste0(sum(nEff<crit), " parameters for which nEff < ",crit, "\n"))
which(nEff<crit)

# ==================================================================================== #
#### Plot posterior densities: ####

# stan_dens(fit, pars = "transf_mu_beta") # next to each other (horizontally)

# -------------------------------------------------- #
## Multiple densities in 1 panel: loop and save (all in one plot):

# List of all parameters based on model:
cat(paste0("Select all possible parameters for model ", model2fit, "\n"))
if (suffix == "_standard"){
  if (model2fit == 1){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[")
  } else if(model2fit == 2){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "epsilon[")
  } else if(model2fit %in% c(3)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "betaWin[", "betaAvoid[", "delta[", "epsilon[")
  } else if(model2fit %in% c(4)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "deltaWin[", "deltaAvoid[", "epsilon[")
  } else {
    stop("Unknown model2fit")
  }
} else if (suffix == "_deltaIntSlope"){
  if (model2fit == 1){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaInt[")
  } else if(model2fit == 2){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaInt[", "deltaSlope[", "epsilon[")
  } else if(model2fit %in% c(3)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "betaWin[", "betaAvoid[", "deltaInt[", "deltaSlope[", "epsilon[")
  } else if(model2fit %in% c(4)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[")
  } else if(model2fit %in% c(5, 6, 7, 8)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[", "pi[")
  } else if(model2fit %in% c(9, 10, 11)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[", "pi[", "theta[")
  } else if(model2fit %in% c(12)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[", "piCong[", "piIncong[")
  } else {
    stop("Unknown model2fit")
  }
} else if (suffix == "_sepParam"){
  if (model2fit == 1){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[")
  } else if(model2fit == 2){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "epsilon[")
  } else if(model2fit %in% c(3)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "betaWin[", "betaAvoid[", "delta[", "epsilon[")
  } else if(model2fit %in% c(4)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "deltaWin[", "deltaAvoid[", "epsilon[")
  } else if(model2fit %in% c(5:8)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "deltaWin[", "deltaAvoid[", "epsilon[", "pi[")
  } else if(model2fit %in% c(9:14)){
    parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta[", "deltaWin[", "deltaAvoid[", "epsilon[", "pi[", "theta[")
  } else {
    stop("Unknown model2fit")
  }
} else {
  stop("Unknown suffix")
}
cat(paste0("Parameters are ", paste(parNameVec,collapse = ", "), "\n"))

# ------------------------------------------------------ #
## Loop over parameter types (transf_mu, sd, subject ones, plot as panels):

if (length(dev.list() != 0)){dev.off()}
WIDTH = 960; HEIGHT = 640
# WIDTH = 480; HEIGHT = 320
for (iPar in 1:length(parNameVec)){
  parName <- parNameVec[iPar]
  cat(paste0("Start plotting ", parName, "\n"))
  selParamVec <- paramVec[grep(parName, paramVec, fixed = T)] # sub-selection based on pattern
  print(stan_dens(fit, pars = selParamVec))
  readline(prompt = "Press [enter] to continue")
  if(iPar == length(parNameVec)){cat("Finished :-)\n")}
}

## Densities of selected parameters:
# selParamVec <- paramVec[grep("transf", paramVec, fixed = T)] # sub-selection based on pattern
# print(stan_dens(fit, pars = selParamVec))
# print(stan_dens(fit, pars = "transf_mu_piCong"))
# print(stan_dens(fit, pars = "transf_mu_piIncong"))
print(stan_dens(fit, pars = "transf_mu_pi"))
print(stan_dens(fit, pars = "transf_mu_theta"))

# ==================================================================================== #
#### Inspect parameter values: ####

library(pastecs)

## Select data:
paramData <- data.frame(posterior) # cast into data frame
selData <- paramData[, grepl("transf", names(paramData))]
# selData <- paramData[, grepl("sd", names(paramData))]
# selData <- paramData[, grepl("sub_piCong", names(paramData))]
# selData <- paramData[, grepl("sub_piIncong", names(paramData))]
names(selData)

## Group means:
# densityplot(selData$transf_mu_piCong)
# densityplot(selData$transf_mu_piIncong)
# t(stat.desc(selData$transf_mu_piCong))
# t(stat.desc(selData$transf_mu_piIncong))
densityplot(selData$transf_mu_pi)
densityplot(selData$transf_mu_pi)
t(stat.desc(selData$transf_mu_theta))
t(stat.desc(selData$transf_mu_theta))

## Difference:
difVec <- selData$transf_mu_pi - selData$transf_mu_tau
difVec <- selData$transf_mu_theta - selData$transf_mu_pi
difVec <- selData$transf_mu_theta - selData$transf_mu_pi

densityplot(difVec)
round(t(stat.desc(difVec)), 3)
round(quantile(difVec, seq(0, 1, 0.05), na.rm = T), 3)
round(quantile(difVec, c(0.025, 0.50, 0.975), na.rm = T), 3)

## Groups SDs:
# densityplot(selData$sd_piCong)
# densityplot(selData$sd_piIncong)
densityplot(selData$sd_pi)
densityplot(selData$sd_theta)

stat.desc(selData$transf_mu_theta)
plot(selData$transf_mu_theta, selData$transf_mu_pi)

## Correlation sub means:
paramData <- data.frame(posterior)
# selData <- paramData[, grepl("sub_piCong", names(paramData))]
selData <- paramData[, grepl("sub_pi", names(paramData))]
a <- colMeans(selData)
plot(a)
# selData <- paramData[, grepl("sub_piIncong", names(paramData))]
selData <- paramData[, grepl("sub_theta", names(paramData))]
b <- colMeans(selData)
plot(b)

## Plot against each other:
plot(a, b)

## Pretty plot:
maxLim <- 0.4
plot(a, b, xlim = c(-1*maxLim, maxLim), ylim = c(-1*maxLim, maxLim),
     main = paste0("r = ", round(cor(a, b), 3)))
abline(lm(b ~ a), col = "red", lwd = 3)
abline(a = 0, b = 1, lty = 2)


# ==================================================================================== #
#### Plot bi-variate distributions (scatter plots): ####

## Select all parameters with this character subset:
parName <- "transf"
selParamVec <- paramVec[grep(parName, paramVec)] # sub-selection based on pattern

## Loop and inspect:
if (length(dev.list() != 0)){dev.off()}
for (iPar1 in 1:(length(selParamVec)-1)){ # loop
  for (iPar2 in min(length(selParamVec)-1, iPar1+1):length(selParamVec)){ # loop
    if (iPar1 != iPar2){
      
      cat("Create new plot ...\n")
      cat(paste0("Bivariate density for parameters ", selParamVec[iPar1], " and ", selParamVec[iPar2], ":\n"))
      
      ## Print correlation:
      corVal <- cor(as.numeric(unlist(posterior[selParamVec[iPar1]])), 
                    as.numeric(unlist(posterior[selParamVec[iPar2]])))
      cat(paste0("Correlation r = ", corVal, "\n"))
      if (abs(corVal) > 0.30){cat(paste0("ATTENTION: correlation of ", corVal, "!!!!!!!!!!!!\n"))}
      
      ## Create plot:
      # bivariate normal distribution---should be circle-shaped (independent normal distributions)
      print(stan_scat(fit, pars = c(selParamVec[iPar1], selParamVec[iPar2])))
      readline(prompt = "Press [enter] to continue:")
      if(iPar1 == (length(selParamVec) - 1) & iPar2 == length(selParamVec)){cat("Finished :-)\n"); break}
    }
  }
}
# cat("Finished!\n") # skips first plot

# print(stan_scat(fit, pars = c("transf_mu_beta", "transf_mu_deltaAvoid")))
# cor(posterior$transf_mu_beta, posterior$transf_mu_deltaAvoid)

# print(stan_scat(fit, pars = c("transf_mu_tau", "transf_mu_piCong")))
# cor(posterior$transf_mu_tau, posterior$transf_mu_piCong)

# ------------------------------------------------- #
### Plot correlation matrix:
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

## Select variables, turn into data frame:
paramData <- data.frame(posterior)
corData <- paramData[, grepl("transf", names(paramData))]
# corData <- paramData[, grepl("sd", names(paramData))]
M <- cor(corData) # compute correlations

library(corrplot)
# corrplot(M)
corrplot::corrplot(M, 
                   addCoef.col = "black", 
                   col = rev(COL2(diverging = "RdBu")),
                   tl.col = "black", tl.pos = "lt")


# ==================================================================================== #
#### Average autocorrelation per lag: ####

## Select all parameters with this character subset:
parName <- "transf"
# parName <- "beta"
selParamVec <- paramVec[grep(parName, paramVec, fixed = T)] # sub-selection based on pattern

## Plot for selected parameter:
stan_ac(fit, pars = selParamVec) # average auto-correlation per lag

## END OF FILE.