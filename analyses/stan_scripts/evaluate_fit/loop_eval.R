# ============================================================================================================= #
# loop_eval.R
# Simply loop over models to extract a certain metric.

cat("\014")
rm(list = ls())

# ============================================================================================================= #
#### Load packages: ####

require(stringr)
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
#### Loop over models: ####

WAIC <- c()
LOOIC <- c()

nMod <- 12
# suffix <- "_standard"
# suffix <- "_deltaIntSlope"
suffix <- "_sepParam"

# for (model2fit in 1:nMod){
for (model2fit in 1:nMod){
    
  # ---------------------------------------------------------------------- #
  ## Load data:
  
  cat("* -------------------------------------------------------- * \n")
  
  model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
  stanfitsfile <- paste0(dirs$results, "M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2
  load(stanfitsfile)
  
  # ---------------------------------------------------------------------- #
  ## Duration:
  
  cat(paste0(suffix, ": M", model2fit_str, ": ", capture.output(duration), "\n"))
  
  # ---------------------------------------------------------------------- #
  ## WAIC:
  
  cat(paste0(suffix, ": M", model2fit_str, ": WAIC = ", round(loo::waic(extract_log_lik(f2))$waic,1), "\n"))
  WAIC[model2fit] <- round(loo::waic(extract_log_lik(f2))$waic,1)
  
  # ---------------------------------------------------------------------- #
  ## LOOIC:
  
  cat(paste0(suffix, ": M", model2fit_str, ": LOOIC = ", round(loo(f2)$looic,1), "\n"))
  LOOIC[model2fit] <- round(loo(f2)$looic,1)
  
  # ---------------------------------------------------------------------- #
  ## Medium log-likelihood:
  
  posterior <- rstan::extract(f2, inc_warmup = FALSE, permuted = TRUE) # collapse across chains
  log_lik <- posterior$log_lik; median(log_lik)
  cat(paste0(suffix, ": M", model2fit_str, ": Median logLik = ", round(median(log_lik), 1), "\n"))

  # ---------------------------------------------------------------------- #
  ## Divergent transitions:
  # print(check_hmc_diagnostics(f2))

  sampler_params <- get_sampler_params(f2, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']; n = sum(divergent); N = length(divergent)
  cat(paste0(suffix, ": M", model2fit_str, ": Divergent transitions: ",
             str_pad(n, width = 3, side = "left", pad = "0"), " out of ", N, "\n"))

  # ---------------------------------------------------------------------- #
  ## Maximal tree depth reached:
  
  sampler_params <- get_sampler_params(f2, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__'];
  max_depth = 10; n = length(treedepths[sapply(treedepths, function(x) x == max_depth)]); N = length(treedepths)
  cat(paste0(suffix, ": M", model2fit_str, ": Maximal tree  depth: ",
             str_pad(n, width = 3, side = "left", pad = "0"), " out of ", N, "\n"))

  # ---------------------------------------------------------------------- #
  ## Rhat:
  
  parAvgTable <- summary(f2)$summary
  Rhat <- parAvgTable[, 10]
  cat(paste0(suffix, ": M", model2fit_str, ": max. Rhat = ", round(max(Rhat), 3), " by ", names(sort(Rhat, decreasing = T)[1]), "\n"))

  # ---------------------------------------------------------------------- #
  ## nEff:
  
  parAvgTable <- summary(f2)$summary
  nEff <- parAvgTable[, 9]
  cat(paste0(suffix, ": M", model2fit_str, ": min. nEff = ", round(min(nEff), 3), " by ", names(sort(nEff, decreasing = F)[1]), "\n"))

}

# ============================================================================================================= #
#### Load single data set: ####

## Settings:
model2fit <- 12
# suffix <- ""
suffix <- "_sepParam"

## Load and extract data:
model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
stanfitsfile <- paste0(dirs$results, "M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2
load(stanfitsfile) # load data
fit <- f2 # reassign
posterior <- rstan::extract(fit, inc_warmup = FALSE, permuted = TRUE) # extract all samples of each parameter
parAvgTable <- summary(fit)$summary # summary averaged over samples
paramVec <- rownames(parAvgTable) # extract names of all parameters
paramVec[grep("transf", paramVec)] 

# round(summary(fit)$summary[grep("transf", paramVec), ], 3)
round(summary(fit)$summary[grep("transf", paramVec), c(1, 5, 7)], 3)

# ============================================================================================================= #
#### Trace plots: ####

parName <- "transf_mu_"
# parName <- "sd_"

selVec <- paramVec[grep(parName,paramVec, fixed = T)] # sub-selection based on pattern
for (iPar in 1:length(selVec)){ # loop
  cat(paste0("Trace plot for parameter ",selVec[iPar], "\n"))
  print(traceplot(fit,pars = selVec[iPar])) # orange and purple
  # print(mcmc_trace(fit,  pars = paramVec[iPar])) # blue
  readline(prompt="Press [enter] to continue")
  if(iPar == length(selVec)){cat("Finished :-)\n")}
}

# ============================================================================================================= #
#### Posterior densities: ####

## _full:
parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaInt[")
parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaInt[", "deltaSlope[", "epsilon[")
parNameVec <- c("transf", "sd", "alpha[", "tau[", "betaWin[", "betaAvoid[", "deltaInt[", "deltaSlope[", "epsilon[")
parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "deltaSlope[", "deltaWin[", "deltaAvoid[", "epsilon[")

## _noInt:
parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta")
parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta", "epsilon[")
parNameVec <- c("transf", "sd", "alpha[", "tau[", "betaWin[", "betaAvoid[", "delta", "epsilon[")
parNameVec <- c("transf", "sd", "alpha[", "tau[", "beta[", "delta", "deltaWin[", "deltaAvoid[", "epsilon[")

# --------------------- #
## Loop and save:
WIDTH = 960; HEIGHT = 640
# WIDTH = 480; HEIGHT = 320
for (iPar in 1:length(parNameVec)){
  parName <- parNameVec[iPar]
  cat(paste0("Start plotting ",parName, "\n"))
  selVec <- paramVec[grep(parName,paramVec, fixed = T)] # sub-selection based on pattern
  png(paste0(plotdir,dataName, "_density_mod",model2fit_str,suffix, "_",parName, ".png"),width = WIDTH, height = HEIGHT, type ="cairo")
  print(stan_dens(fit, pars = selVec))
  # stan_hist(f2, pars = "X[1]") # Histogram - rather use density
  dev.off()
}

# --------------------- #
## Loop and inspect:
WIDTH = 960; HEIGHT = 640
# WIDTH = 480; HEIGHT = 320
parNameVec <- c("transf", "sd")
for (iPar in 1:length(parNameVec)){
  parName <- parNameVec[iPar]
  cat(paste0("Start plotting ",parName, "\n"))
  selVec <- paramVec[grep(parName,paramVec, fixed = T)] # sub-selection based on pattern
  # selVec <- selVec[(length(selVec)/2+1):length(selVec)]
  print(stan_dens(fit, pars = selVec))
  # stan_hist(f2, pars = "X[1]") # Histogram - rather use density
  readline(prompt="Press [enter] to continue")
}

# ============================================================================================================= #
#### Trade-offs: Bi-variate densities: ####

# ------------------------------------------------------ #
## Loop over parameters:

selVec <- paramVec[grep("transf",paramVec)] # sub-selection based on pattern
# selVec <- paramVec[grep("sd",paramVec)] # sub-selection based on pattern
# ------------------------------ #
## Load and save:
for (iPar1 in 1:length(selVec)){ # loop
  for (iPar2 in min(length(selVec),iPar1+1):length(selVec)){ # loop
    if (iPar1 != iPar2 ){
      cat("Create new plot...\n")
      cat(paste0("Bivariate density for parameters ",selVec[iPar1], " and ",selVec[iPar2], ":\n"))
      # Plot:
      png(paste0(plotdir,dataName, "_bivariate_density_mod",model2fit_str,suffix, "_",selVec[iPar1], "_",selVec[iPar2], ".png"),width = WIDTH, height = HEIGHT, type ="cairo")
      print(stan_scat(fit, pars = c(selVec[iPar1],selVec[iPar2])))
      # readline(prompt="Press [enter] to continue:")
      dev.off()
      if(iPar1 == length(selVec) & iPar2 == (length(selVec)-1)){cat("Finished :-)\n"); break}
    }
  }
}

# ------------------------------ #
## Load and inspect:
for (iPar1 in 1:length(selVec)){ # loop
  for (iPar2 in min(length(selVec),iPar1+1):length(selVec)){ # loop
    if (iPar1 != iPar2 ){
      cat("Create new plot...\n")
      cat(paste0("Bivariate density for parameters ",selVec[iPar1], " and ",selVec[iPar2], ":\n"))
      # Plot:
      # png(paste0(plotdir,dataName, "_bivariate_density_mod",model2fit_str,suffix, "_",selVec[iPar1], "_",selVec[iPar2], ".png"),width = WIDTH, height = HEIGHT, type ="cairo")
      print(stan_scat(fit, pars = c(selVec[iPar1],selVec[iPar2])))
      readline(prompt="Press [enter] to continue:")
      # dev.off()
      if(iPar1 == length(selVec) & iPar2 == (length(selVec)-1)){cat("Finished :-)\n"); break}
    }
  }
}

##  Single pairs:
# stan_scat(fit, pars = c("transf_mu_alpha", "transf_mu_beta")) # multivariate distribution---should be circle-shaped (independent normal distributions)

## Correlations:
# Cast posterior into data frame:
posterior <- as.data.frame(posterior)
# All correlations of group-level means in matrix:
library(ltm)
selVec <- paramVec[grep("transf",paramVec)] # sub-selection based on pattern
rcor.test(posterior[,selVec])
tmp <- as.matrix(unlist(rcor.test(posterior[,selVec])))
tmp[(1:length(tmp)/2)]
# Single correlations:
cor(posterior$transf_mu_alpha,posterior$transf_mu_tau)
cor(posterior$transf_mu_alpha,posterior$transf_mu_beta)
cor(posterior$transf_mu_alpha,posterior$transf_mu_deltaInt)
cor(posterior$transf_mu_alpha,posterior$transf_mu_deltaSlope)
cor(posterior$transf_mu_alpha,posterior$transf_mu_epsilon)
cor(posterior$transf_mu_beta,posterior$transf_mu_deltaSlope)
# ============================================================================================================= #
#### Auto-correlation: ####

#  --------------------- #
parNameVec <- c("transf", "sd")
## Loop and save:
WIDTH = 960; HEIGHT = 640
# WIDTH = 480; HEIGHT = 320
for (iPar in 1:length(parNameVec)){
  parName <- parNameVec[iPar]
  cat(paste0("Start plotting ",parName, "\n"))
  selVec <- paramVec[grep(parName,paramVec, fixed = T)] # sub-selection based on pattern
  # png(paste0(plotdir,dataName, "_autocorrelation_mod",model2fit_str,suffix, "_",parName, ".png"),width = WIDTH, height = HEIGHT, type ="cairo")
  print(stan_ac(fit, pars = selVec)) # average auto-correlation per lag
  readline(prompt="Press [enter] to continue:")
  # dev.off()
}

## Single plot:
# parName <- "transf"
# selVec <- paramVec[grep(parName,paramVec, fixed = T)] # sub-selection based on pattern
# stan_ac(fit, pars = selVec) # average auto-correlation per lag

# END