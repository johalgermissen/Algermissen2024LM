#!/usr/bin/env Rscript
# ======================================================================================================== #
# simulate_priors.R
# Generate normal distributions, transform, plot.
# Johannes Algermissen, 2023.

rm(list = ls())

## Load libraries:

library(lattice)
library(boot)

plotdir <- "<insert root dir here>"
if(!dir.exists(plotdir)){dir.create(plotdir)}

# ======================================================================================================= #
#### Functions: ####

# Phi_approx = pnorm

log1p_exp <- function(x){
  return(log(1+exp(x)))
}

# ======================================================================================================= #
#### Plot input vs. output after transformation: ####

# Define input space:
x <- seq(-10, 10, 0.01)
x <- seq(-1, 1, 0.01)
x <- seq(0, 1, 0.01)

## Linear and exponential transformations:
plot(x, x) # linaer
plot(x, exp(x)) # exponential
plot(x, log1p_exp(x)) # 1pexp
plot(x, inv.logit(x)) # inverse logit/ softmax

# ======================================================================================================= #
#### Simulate normal distribution: ####

## Parameters:
nSample <- 100000
nMean <- 0
nSD <- 3

## Simulate:
tmp <- rnorm(nSample, nMean, nSD)

# ======================================================================================================= #
#### Plot standard distribution: ####

png(paste0(plotdir,"density_M", nMean,"_SD", nSD, "_untranformed.png"), width = 480, height = 480)
plot(density(tmp), lwd = 2,
     xlab = "X", main = paste0("Normal distribution, M = ", nMean,", SD = ", nSD,", untransformed"))
dev.off()

densityplot(tmp, lwd = 2,
            xlab="X", main = paste0("Normal distribution, M = ",nMean,", SD = ", nSD,", untransformed"))

# ======================================================================================================= #
#### log1p_exp: ####
# https://mc-stan.org/docs/2_21/functions-reference/composed-functions.html

plot(density(exp(tmp)))

tmp_transf <- log1p_exp(tmp)

png(paste0(plotdir,"density_M",nMean,"_SD",nSD,"_log1p_exp.png"), width = 480, height = 480)
plot(density(tmp_transf), lwd = 2,
     xlab = "X", main = paste0("Normal distribution, M = ",nMean,", SD = ", nSD,", transformed via log1p_exp"))
dev.off()

densityplot(tmp_transf)

# ======================================================================================================= #
#### Phi-approx: ####
# https://mc-stan.org/docs/2_21/functions-reference/Phi-function.html
# https://www.r-bloggers.com/2013/02/normal-distribution-functions/

tmp_transf <- pnorm(tmp)

png(paste0(plotdir,"density_M",nMean,"_SD",nSD,"_phiapprox.png"), width = 480, height = 480)
plot(density(tmp_transf), lwd = 2,
     xlab = "X", main = paste0("Normal distribution, M = ",nMean,", SD = ", nSD,", transformed via Phi_approx (pnorm)"))
dev.off()

densityplot(tmp_transf)

# ======================================================================================================= #
#### inv.logit: ####
# https://mc-stan.org/docs/2_21/functions-reference/link-functions.html

tmp_transf <- inv.logit(tmp)

png(paste0(plotdir,"density_M",nMean,"_SD",nSD,"_invlogit.png"), width = 480, height = 480)
plot(density(tmp_transf), lwd = 2,
     xlab = "X", main = paste0("Normal distribution, M = ",nMean,", SD = ", nSD,", transformed via inv.logit (softmax)"))
dev.off()

densityplot(tmp_transf)

## Control: uniform distribution
plot(density(runif(nSample, 0, 1)))
densityplot(runif(nSample, 0, 1))

# ======================================================================================================= #
# ======================================================================================================= #
# ======================================================================================================= #
# ======================================================================================================= #
#### Non-centered parameterization: ####

plotdir <- "C:/Users/Johannes Algermissen/Documents/AAANijmegen/3_PhD_Hanneke/ActiveSampleBias/RLDDM/RLDDMs_Vanessa/priors_noncentered/"
if(!dir.exists(plotdir)){dir.create(plotdir)}

nSample <- 10000

muMean <- 5
muSD <- 10
sigmaMean <- 0
sigmaSD <- 1

## Simulate:
mu <- rnorm(nSample,muMean,muSD)
sigma <- rnorm(nSample,sigmaMean,sigmaSD)
z <- rnorm(nSample,0,1)

## Check:
plot(density(mu))
plot(density(sigma))
plot(density(z))

# ------------------------------------------------- #
### exp:

param <- mu + z * sigma # exp
png(paste0(plotdir,"density_mu_M",muMean,"_SD",muSD,"_sigma_M",sigmaMean,"_SD",sigmaSD,"_noncentered_untransformed.png"), width = 480, height = 480)
plot(density(param), lwd = 2, 
     xlab = "X", main = paste0("Normal distribution, Mu: M = ",muMean,", SD = ", muSD,", Sigma: M = ",sigmaMean,", SD = ", sigmaSD,",\nnon-centered parameterization, untransformed"))
dev.off()

# ------------------------------------------------- #
### exp:

param <- exp(mu + z * sigma) # exp
png(paste0(plotdir,"density_mu_M",muMean,"_SD",muSD,"_sigma_M",sigmaMean,"_SD",sigmaSD,"_noncentered_exp.png"), width = 480, height = 480)
plot(density(param), lwd = 2, xlim = c(0,1000000000),
     xlab = "X", main = paste0("Normal distribution, Mu: M = ",muMean,", SD = ", muSD,", Sigma: M = ",sigmaMean,", SD = ", sigmaSD,",\nnon-centered parameterization, transformed via exp"))
dev.off()

# ------------------------------------------------- #
### log1p_exp:

param <- log1p_exp(mu + z * sigma) # log1p_exp
png(paste0(plotdir,"density_mu_M",muMean,"_SD",muSD,"_sigma_M",sigmaMean,"_SD",sigmaSD,"_noncentered_log1p_exp.png"), width = 480, height = 480)
plot(density(param), lwd = 2,
     xlab = "X", main = paste0("Normal distribution, Mu: M = ",muMean,", SD = ", muSD,", Sigma: M = ",sigmaMean,", SD = ", sigmaSD,",\nnon-centered parameterization, transformed via log1p_exp"))
dev.off()

# ------------------------------------------------- #
### Phi_approx:

param <- pnorm(mu + z * sigma) # exp
png(paste0(plotdir,"density_mu_M",muMean,"_SD",muSD,"_sigma_M",sigmaMean,"_SD",sigmaSD,"_noncentered_phiapprox.png"), width = 480, height = 480)
plot(density(param), lwd = 2,
     xlab = "X", main = paste0("Normal distribution, Mu: M = ",muMean,", SD = ", muSD,", Sigma: M = ",sigmaMean,", SD = ", sigmaSD,",\nnon-centered parameterization, transformed via Phi_approx (pnorm)"))
dev.off()

# plot(density(inv.logit(rnorm(nSample,0,1))))

# ------------------------------------------------- #
### inv.logit:

param <- inv.logit(mu + z * sigma) # inv.logit
png(paste0(plotdir,"density_mu_M",muMean,"_SD",muSD,"_sigma_M",sigmaMean,"_SD",sigmaSD,"_noncentered_invlogit.png"), width = 480, height = 480)
plot(density(param), lwd = 2,
     xlab = "X", main = paste0("Normal distribution, Mu: M = ",muMean,", SD = ", muSD,", Sigma: M = ",sigmaMean,", SD = ", sigmaSD,",\nnon-centered parameterization, transformed via inv.logit (softmax)"))
dev.off()

# plot(density(inv.logit(rnorm(nSample,0,1))))

# END
