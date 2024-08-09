#!/usr/bin/env Rscript
# ======================================================================================================== #
# eval_parameter_recovery.R
# Evaluate parameter recovery.
# Johannes Algermissen, 2023.
# ======================================================================================================== #
#### 00a) Clear: ####

### Clear workspace and console.
cat("\014")
rm(list = ls()); 


# ======================================================================================================== #
#### 00b) Load libraries: ####

library(boot) # for inv.logit
library(stringr) # for str_pad

library(rstan) # for reading posterior
library(loo) # for extract

library(lattice)
library(pastecs)

library(RColorBrewer)
library(MetBrewer)
library(corrplot)

# ======================================================================================================== #
#### 00c) Set stan directories: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)

# ==================================================================================== #
#### 00d) Set regression directories: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/regression/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs.R")) # Load packages and options settings
regDirs <- set_dirs(rootDir)

## Load custom functions:
source(paste0(regDirs$codeDir, "functions/00_mgngstakes_functions_regression.R")) # Load functions

# ======================================================================================================== #
#### 00e) Possible parameters: ####

parType = c("alpha", "tau", "beta", "betaWin", "betaAvoid", "deltaInt", "deltaSlope", "deltaWin", "deltaAvoid", "epsilon", "pi", "theta") 
paramOut4Mod <- list(c(1,2,3,6), c(1,2,3,6,7,10), c(1,2,4,5,6,7,10), c(1,2,3,7,8,9,10),
                     c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11),
                     c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), 
                     c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12)
)


# ======================================================================================================== #
#### 01) Select model: ####

## Select model:
model2gen <- 12

# model2gen <- 3
model2gen_str <- str_pad(model2gen, width = 2, side = "left", pad = "0")

## Suffix:
suffix <- "_sepParam"

## Check inputs:
cat(paste0("Model to fit is ", model2gen_str), "\n")
cat(paste0("Suffix is ", suffix, "\n"))

# ------------------------------------------------------------------- #
### Input and output parameters:

## Select parameters for selected model:
paramIdx <- paramOut4Mod[[model2gen]] # extract indices for selected model
paramInit <- parType[paramIdx]
paramOut <- c(paramInit, "log_lik") # add log likelihood
paramInit <- paste0("nd_", paramInit) # add nd_

# ======================================================================================================== #
#### 02) Load fitted parameters: ####

## Directory where to load from:
# paramRecovDir <- dirs$paramRecov
# paramRecovDir <- paste0(dirs$paramRecov, "parameterRecovery_2023_12_08/")
# paramRecovDir <- paste0(dirs$paramRecov, "subParMvtnorm500/")
paramRecovDir <- paste0(dirs$paramRecov, "subParMvtnorm1000/")

## File name to load:
paramRecovFile <- paste0("M", model2gen_str, suffix, "_parameterRecovery.Rdata")
cat(paste0("Start loading ", paramRecovFile, "\n"))

# --------------------------------------- #
## Save as R workspace:
load(paste0(paramRecovDir, paramRecovFile))
cat("Finished loading :-)\n")

## Data dimensions:
nIter <- sum(complete.cases(fittedIterParMat[, 1]))
nParamSub <- ncol(fittedIterParMat)
cat(paste0("Found ", nIter, " iterations for ", nParamSub, " parameters\n"))

# ======================================================================================================== #
# ======================================================================================================== #
#### 03a) Corrplot: ####

## Combine:
corData <- data.frame(cbind(sampledIterParMat, fittedIterParMat))

## Assign names:
# https://stackoverflow.com/questions/35965433/r-corrplot-change-data-labels
paramNames <- paramOut[1:(length(paramOut) - 1)]
nParam <- length(paramNames)
input <- paramNames
for (iParam in 1:length(input)){ # iParam <- 1
  if (input[iParam] %in% c("alpha", "tau", "beta", "epsilon")){input[iParam] <- paste0("$", input[iParam])}
  if (grepl("beta", input[iParam]) && model2gen %in% c(3)){input[iParam] <- paste0("$beta[", gsub("beta", "", input[iParam]), "]")}
  if (grepl("delta", input[iParam])){input[iParam] <- paste0("$delta[", gsub("delta", "", input[iParam]), "]")}
  if (input[iParam] %in% c("$alpha") && model2gen %in% c(5, 9, 10)){input[iParam]  <- "$alpha[Low]"}
  if (input[iParam] %in% c("pi") && model2gen %in% c(5, 9, 10)){input[iParam]  <- "$alpha[High]"}
  if (input[iParam] %in% c("$tau") && model2gen %in% c(6, 9, 11, 12)){input[iParam]  <- "$tau[Low]"}
  if (input[iParam] %in% c("pi") && model2gen %in% c(6, 11)){input[iParam]  <- "$tau[High]"}
  if (input[iParam] %in% c("theta") && model2gen %in% c(9)){input[iParam]  <- "$tau[High]"}
  if (input[iParam] %in% c("pi") && model2gen %in% c(12)){input[iParam]  <- "$tau[High/Cong]"}
  if (input[iParam] %in% c("theta") && model2gen %in% c(12)){input[iParam]  <- "$tau[High/Incong]"}
  if (input[iParam] %in% c("$beta") && model2gen %in% c(7)){input[iParam]  <- "$beta[Low]"}
  if (input[iParam] %in% c("pi") && model2gen %in% c(7)){input[iParam]  <- "$beta[High]"}
  if (input[iParam] %in% c("pi") && model2gen %in% c(8)){input[iParam]  <- "$delta[High]"}
  if (input[iParam] %in% c("theta") && model2gen %in% c(10, 11)){input[iParam]  <- "$delta[High]"}
}
paramNames <- input
cat(paste0("Use parameter names ", paste0(paramNames, collapse = ", "), "\n"))
names(corData) <- c(paramNames, paramNames)

## Correlation matrix on complete cases:
# M <- cor(corData, use = "pairwise.complete.obs")
M <- cor(corData[, 1:nParam], corData[, (nParam + 1):ncol(corData)], use = "pairwise.complete.obs") # coefs in rows, questionnaires in columns: relaxed, pairwise complete cases
# M

# ------------------------------------------------------------------------------ #
### Plot:
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# https://cran.r-hub.io/web/packages/corrplot/vignettes/corrplot-intro.html
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf

## Standard plot (huge font size):
corrplot::corrplot(M, type = "upper",
                   addCoef.col = 'grey50', # darkgrey, grey50, grey60, orange
                   number.cex = 1.2, tl.cex = 2,
                   col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")

## Save:
png(paste0(dirs$plots, "parameter_recovery_corrplot_M", model2gen_str, "_nIter", nIter, ".png"), width = 480, height = 480)
corrplot::corrplot(M, type = "upper",
                   addCoef.col = 'grey50', # darkgrey, grey50, grey60, orange
                   number.cex = 1.2, tl.cex = 2,
                   col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
dev.off()

# ------------------------------------------------------------------------------ #
### Alternative versions:

# corrplot(M, method = 'ellipse', type = 'upper', col = rev(COL2('RdBu')), tl.col = "black")
# corrplot.mixed(M, 
#                lower = "number", upper = "circle", diag = "u",
#                lower.col = "black", upper.col = rev(COL2('RdBu')), tl.col = "black")

## Change font size manually, still print correlations into tiles:
CEX <- 0.6
corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt",
                   number.cex = CEX, tl.cex = CEX)

## Color only:
corrplot::corrplot(M, col = rev(COL2('RdBu')), tl.col = "black")

# ======================================================================================================== #
#### 03b) Plot scatterplot for each parameter: ####

## Inspect:
# sampledIterParMat[1:nIter, ]
# fittedIterParMat[1:nIter, ]

## Load custom functions:
source(paste0(regDirs$codeDir, "00_mgngstakes_functions_regression.R")) # Load functions

## Normal scatter plot:
for (iParam in 1:nParamSub){
  
  ## Extract parameter:
  x <- sampledIterParMat[1:nIter, iParam]
  y <- fittedIterParMat[1:nIter, iParam]
  xMin <- min(c(x, y))
  xMax <- max(c(x, y))
  xLim <- c(xMin, xMax)
  cat(paste0("Start parameters no. ", iParam, ": ", paramOut[iParam], "\n"))
  
  # ------------------------------------------------------------------ #
  ### Standard plot:
  
  ## Plot:
  # plot(x, y,
  #      xlim = c(xMin, xMax), ylim = c(xMin, xMax),
  #      xlab = "True parameter", ylab = "Fitted parameter",
  #      main = paste0(paramOut[iParam], ": r = ", cor(x, y)))
  # lines(c(xMin, xMax), c(xMin, xMax), lty = 2) # identity line
  # abline(lm(y ~ x), lwd = 3, col = "red") # regression line

  # ------------------------------------------------------------------ #
  ### Pretty plot:
  
  ## Print to console:
  paramLab <- rename_parameter(paramOut[iParam], model2gen, suffix) # rename
  paramLab <- string2Greek(paramLab) # prepare for LaTeX printing
  
  corData <- data.frame(x  = x, y = y)
  p <- plot_correlation(data = corData, xVar = "x", yVar = "y", 
                        # xLab = paste0("true $", paramNamesLaTeX[iParam], "$"), yLab = paste0("$fitted ", paramNamesLaTeX[iParam], "$"),
                        xLab = paste0("true $", paramLab, "$"), yLab = paste0("fitted $", paramLab, "$"),
                        # xLab = paste0(paramNamesLaTeX[iParam], "(true)"), yLab = paste0(paramNamesLaTeX[iParam], "(fitted)"),
                        # xLab = paramNamesLaTeX[iParam], yLab = paramNamesLaTeX[iParam], 
                        xLim = xLim, yLim = xLim, FTS = 28, useLaTeX = T)
  
  ## Save:
  plotName <- paste0("parameter_recovery_M", model2gen_str, "_", paramOut[iParam])
  cat(paste0("Save as ", plotName, "\n"))
  png(paste0(dirs$plots, plotName, ".png"), width = 480, height = 480)
  print(p)
  dev.off()
  
  # ------------------------------------------------------------------ #
  ### Print to console:
  cat("True parameters:\n")
  print(t(stat.desc(x)))
  cat("Fitted parameters:\n")
  print(t(stat.desc(y)))
  
  readline(prompt = "Press [enter] to continue")
  if(iParam == nParamSub){cat("Finished all parameters! :-)\n")}
}

# ======================================================================================================== #
#### 04a) Inspect diagonal of matrix: ####

# M
diagVec <- diag(M)
round(t(stat.desc(diagVec)), 2)

nIterPerm <- nIter
maxPermVec <- rep(NA, nIterPerm)
for (iIterPerm in 1:nIterPerm){
  permIdx <- sample(1:nIter, nIter, replace = F) # permute rows
  # corData <- data.frame(cbind(sampledIterParMat, fittedIterParMat))
  permData <- cbind(sampledIterParMat, fittedIterParMat[permIdx, ]) # permute only rows of fitted parameters
  M <- cor(permData[, 1:nParam], permData[, (nParam + 1):ncol(corData)], use = "pairwise.complete.obs") # coefs in rows, questionnaires in columns: relaxed, pairwise complete cases
  diagVec <- diag(M) # take only diagonal
  maxPermVec[iIterPerm] <- max(diagVec, na.rm = T)
}

## Inspect:
max(maxPermVec)
round(t(stat.desc(c(maxPermVec))), 3)
#      nbr.val nbr.null nbr.na    min   max range    sum median  mean SE.mean CI.mean.0.95 var std.dev coef.var
# [1,]    1000        0      0 -0.008 0.152 0.161 44.704  0.043 0.045   0.001        0.001   0    0.02    0.448
densityplot(c(maxPermVec))
round(quantile(c(maxPermVec), seq(0, 1, 0.05), na.rm = T), 3)
#     0%     5%    10%    15%    20%    25%    30%    35%    40%    45%    50%    55%    60%    65%    70%    75%    80%    85%    90%    95%   100% 
# -0.008  0.014  0.021  0.025  0.028  0.031  0.034  0.036  0.039  0.041  0.043  0.046  0.048  0.050  0.053  0.057  0.060  0.065  0.071  0.079  0.152 


# ======================================================================================================== #
#### 04b) Inspect convergence criteria: ####

### Rhat:
rhatMat
sum(c(rhatMat) > 1.1)
mean(c(rhatMat) > 1.1)
which(rhatMat > 1.1)

sum(c(rhatMat) > 1.01)
mean(c(rhatMat) > 1.01)
which(rhatMat > 1.01)

sum(c(rhatMat) > 1.001)
mean(c(rhatMat) > 1.001)
which(rhatMat > 1.001)

max(rhatMat)
sort(rhatMat, decreasing = T)[1:10]

### nEff:
nStanIter <- 1000 # at least 0.1% of iterations
neffMat
sum(c(neffMat) < nStanIter * 0.01)
mean(c(neffMat) < nStanIter * 0.01)
which(neffMat < nStanIter * 0.01)
min(neffMat)
sort(neffMat, decreasing = F)[1:10]

### WAIC:
plot(waicVec)
densityplot(waicVec)

### LOOIC:
plot(looVec)
densityplot(looVec)

### Duration:
plot(durationVec) # in seconds
sum(durationVec > 60)
mean(durationVec > 60)

### Divergent transitions:
plot(divergentVec)
sum(divergentVec > 0)
mean(divergentVec > 0)

### Max tree depth reached:
plot(maxTreeVec)
sum(maxTreeVec > 0)
mean(maxTreeVec > 0)

### Energy:
energyMat
sum(c(energyMat) < 0.2)
which(energyMat < 0.2)
min(energyMat)
sort(energyMat, decreasing = F)[1:10]

# END OF FILE.
