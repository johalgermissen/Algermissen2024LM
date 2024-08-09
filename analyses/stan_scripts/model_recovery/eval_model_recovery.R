#!/usr/bin/env Rscript
# ======================================================================================================== #
# eval_model_recovery.R
# Evaluate model recovery.
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
library(viridis)
library(corrplot)

# ======================================================================================================== #
#### 00c) Set directories: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
dirs <- set_dirs_stan(rootDir)
targetDir <- paste0(dirs$modelRecov, "modelRecov_nIter1000/"); nIter <- 1000

# ======================================================================================================== #
#### 00d) Possible parameters: ####

parType = c("alpha", "tau", "beta", "betaWin", "betaAvoid", "deltaInt", "deltaSlope", "deltaWin", "deltaAvoid", "epsilon", "pi", "theta") 
paramOut4Mod <- list(c(1,2,3,6), c(1,2,3,6,7,10), c(1,2,4,5,6,7,10), c(1,2,3,7,8,9,10),
                     c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11), c(1,2,3,7,8,9,10,11),
                     c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), 
                     c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12), c(1,2,3,7,8,9,10,11,12)
)

# ======================================================================================================== #
#### 01) Load and sort statistics into data frame: ####

### Data dimensions:
nMod <- 12
nMod2Sim <- nMod
nMod2Fit <- nMod
# nIter <- 1000
nRow <- nMod2Sim * nMod2Fit * nIter
suffix <- "_sepParam"

# ----------------------------------------------------- #
### Initialize data:

simData <- data.frame(
  model2gen = rep(1:nMod2Sim, each = nMod2Fit * nIter),
  model2fit = rep(1:nMod2Fit, each = nIter, times = nMod2Sim),
  iteration = rep(1:nIter, times = nMod2Sim * nMod2Fit)
) 
simData$WAIC <- NA
simData$LOOIC <- NA
simData$maxRhat <-NA
simData$minNeff <- NA
simData$divergent <- NA
simData$maxTree <- NA
simData$minEnergy <- NA

# ----------------------------------------------------- #
### Loop over models and iterations, read:

start2Gen <- 1; start2Fit <- 1
for (model2gen in start2Gen:nMod2Sim){ # model2gen <- 1
  for (model2fit in start2Fit:nMod2Fit){ # model2fit <- 1

    ## Model IDs as strings:
    model2gen_str <- str_pad(model2gen, width = 2, side = "left", pad = "0")
    model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
    
    ## Load file:
    fileName <- paste0("genM", model2gen_str, "_fitM", model2fit_str, suffix, "_modelRecovery_fit.Rdata")
    cat(paste0("Load ", fileName, ": "))
    fullFileName <- paste0(targetDir, fileName)
    load(fullFileName)
    
    ## Data dimensions:
    nIter <- sum(complete.cases(fittedIterParMat[, 1]))
    nParamSub <- ncol(fittedIterParMat)
    cat(paste0("Found ", nIter, " iterations for ", nParamSub, " parameters\n"))
    
    ## Locate where to store data:
    rowIdx <- which(simData$model2gen == model2gen & simData$model2fit == model2fit)
    
    ## Sort WAIC and LOOIC:
    simData$WAIC[rowIdx] <- waicVec
    simData$LOOIC[rowIdx] <- looVec
    
    ## Sort Rhat, nEff, divergent transitions:
    simData$maxRhat[rowIdx] <- apply(rhatMat, 1, max, na.rm = T)
    simData$minNeff[rowIdx] <- apply(neffMat, 1, min, na.rm = T)
    simData$divergent[rowIdx] <- divergentVec
    simData$maxTree[rowIdx] <- maxTreeVec
    simData$minEnergy[rowIdx] <- apply(energyMat, 1, min, na.rm = T)
    
  } # end model2fit
  if (model2gen == nMod2Sim){cat("Finished loading all fits! :-)\n")}
} # end model2gen

## Save:
write.csv(simData, paste0(targetDir, "summaryData_nIter,", nIter, ".csv"), row.names = F)

## Load:
# simData <- read.csv(paste0(targetDir, "summaryData.csv"))

# -------------------------------------------------------------- #
### Inspect:

## Missing simulations:
sum(is.na(simData$WAIC))
sum(is.na(simData$LOOIC))
simData[is.na(simData$LOOIC), c("model2gen", "model2fit")] # where are NAs?
table(simData$model2gen[is.na(simData$LOOIC)])
table(simData$model2fit[is.na(simData$LOOIC)])


## Model fit criteria:
densityplot(simData$WAIC)
densityplot(simData$LOOIC)

## Convergence checks:
densityplot(simData$maxRhat)
densityplot(simData$minNeff)
table(simData$divergent)
table(simData$maxTree)
densityplot(simData$minEnergy)

# ======================================================================================================== #
#### 02a) Compute forward matrix: ####

## Option A: Use all models:
selModVec <- c(1:12)

## Option B: Use only selected models:
# selModVec <- c(1, 2, 4, 6, 12)
# selModVec <- c(1:4, 6, 12)

## Select data:
aggrData <- subset(simData, model2gen %in% selModVec & model2fit %in% selModVec)
unique(aggrData$model2gen)
unique(aggrData$model2fit)
# table(aggrData$iteration)

## Initialize:
model2genVec <- sort(unique(aggrData$model2gen)); nModGenSel <- length(model2genVec)
model2fitVec <- sort(unique(aggrData$model2fit)); nModFitSel <- length(model2fitVec)
forwardMat <- matrix(NA, nrow = nModGenSel, ncol = nModFitSel)
for (iModGen in 1:nModGenSel){ # model2gen <- 12

  model2gen <- model2genVec[iModGen]
  
  cat(paste0(">>> Start extracting for model M", model2gen, "\n"))
  
  ## Winning model for data generated with model2gen:
  winModVec <- rep(NA, nIter) # iIter <- 1
  
  ## Determine best fitting model for each data set:
  for (iIter in 1:nIter){ # model2fit <- 1
    
    ## Locate where statistics for this data set of this generating model are stored:
    rowIdx <- which(aggrData$model2gen == model2gen & aggrData$iteration == iIter)
    if (length(rowIdx) == 0){stop(paste0("No fit yet for model2gen = ", model2gen, "; iteration = ", iIter, "\n"))}
  
    ### Determine winning model for this data set:
    ## WAIC:
    # minVal <- min(aggrData$WAIC[rowIdx], na.rm = T)
    # minIdx <- which(aggrData$WAIC[rowIdx] == minVal)
    
    ## LOOIC:
    minVal <- min(aggrData$LOOIC[rowIdx], na.rm = T)
    minIdx <- which(aggrData$LOOIC[rowIdx] == minVal)
    
    ## Save ID of winning model:
    minIdx <- minIdx[length(minIdx)] # if multiple: take last one
    winModVec[iIter] <- aggrData$model2fit[rowIdx[minIdx]]
  
  }
  
  ## For each fitted model, count how often it won for this generating model: 
  cat(paste0(">>> Compute p(fitted|given) for model M", model2gen, "\n"))
  for (iModFit in 1:nModFitSel){ # model2fit <- 1
    model2fit <- model2fitVec[iModFit]
    forwardMat[iModGen, iModFit] <- sum(winModVec == model2fit)/nIter
  }
}

rowSums(forwardMat) # should be all 1
stopifnot(all(rowSums(forwardMat) == 1))
# colSums(forwardMat)

## Inspect diagonal:
diagVec <- diag(forwardMat); diagVec
round(t(stat.desc(diagVec)), 2)

## Prepare matrix:
forwM <- round(forwardMat, 2) 
rownames(forwM) <- paste0("M", selModVec)
colnames(forwM) <- paste0("M", selModVec)

## Plot:
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# https://cran.r-hub.io/web/packages/corrplot/vignettes/corrplot-intro.html
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
corrplot::corrplot(forwM, # type = "upper",
                   addCoef.col = 'black', # darkgrey, grey50, grey60, orange
                   number.cex = 1.2, tl.cex = 2, is.corr = FALSE, col.lim = c(0, 1),
                   col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
corrplot::corrplot(forwM, method = "color",
                   addCoef.col = 'grey50', # darkgrey, grey50, grey60, orange
                   number.cex = 1.2, tl.cex = 2, is.corr = FALSE, col.lim = c(0, 1),
                   col = viridis(100), tl.col = "black", tl.pos = "lt")

png(paste0(dirs$plots, "model_recovery_forwardMatrix_nIter", nIter, ".png"), width = 480, height = 480)
png(paste0(dirs$plots, "model_recovery_forwardMatrix_nIter", nIter, "_", paste0(paste0("M", selModVec), collapse = ""), ".png"), width = 480, height = 480)
dev.off()

# ======================================================================================================== #
#### 02b) Compute inverse matrix: ####

## Initialize:
forwardMatCount <- forwardMat * nIter # convert back in raw counts
model2genVec <- sort(unique(aggrData$model2gen)); nModGenSel <- length(model2genVec)
model2fitVec <- sort(unique(aggrData$model2fit)); nModFitSel <- length(model2fitVec)
inverseMat <- matrix(NA, nrow = nModGenSel, ncol = nModFitSel)
for (iModFit in 1:nModFitSel){ # iModFit <- 1

  ## Among all data sets where model2fit wins: how often was it each generative model?
  inverseMat[, iModFit] <- forwardMatCount[, iModFit] / sum(forwardMatCount[, iModFit])
}

colSums(inverseMat) # should be all 1
stopifnot(all(round(colSums(inverseMat), 6) == 1))

## Inspect diagonal:
diagVec <- diag(inverseMat); diagVec
round(t(stat.desc(diagVec)), 2)

## Prepare matrix:
invM <- round(inverseMat, 2) 
rownames(invM) <- paste0("M", selModVec)
colnames(invM) <- paste0("M", selModVec)

## Plot:
corrplot::corrplot(invM, # type = "upper",
                   addCoef.col = 'black', # darkgrey, grey50, grey60, orange
                   number.cex = 1.2, tl.cex = 2, is.corr = FALSE, col.lim = c(0, 1),
                   col = COL1("YlOrRd"), tl.col = "black", tl.pos = "lt") # YlOrRd brewer.pal(9, "YlOrRd")
corrplot::corrplot(invM, method = "color",
                   addCoef.col = 'grey50', # darkgrey, grey50, grey60, orange
                   number.cex = 1.2, tl.cex = 2, is.corr = FALSE, col.lim = c(0, 1),
                   col = viridis(100), tl.col = "black", tl.pos = "lt")

png(paste0(dirs$plots, "model_recovery_inverseMatrix_nIter", nIter, ".png"), width = 480, height = 480)
png(paste0(dirs$plots, "model_recovery_inverseMatrix_nIter", nIter, "_", paste0(paste0("M", selModVec), collapse = ""), ".png"), width = 480, height = 480)
dev.off()

# ======================================================================================================== #
#### 03) Compute permutation distribution: ####

# --------------------------------------------------------- #
### Option A: mimick entire process:

model2genVec <- sort(unique(aggrData$model2gen)); nModGenSel <- length(model2genVec)
model2fitVec <- sort(unique(aggrData$model2fit)); nModFitSel <- length(model2fitVec)

## Initialize:
nIterPerm <- 100
cat(paste0("Compute permutation null distribution for matrix of size ", nModGenSel, " x ", nModFitSel, " with ", nIterPerm, " iterations\n"))
maxdiagPermVec <- rep(NA, nIterPerm)
diagPermMat <- matrix(NA, nrow = nIterPerm, ncol = nModFitSel)

## Start iterations:
start.time.loop <- Sys.time()
for (iIterPerm in 1:nIterPerm){
  
  cat(paste0("Start iteration ", iIterPerm, "\n"))
  
  ## Initialize:
  permMat <- matrix(NA, nrow = nModGenSel, ncol = nModFitSel)
  
  for (iModGen in 1:nModGenSel){ # model2gen <- 12
    
    model2gen <- model2genVec[iModGen]
    
    ## Winning model for data generated with model2gen:
    winModVec <- rep(NA, nIter) # iIter <- 1
    
    ## Determine best fitting model for each data set:
    for (iIter in 1:nIter){ # model2fit <- 1
  
      ## Identify rows:
      rowIdx <- which(aggrData$model2gen == model2gen & aggrData$iteration == iIter)
      if (length(rowIdx) == 0){stop(paste0("No fit yet for model2gen = ", model2gen, "; iteration = ", iIter, "\n"))}
      
      ## Sample:
      tmp <- sample(aggrData$LOOIC[rowIdx], length(rowIdx), replace = F) # extract, randomly permute

      ## Find minimum:
      minVal <- min(tmp, na.rm = T)
      minIdx <- which(tmp == minVal)
  
      ## Save ID of winning model:
      minIdx <- minIdx[length(minIdx)] # if multiple: take last one
      winModVec[iIter] <- minIdx
  
    } # end iIter
    
    ## For each fitted model, count how often it won for this generating model: 
    for (iModFit in 1:nModFitSel){ # model2fit <- 1
      model2fit <- model2fitVec[iModFit]
      permMat[iModGen, iModFit] <- sum(winModVec == model2fit)/nIter
    } # end model2fit
  } # end model2gen
    
  ## Save maximum on-diagonal element:
  maxdiagPermVec[iIterPerm] <- max(diag(permMat), na.rm = T)
  diagPermMat[iIterPerm, ] <- diag(permMat)
} # end iIterPerm
end.time.loop <- Sys.time()
duration <- difftime(end.time.loop, start.time.loop); duration

# --------------------------------------------------------- #
### Option B:
nIterPerm <- 1000
maxPermVec <- rep(NA, nIterPerm)

for (iIterPerm in 1:nIterPerm){
  for (model2gen in 1:nMod2Sim){ # model2gen <- 12
    ## For each fitted model, count how often it won for this generating model: 
    for (model2fit in 1:nMod2Fit){ # model2fit <- 1
      permMat[model2gen, model2fit] <- sum(winModVec == model2fit)/nIter
    } # end model2fit
  } # end model2gen
} # end iIterPerm


# END OF FILE.
