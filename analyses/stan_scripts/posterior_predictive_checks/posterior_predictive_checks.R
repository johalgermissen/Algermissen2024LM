# ==================================================================================== #
#' posterior_predictive_checks.R
#' - Load raw data and simulations for given model for given settings.
#' - Average over iterations of simulations, assign to data frame.
#' - Create plots of simulated responses/ RTs given task factors.
## Johannes Algermissen, 2023.

cat("\014")
rm(list = ls())

# ==================================================================================== #
#### Load packages: ####

require(ggplot2)
require(lattice)
require(ltm)
require(pastecs)
require(stringr)

require(lme4)
require(afex)
require(effects)

require(rstan)
require(bayesplot)
require(loo)

options(scipen = 20)

# ==================================================================================== #
#### 00a) Select model, suffix, simType: ####

## Model ID:
model2fit <- 12
model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")

## Model type:
suffix <- "_sepParam"

## Simulation type:
simType <- "modSim"

## Check inputs:
cat(paste0("Model to fit is ", model2fit_str), "\n")
cat(paste0("Suffix is ", suffix, "\n"))
cat(paste0("Simulation type is ", simType, "\n"))

# ==================================================================================== #
#### 00b) Set directories for stsan: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/stan_scripts/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs_stan.R")) # Load packages and options settings
stanDirs <- set_dirs_stan(rootDir)

# ==================================================================================== #
#### Load directories for regression:

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/regression/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs.R")) # Load packages and options settings
dirs <- set_dirs(rootDir)

# ------------------------------------------------- #
## Load custom functions:

source(paste0(dirs$codeDir, "functions/00_mgngstakes_functions_regression.R")) # Load functions

# ==================================================================================== #
#### 01) Load simulations: ####

simFile <- paste0(stanDirs$sims, "M", model2fit_str, suffix, "_", simType, "_simulations.Rdata")
cat(paste0("Start loading ", simFile))
load(simFile) # loads data, simQGo, simQNoGo, simResp, simRT
cat(paste0("Finished loading ", simFile, "! :-)\n"))

# ==================================================================================== #
#### 02) Assign to data: ####

## Dimensions of simulated data:
dim(simResp) # iterations x trials
dim(simRT) # iterations x trials

max(simRT, na.rm = T)
sum(simRT > 1.3, na.rm = T)
mean(simRT > 1.3, na.rm = T)
# densityplot(c(simRT)) # vector too long to be plotted

## Original data:
names(data)

## Add original responses:
data$response_n <- data$resp
data$ACC_n <- ifelse(data$reqAction == 1, data$response_n, 1 - data$response_n)
data$RT_n <- data$rt

## Compute Accuracy:
data$simACC_n <- ifelse(data$reqAction == 1, data$simResp_n, 1 - data$simResp_n)

## Average over iterations, assign to data:
# dim(simResp)
data$simResp_n <- colMeans(simResp, na.rm = T)
data$simRT_n <- colMeans(simRT, na.rm = T) # 

## Set NoGo RTs to NA:
data$rt[which(data$rt == 0)] <- NA
data$simRT_n[which(data$simRT_n == 0)] <- NA

## Set too slow RTs to NA:
threshRT <- 1.3
simResp_cleaned <- simResp # copy over
simRT_cleaned <- simRT # copy over
simRT_cleaned[simResp_cleaned == 0] <- NA # set NoGo RTs to NA
simResp_cleaned[simRT_cleaned > threshRT] <- 0 # set responses with too long RTs to 0
simRT_cleaned[simRT_cleaned > threshRT] <- NA # delete too long RTs
data$simResp_cleaned_n <- colMeans(simResp_cleaned, na.rm = T)
data$simACC_cleaned_n <- ifelse(data$reqAction == 1, data$simResp_cleaned_n, 1 - data$simResp_cleaned_n)
data$simRT_cleaned_n <- colMeans(simRT_cleaned, na.rm = T) # 

## Plot against each other:
# plot(data$simResp_n, data$simResp_cleaned_n)
# plot(data$simRT_n, data$simRT_cleaned_n)

## Standardize:
data$RT_z <- as.numeric(scale(data$RT_n))
data$simResp_z <- as.numeric(scale(data$simResp_n))
data$simResp_cleaned_z <- as.numeric(scale(data$simResp_cleaned_n))
data$simRT_z <- as.numeric(scale(data$simRT_n))
data$simRT_cleaned_z <- as.numeric(scale(data$simRT_cleaned_n))

## Add factors:
data$subject_n <- data$subject
data$subject_f <- factor(data$subject_n)

data$reqAction_n <- data$reqAction
data$reqAction_f <- factor(data$reqAction_n, levels = c(1, 0), labels = c("Go", "NoGo"))
data$valence_n <- data$valence
data$valence_f <- factor(data$valence_n, levels = c(1, 0), labels = c("Win", "Avoid"))
data$stakes_n <- data$stakes
data$stakes_f <- factor(data$stakes_n, levels = c(1, 0), labels = c("high", "low"))
data$reqCongruency_n <- data$congruency
data$reqCongruency_f <- factor(data$reqCongruency_n, levels = c(1, 0), labels = c("congruent", "incongruent"))
data$reqCongruency_short1_f <- factor(data$reqCongruency_n, levels = c(1, 0), labels = c("cong", "incong"))
data$reqCongruency_short2_f <- factor(data$reqCongruency_n, levels = c(1, 0), labels = c("con", "inc"))
data$condition_n <- 2 * (1 - data$reqAction_n) + (1 - data$valence_n) + 1
data$condition_f <- factor(data$condition_n, levels = 1:4, labels = c("G2W", "G2A", "NG2W", "NG2A"))
data$condition_short1_f <- factor(data$condition_n, levels = 1:4, labels = c("G2W", "G2A", "N2W", "N2A"))

# =================================================================================================================== #
#### 03) RT distributions: ####

xLimRT <- c(0, 1.3)

## Global distribution:
densityplot(data$RT_n, xlim = xLimRT) # original
densityplot(data$simRT_n, xlim = xLimRT) # simulated
densityplot(data$simRT_cleaned_n, xlim = xLimRT) # simulated

## Effect of valence:
densityplot(~ RT_n, data = data, groups = valence_f, auto.key = TRUE, xlim = xLimRT) # simulated
densityplot(~ simRT_n, data = data, groups = valence_f, auto.key = TRUE, xlim = xLimRT) # simulated
densityplot(~ simRT_cleaned_n, data = data, groups = valence_f, auto.key = TRUE, xlim = xLimRT) # simulated

## Effect of stakes:
densityplot(~ RT_n, data = data, groups = stakes_f, auto.key = TRUE, xlim = xLimRT) # simulated
densityplot(~ simRT_n, data = data, groups = stakes_f, auto.key = TRUE, xlim = xLimRT) # simulated
densityplot(~ simRT_cleaned_n, data = data, groups = stakes_f, auto.key = TRUE, xlim = xLimRT) # simulated


# =================================================================================================================== #
#### 04a) Plot response/ RT ~ reqAction x valence (square): ####

plotData <- data
length(unique(data$subject))

source(paste0(dirs$codeDir, "00_mgngstakes_functions_regression.R")) # Load functions

## Responses:
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "simResp_n", zVar = "valence_f", subVar = "subject_f")
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "simResp_cleaned_n", zVar = "valence_f", subVar = "subject_f")

## RTs:
yLimRT <- c(0.4, 1.0)
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "RT_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "simRT_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "simRT_cleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)

# =================================================================================================================== #
#### 04b) Plot ACC/ RT ~ stakes (square): ####

## Accuracy:
yLimpGo <- c(0, 1)
p <- custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimpGo)
p <- custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "simACC_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimpGo)
p <- custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "simACC_cleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimpGo)
p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimpGo)
p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "simACC_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimpGo)
p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "simACC_cleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimpGo)

## RTs:
yLimRT <- c(0.4, 1.0)
p <- custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "RT_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
p <- custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "simRT_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
p <- custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "simRT_cleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "RT_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "simRT_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "simRT_cleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)

# =================================================================================================================== #
#### 04c) Save (narrow): ####

yLimRT <- c(0.4, 1.0)
yLimpGo <- c(0, 1)
plotWidth <- 240
plotHeight <- 480

# ------------------------------------------------ #
### Valence:
## Responses:

yVar <- "response_n"; xVar <- "reqAction_f"; zVar <- "valence_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()
yVar <- "simResp_n"; xVar <- "reqAction_f"; zVar <- "valence_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()
yVar <- "simResp_cleaned_n"; xVar <- "reqAction_f"; zVar <- "valence_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()

## RTs:
yVar <- "RT_n"; xVar <- "reqAction_f"; zVar <- "valence_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()
yVar <- "simRT_n"; xVar <- "reqAction_f"; zVar <- "valence_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()
yVar <- "simRT_cleaned_n"; xVar <- "reqAction_f"; zVar <- "valence_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()

# ------------------------------------------------ #
### Stakes:
## Responses:
yVar <- "response_n"; xVar <- "reqCongruency_short2_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()
yVar <- "simResp_n"; xVar <- "reqCongruency_short2_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()
yVar <- "simResp_cleaned_n"; xVar <- "reqCongruency_short2_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()

yVar <- "response_n"; xVar <- "condition_short1_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()
yVar <- "simResp_n"; xVar <- "condition_short1_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()
yVar <- "simResp_cleaned_n"; xVar <- "condition_short1_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimpGo; savePlot()

## RTs:
yVar <- "RT_n"; xVar <- "reqCongruency_short2_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()
yVar <- "simRT_n"; xVar <- "reqCongruency_short2_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()
yVar <- "simRT_cleaned_n"; xVar <- "reqCongruency_short2_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()

yVar <- "RT_n"; xVar <- "condition_short1_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()
yVar <- "simRT_n"; xVar <- "condition_short1_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()
yVar <- "simRT_cleaned_n"; xVar <- "condition_short1_f"; zVar <- "stakes_f"; subVar <- "subject_f"; yLimPlot <- yLimRT; savePlot()

## Define function:
savePlot <- function(){
  ## Plot:
  p <- custom_barplot2(plotData, yVar = yVar, xVar = xVar, zVar = zVar, subVar = subVar, yLim = yLimPlot, isConnect = F)
  
  ## Save:
  figName <- paste0("custombarplot2_M", model2fit_str, "_", yVar, "~", xVar, "_", zVar, ".png")
  png(paste0(dirs$plotDir, figName), width = plotWidth, height = plotHeight)
  print(p)
  dev.off()
  print(p)  
}

# END
