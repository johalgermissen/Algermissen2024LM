#!/usr/bin/env Rscript
# ============================================================================ #
## MGNGStakes plots of behaviour.
# Johannes Algermissen, 2023.

rm(list = ls())

# ============================================================================ #
#### Set directories, load packages and custom functions: ####

## Set codeDir:
codeDir    <- "/project/2420133.01/MGNGStakes/2019_mgngstakes_johalg/regression/"
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs.R")) # Load packages and options settings
dirs <- set_dirs(rootDir)

## Load colours:
source(paste0(helperDir, "load_colours.R")) # Load colours
colours <- load_colours()

## Load packages:
source(paste0(helperDir, "package_manager.R")) # Load packages and options settings

# ------------------------------------------------- #
## Load custom functions:

source(paste0(dirs$codeDir, "functions/00_mgngstakes_functions_regression.R")) # Load functions

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### Read behavior: ####

## Find all raw data files:
data <- read_behavior()

# ----------------------------------- #
## Preprocessing:

data <- wrapper_preprocessing(data)
# head(data)
cat("Finished pre-processing :-)\n"); # beep()

names(data)
table(data$subject_n)

# ----------------------------------- #
## Select data, standardize variables:

modData <- select_standardize(data)
plotData <- modData
table(plotData$subject_f)
cat("Finished standardization :-)\n")

# ============================================================================ #
#### SELECTION OF PLOTS FOR PAPER: ####

plotData <- modData

custom_lineplot_gg(plotData, xVar = "cueRep_n", yVar = "response_n", zVar = "condition_f", subVar = "subject_f", 
                   selCol = c("#007174", "#FF654E", "#007174", "#FF654E"), selLineType = c(1, 1, 2, 2), breakVec = c(1, 5, 10, 15, 20), addLegend = F,
                   savePNG = T)

custom_barplot2(plotData, xVar = "valence_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")

custom_barplot2(plotData, xVar = "stakes_f", yVar = "ACC_n", zVar = "reqCongruency_short2_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f")

yLimRT <- c(0.4, 1.0)
custom_barplot2(plotData, xVar = "valence_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)

custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "stakes_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)

## Densityplots:
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "stakes_f", addLegend = F) 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "stakes_f", addLegend = T) 

# ============================================================================ #
#### DENSITY OF RTs: ####

p <- customplot_density2(plotData, xVar = "RT_n", zVar = "reqAction_f") 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "valence_f") 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "stakes_f") 

p <- customplot_density2(plotData, xVar = "RT_n", zVar = "reqAction_f", addLegend = T) 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "valence_f", addLegend = T) 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "stakes_f", addLegend = T) 

# ============================================================================ #
#### CUSTOM BAR PLOTS averaged over time: ####

# ---------------------------------------------------------------------------- #
#### Bar plots RESPONSE: ####

plotData <- modData

## 1D bar plot: Response per block:
custom_barplot1(plotData, xVar = "block_f", yVar = "response_n", subVar = "subject_f")
custom_barplot1(plotData, xVar = "cueRep_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per required action per required action/ valence/ stakes:
custom_barplot2(plotData, xVar = "block_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "response_n", zVar = "stakes_f", subVar = "subject_f")

custom_barplot2(plotData, xVar = "cueRep_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "response_n", zVar = "stakes_f", subVar = "subject_f", isPoint = F, isConnect = F)

## 1D bar plot: Response per required action:
custom_barplot1(plotData, xVar = "reqAction_f", yVar = "response_n", subVar = "subject_f")

## 1D bar plot: Response per valence:
custom_barplot1(plotData, xVar = "valence_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per required action per valence:
custom_barplot2(plotData, xVar = "valence_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")

## 1D bar plot: Response per condition:
custom_barplot1(plotData, xVar = "condition_f", yVar = "response_n", subVar = "subject_f")

## 1D bar plot: Response per stakes:
custom_barplot1(modData, xVar = "stakes_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per stakes per required action/valence:
custom_barplot2(plotData, xVar = "stakes_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "stakes_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "stakes_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "valence_f", yVar = "response_n", zVar = "stakes_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "condition_f", yVar = "response_n", zVar = "stakes_f", subVar = "subject_f")

## 1D bar plot: Response per congruency (required action):
custom_barplot1(plotData, xVar = "reqCongruency_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per congruency (required action) per stakes:
custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "response_n", zVar = "stakes_f", subVar = "subject_f")

## 1D bar plot: Response per outcome last trial:
custom_barplot1(plotData, xVar = "outcome_lag1_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per required action per outcome last trial:
custom_barplot2(plotData, xVar = "outcome_lag1_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "outcome_lag1_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "outcome_lag1_short_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "valence_f", yVar = "response_n", zVar = "outcome_lag1_short_f", subVar = "subject_f")

## 3-way interaction: response per required action x valence x stakes:
plotData <- subset(modData, !is.na(response_n) & !is.na(reqAction_f) & !is.na(valence_f) & !is.na(stakes_f))
ggplot(plotData, aes(y = response_n, x = reqAction_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(stakes_f))
custom_barplot3(data = modData, yVar = "response_n", xVar = "reqAction_f", zVar = "valence_f", splitVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)

# ---------------------------------------------------------------------------- #
#### Bar plots ACCURACY: ####

plotData <- modData

## 1D bar plot: Accuracy per block:
custom_barplot1(plotData, xVar = "block_f", yVar = "ACC_n", subVar = "subject_f")
custom_barplot1(plotData, xVar = "cueRep_f", yVar = "ACC_n", subVar = "subject_f")

## 2D bar plot: Response per required action per required action/ valence/ stakes:
custom_barplot2(plotData, xVar = "block_f", yVar = "ACC_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "ACC_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "ACC_n", zVar = "reqCongruency_short_f", subVar = "subject_f")

custom_barplot2(plotData, xVar = "cueRep_f", yVar = "ACC_n", zVar = "reqAction_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "ACC_n", zVar = "valence_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f", isPoint = F, isConnect = F)

## 1D bar plot: Accuracy per required action:
custom_barplot1(plotData, xVar = "reqAction_f", yVar = "ACC_n", subVar = "subject_f")

## 1D bar plot: Accuracy per valence:
custom_barplot1(plotData, xVar = "valence_f", yVar = "ACC_n", subVar = "subject_f")

## 2D bar plot: Accuracy per required action per valence:
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "ACC_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "valence_f", yVar = "ACC_n", zVar = "reqAction_f", subVar = "subject_f")

## 1D bar plot: Accuracy per condition:
custom_barplot1(plotData, xVar = "condition_f", yVar = "response_n", subVar = "subject_f")

## 1D bar plot: Accuracy per stakes:
custom_barplot1(plotData, xVar = "stakes_f", yVar = "ACC_n", subVar = "subject_f")

## 2D bar plot: Accuracy per stakes per required action/ valence/ congruency/ condition:
custom_barplot2(plotData, xVar = "stakes_f", yVar = "ACC_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "stakes_f", yVar = "ACC_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "valence_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f")

custom_barplot2(plotData, xVar = "stakes_f", yVar = "ACC_n", zVar = "reqCongruency_short2_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f", fontSize = 28)

## 1D bar plot: Accuracy per congruency (required action):
custom_barplot1(plotData, xVar = "reqCongruency_f", yVar = "ACC_n", subVar = "subject_f")

## 3-way interaction: Accuracy per required action x valence x stakes:
plotData <- subset(modData, !is.na(ACC_n) & !is.na(reqAction_f) & !is.na(valence_f) & !is.na(stakes_f))
ggplot(plotData, aes(y = ACC_n, x = reqAction_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(stakes_f))
custom_barplot3(data = modData, yVar = "ACC_n", xVar = "reqAction_f", zVar = "valence_f", splitVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)

# ---------------------------------------------------------------------------- #
#### Bar plots RTs (raw): ####

plotData <- modData

# yLimRT <- c(0, 1.0)
yLimRT <- c(0.4, 1.0)

## 1D bar plot: RT per block:
custom_barplot1(plotData, xVar = "block_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)
custom_barplot1(plotData, xVar = "cueRep_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per congruency (required action) per block:
custom_barplot2(plotData, xVar = "block_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "block_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "block_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short_f", subVar = "subject_f", yLim = yLimRT)

custom_barplot2(plotData, xVar = "cueRep_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", isPoint = F, isConnect = F)

## 1D bar plot: RT per required action:
custom_barplot1(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 1D bar plot: RT per valence:
custom_barplot1(plotData, xVar = "valence_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per valence per required action:
custom_barplot2(plotData, xVar = "valence_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)

## 1D bar plot: RT per cue condition:
custom_barplot1(modData, xVar = "condition_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 1D bar plot: RT per stakes:
custom_barplot1(plotData, xVar = "stakes_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per required action per stakes:
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "stakes_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "valence_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "stakes_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)

custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "stakes_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT, fontSize = 28)

## 1D bar plot: RT per congruency (required action):
custom_barplot1(modData, xVar = "reqCongruency_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

# -------------------------------------------- #
## 1D bar plot: RT per outcome last trial:
custom_barplot1(modData, xVar = "outcome_last_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)
custom_barplot1(modData, xVar = "outcome_last_rel_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per required action per outcome last trial:
custom_barplot2(modData, xVar = "outcome_last_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(modData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "outcome_last_short_f", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per valence per outcome last trial:
# custom_barplot2(modData, xVar = "outcome_last_short_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)
# custom_barplot2(modData, xVar = "valence_f", yVar = "RTcleaned_n", zVar = "outcome_last_short_f", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per congruency (required action) per outcome last trial:
custom_barplot2(modData, xVar = "outcome_last_short_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short1_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(modData, xVar = "reqCongruency_short1_f", yVar = "RTcleaned_n", zVar = "outcome_last_short_f", subVar = "subject_f", yLim = yLimRT)

custom_barplot2(modData, xVar = "outcome_last_rel_short_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short1_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(modData, xVar = "reqCongruency_short1_f", yVar = "RTcleaned_n", zVar = "outcome_last_rel_short_f", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per repeat/switch per outcome last trial:
custom_barplot2(modData, xVar = "outcome_last_short_f", yVar = "RTcleaned_n", zVar = "repeat_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(modData, xVar = "repeat_f", yVar = "RTcleaned_n", zVar = "outcome_last_short_f", subVar = "subject_f", yLim = yLimRT)

## 3-way interaction: RT per required action x valence x stakes:
plotData <- subset(modData, !is.na(RTcleaned_n) & !is.na(reqAction_f) & !is.na(valence_f) & !is.na(stakes_f))
ggplot(plotData, aes(y = RTcleaned_n, x = reqAction_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(stakes_f))
custom_barplot3(data = modData, yVar = "RTcleaned_n", xVar = "reqAction_f", zVar = "valence_f", splitVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)

# ---------------------------------------------------------------------------- #
#### Bar plots log(RTs): ####

yLimRTLog <- c(4, 7)

## 1D bar plot: log(RT) per emotion:
custom_barplot1(modData, xVar = "valence_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 1D bar plot: log(RT) per valence:
custom_barplot1(modData, xVar = "valence_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 1D bar plot: log(RT) per required action:
custom_barplot1(modData, xVar = "reqAction_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 1D bar plot: log(RT) per congruency (required action):
custom_barplot1(modData, xVar = "reqCongruency_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 1D bar plot: log(RT) per congruency (actual response):
custom_barplot1(modData, xVar = "actCongruency_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 2D bar plot: log(RT) per required action per emotion:
custom_barplot2(modData, xVar = "valence_f", yVar = "RT_ms_log_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRTLog)

# END
