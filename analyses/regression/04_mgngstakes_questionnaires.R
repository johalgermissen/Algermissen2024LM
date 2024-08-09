#!/usr/bin/env Rscript
# ============================================================================ #
## 04_mgngstakes_questionnaires.R
## MGNGStakes correlations between regression coefficients from behaviour and questionnaires.
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
#### 01) Read behavior: ####

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
# modData <- select_standardize(data, sub2excl = c())
length(unique(modData$subject_n))

# ============================================================================ #
#### 02) Read in questionnaires: ####

questData <- load_preprocess_questionnaires()
questData <- questData[which(questData$subject_n != 52), ]

# ============================================================================ #
#### 03) Fit models and plot correlations: ####

## Specify formula:
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "valence"; dvName <- "responses"; coefIdx <- 3
formula <- "ACC_n ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"; yLab <- "stakes"; dvName <- "accuracy"; coefIdx <- 3

formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "valence"; dvName <- "RTs"; coefIdx <- 3
formula <- "RTcleaned_z ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"; yLab <- "stakes"; dvName <- "RTs"; coefIdx <- 3

## Load model:
mod <- fit_lmem(formula)

## Extract coefficients, add to questData:
cat(paste0("Extract term ", names(coef(mod)$subject_f)[coefIdx], "\n"))
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]

## Select questionnaire data:
names(questData)
xVar <- "forward_span_n"; xLab <- "Forward span"
xVar <- "backward_span_n"; xLab <- "Backward span"
xVar <- "mean_span_n"; xLab <- "Memory span"
xVar <- "BIS_total_n"; xLab <- "Barratt Impulsiveness (non-planning)"
xVar <- "neuroticism_total_n"; xLab <- "Big Five Neuroticism"
# questData$neuroticism_total_n <- questData$neuroticism_total

## Plot:
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = paste0(xLab, " score"), 
                      yLab = paste0("Effect of ", yLab, " on ", dvName))

## Save:
plotName <- paste0("correlation_", dvName, "~", yLab, "_", xVar)
cat(paste0("Save as ", plotName, "\n"))
png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)
print(p)
dev.off()

# END OF FILE.
