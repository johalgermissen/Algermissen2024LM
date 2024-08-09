#!/usr/bin/env Rscript
# ============================================================================ #
## 02_mgngstakes_regression.R
## MGNGStakes regression models on behaviour.
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
#### 02) Check learning: ####

nSub <- length(unique(data$subject_n))
pVec <- rep(NA, nSub)
for (iSub in 1:nSub){ # iSub <- 20
  subData <- subset(data, subject_n == iSub)
  mod <- glm("response_n ~ reqAction_f", data = subData, family = "binomial")
  pVec[iSub] <- summary(mod)$coefficients[2, 4]
  cat(paste0("Subject ", iSub, ": p = ", pVec[iSub], "\n"))
} 
which (pVec >= 0.05) # 52

# ============================================================================ #
#### 03) Means & SDs per condition: ####

out <- aggregate_sub_cond(modData, yVar = "response_n", xVarVec = "condition_f", subVar = "subject_f",
                          format = "wide", printSumStats = T)
out <- aggregate_sub_cond(modData, yVar = "ACC_n", xVarVec = "condition_f", subVar = "subject_f",
                          format = "wide", printSumStats = T)
out <- aggregate_sub_cond(modData, yVar = "RTcleaned_n", xVarVec = "condition_f", subVar = "subject_f",
                          format = "wide", printSumStats = T)

out <- aggregate_sub_cond(modData, yVar = "response_n", xVarVec = "condition_stakes_f", subVar = "subject_f",
                          format = "wide", printSumStats = T)
out <- aggregate_sub_cond(modData, yVar = "ACC_n", xVarVec = "condition_stakes_f", subVar = "subject_f",
                          format = "wide", printSumStats = T)
out <- aggregate_sub_cond(modData, yVar = "RTcleaned_n", xVarVec = "condition_stakes_f", subVar = "subject_f",
                          format = "wide", printSumStats = T)

# ============================================================================ #
#### 04) Fit mixed-effects models: ####

# ---------------------------------------------------------------------------- #
#### 04a) Response: ####

## Main effects:
formula <- "response_n ~ reqAction_f + (reqAction_f|subject_f)"
formula <- "response_n ~ valence_f + (valence_f|subject_f)"

## 2-way interactions:
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"

formula <- "response_n ~ condition_f * stakes_f + (condition_f * stakes_f|subject_f)"

## 3-way interaction:
formula <- "response_n ~ reqAction_f * valence_f * stakes_f + (reqAction_f * valence_f * stakes_f|subject_f)"

## Read existing model back in:
mod <- fit_lmem(formula)
quickCI(mod, nRound = 3)
mod <- fit_lmem(formula, useLRT = T)

### Fit logistic regression:
mod <- glmer(formula = formula, data = modData, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)
quickCI(mod, nRound = 3)
mod_LRT <- mixed(formula = formula, data = modData, method = "LRT", type = "III", family = binomial(),
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
anova(mod_LRT)

# ---------------------------------------------------------------------------- #
#### 04b) Accuracy: ####

## Main effects:
formula <- "ACC_n ~ reqAction_f + (reqAction_f|subject_f)"
formula <- "ACC_n ~ valence_f + (valence_f|subject_f)"
formula <- "ACC_n ~ reqCongruency_f + (reqCongruency_f|subject_f)"

## 2-way interactions:
formula <- "ACC_n ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"
formula <- "ACC_n ~ stakes_f * reqCongruency_f + (stakes_f * reqCongruency_f|subject_f)"

formula <- "ACC_n ~ condition_f * stakes_f + (condition_f * stakes_f|subject_f)"

formula <- "ACC_n ~ stakes_f * cueRep_z + (stakes_f * cueRep_z|subject_f)"
formula <- "ACC_n ~ stakes_f * trialnr_z + (stakes_f * trialnr_z|subject_f)"

## 3-way interactions:
formula <- "ACC_n ~ reqAction_f * valence_f * stakes_f + (reqAction_f * valence_f * stakes_f|subject_f)"

## Read existing model back in:
mod <- fit_lmem(formula)
quickCI(mod, nRound = 3)
mod <- fit_lmem(formula, useLRT = T)

### Fit logistic regression:
mod <- glmer(formula = formula, data = modData, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)
quickCI(mod, nRound = 3)
mod_LRT <- mixed(formula = formula, data = modData, method = "LRT", type = "III", family = binomial(),
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
anova(mod_LRT)

# ---------------------------------------------------------------------------- #
#### 04c) RTs: ####

## Main effects:
formula <- "RTcleaned_z ~ reqAction_f + (reqAction_f|subject_f)"
formula <- "RTcleaned_z ~ valence_f + (valence_f|subject_f)"
formula <- "RTcleaned_z ~ reqCongruency_f + (reqCongruency_f|subject_f)"
formula <- "RTcleaned_z ~ actCongruency_f + (actCongruency_f|subject_f)"
formula <- "RTcleaned_z ~ stakes_f + (stakes_f|subject_f)"

## 2-way interactions:
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"

formula <- "RTcleaned_z ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"

formula <- "RTcleaned_z ~ actCongruency_f * stakes_f + (actCongruency_f * stakes_f|subject_f)"
formula <- "RTcleaned_z ~ valence_f * stakes_f + (valence_f * stakes_f|subject_f)"

formula <- "RTcleaned_z ~ condition_f * stakes_f + (condition_f * stakes_f|subject_f)"

formula <- "RTcleaned_z ~ stakes_f * cueRep_z + (stakes_f * cueRep_z|subject_f)"
formula <- "RTcleaned_z ~ stakes_f * trialnr_z + (stakes_f * trialnr_z|subject_f)"

## 3-way interaction:
formula <- "RTcleaned_z ~ reqAction_f * valence_f * stakes_f + (reqAction_f * valence_f * stakes_f|subject_f)"
formula <- "RTcleaned_z ~ stakes_f * reqAction_f * valence_f + (stakes_f * reqAction_f * valence_f|subject_f)"

## Unstandardized RTs:
formula <- "RTcleaned_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
formula <- "RTcleaned_n ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"
# formula <- "RT_z ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"

## Uncleaned RTs:
formula <- "RT_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
formula <- "RT_z ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"

mod <- fit_lmem(formula)
quickCI(mod, nRound = 3)
mod <- fit_lmem(formula, useLRT = T)

### Fit linear regression:
mod <- lmer(formula = formula, data = modData,
             control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)
quickCI(mod, nRound = 3)
mod_LRT <- mixed(formula = formula, data = modData, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
anova(mod_LRT)

# ============================================================================ #
#### 05) Follow-up with emmeans: ####
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/

emmeans(mod, ~ stakes_f)
emmeans(mod, specs = pairwise ~ stakes_f, 
        regrid = "response")

emmeans(mod, specs = pairwise ~ stakes_f | condition_f, 
        regrid = "response", interaction = "pairwise")
emmeans(mod, specs = pairwise ~ stakes_f | condition_f, 
        interaction = "pairwise")

emmeans(mod, specs = pairwise ~ condition_f | stakes_f, 
        regrid = "response", interaction = "pairwise")


emmeans(mod, ~ DBS_f | effort)
emmeans(mod, ~ DBS_f * valence_f)
emmeans(mod, specs = pairwise ~ DBS_f | reqAction_f*valence_f, 
        trans = "response", interaction = "pairwise")
emmeans(mod, specs = pairwise ~ DBS_f | condition_f, trans = "response", interaction = "pairwise")

# ============================================================================ #
#### 06) Plot effects with plot(effect()): ####

## Main effect:
plot(effect("reqAction_f", mod)) # 
plot(effect("valence_f", mod)) # 
plot(effect("stakes_f", mod)) # 
plot(effect("reqCongruency_f", mod)) # 

## 2-way interactions:
plot(effect("reqAction_f:valence_f", mod)) # 
plot(effect("reqAction_f:valence_f", mod, x.var = "valence_f")) # 

plot(effect("reqCongruency_f:stakes_f", mod)) # 
plot(effect("reqCongruency_f:stakes_f", mod, x.var = "reqCongruency_f")) # 
plot(effect("reqCongruency_f:stakes_f", mod, x.var = "stakes_f")) # 

plot(effect("condition_f:stakes_f", mod)) # 
plot(effect("condition_f:stakes_f", mod, x.var = "stakes_f")) # 

plot(effect("stakes_f:trialnr_z", mod)) # 
plot(effect("stakes_f:trialnr_z", mod, x.var = "stakes_f")) # 
plot(effect("stakes_f:cueRep_z", mod)) # 
plot(effect("stakes_f:cueRep_z", mod, x.var = "stakes_f")) # 
plot(effect("reqCongruency_f:stakes_f:trialnr_z", mod)) # 
plot(effect("reqCongruency_f:stakes_f:cueRep_z", mod)) # 

plot(effect("outcome_last_f:reqCongruency_f", mod)) # 
plot(effect("outcome_last_f:reqCongruency_f", mod, x.var = "reqCongruency_f")) # 

## 3-way interactions:
plot(effect("valence_f:reqResp_f:outcome_lag1_f", mod, x.var = "outcome_lag1_f")) # 

# ============================================================================ #
#### 07) Plot effect with custom_regressionbar1: ####

plot(effect("reqAction_f", mod)) # 
plot(effect("valence_f", mod)) # 
plot(effect("reqCongruency_f", mod)) # 

# ---------------------------------------------- #
## response_n ~ reqAction_f
p <- custom_regressionbar1(mod, selEff = "reqAction_f", # 
                      yLim = c(0, 1),
                      # yLim = c(0.4, 0.6),
                      selCol = colours$redblue, xLabels = c("Go", "NoGo"), 
                      xLab = "Req. Action", yLab = "Response")
png(paste0(dirs$plotDir, "regressionbar_response~reqAction.png"), width = 480, height = 480)
print(p)
dev.off()

## response_n ~ valence_f
p <- custom_regressionbar1(mod, selEff = "valence_f", # 
                      yLim = c(0, 1),
                      # yLim = c(0.4, 0.6),
                      selCol = colours$greenred, xLabels = c("Win", "Avoid"), 
                      xLab = "Valence", yLab = "Response")
png(paste0(dirs$plotDir, "regressionbar_response~valence.png"), width = 480, height = 480)
print(p)
dev.off()

# ---------------------------------------------- #
## RTcleaned_n ~ reqAction_f
p <- custom_regressionbar1(mod, selEff = "reqAction_f", # 
                      yLim = c(0.45, 0.8),
                      # yLim = c(0, 1),
                      selCol = colours$redblue, xLabels = c("Go", "NoGo"), 
                      xLab = "Req. Action", yLab = "RT (in sec.)")
png(paste0(dirs$plotDir, "regressionbar_RTcleaned~reqAction.png"), width = 480, height = 480)
print(p)
dev.off()

## RTcleaned_n ~ valence_f
p <- custom_regressionbar1(mod, selEff = "valence_f", # 
                      yLim = c(0.45, 0.8),
                      # yLim = c(0, 1),
                      selCol = colours$greenred, xLabels = c("Win", "Avoid"), 
                      xLab = "Valence", yLab = "RT (in sec.)")
png(paste0(dirs$plotDir, "regressionbar_RTcleaned~valence.png"), width = 480, height = 480)
print(p)
dev.off()

# ---------------------------------------------- #
## ACC_n ~ reqCongruency_f
p <- custom_regressionbar1(mod, selEff = "reqCongruency_f", # 
                      yLim = c(0, 1),
                      # yLim = c(0.4, 0.6),
                      selCol = colours$greenpink1, xLabels = c("cong.", "incong."), 
                      xLab = "Congruency", yLab = "Accuracy")
png(paste0(dirs$plotDir, "regressionbar_ACC~reqCongruency.png"), width = 480, height = 480)
print(p)
dev.off()

custom_regressionbar1(mod, selEff = "stakes_f", # 
                      yLim = c(0, 1),
                      # yLim = c(0.4, 0.6),
                      selCol = c("orange", "purple"), xLabels = c("high", "low"), 
                      xLab = "Stakes", yLab = "Accuracy")
png(paste0(dirs$plotDir, "regressionbar_ACC~stakes.png"), width = 480, height = 480)
print(p)
dev.off()

# ---------------------------------------------- #
## RTcleaned_n ~ reqCongruency_f
p <- custom_regressionbar1(mod, selEff = "reqCongruency_f", # 
                      yLim = c(0.45, 0.80),
                      # yLim = c(0, 1),
                      selCol = colours$greenpink1, xLabels = c("cong.", "incong."), 
                      xLab = "Congruency", yLab = "RT (in sec.)")
png(paste0(dirs$plotDir, "regressionbar_RTcleaned~reqCongruency.png"), width = 480, height = 480)
print(p)
dev.off()

p <- custom_regressionbar1(mod, selEff = "stakes_f", # 
                      yLim = c(0.45, 0.80),
                      # yLim = c(0, 1),
                      selCol = c("orange", "purple"), xLabels = c("high", "low"), 
                      xLab = "Stakes", yLab = "RT (in sec.)")
png(paste0(dirs$plotDir, "regressionbar_RTcleaned~reqCongruency.png"), width = 480, height = 480)
print(p)
dev.off()

# ============================================================================ #
#### 08) Plot coefficients with custom_coefplot: ####

nCoef <- length(fixef(mod)) - 1; nCoef

## Create colours:
# selCol <- rev(jet.colors(length(fixef(mod)) - 1, alpha = 1)); colName <- "Jet" # red first

## Select colours:
library(MetBrewer)
selCol <- met.brewer("Demuth", n = nCoef, type = "continuous"); colName <- "Demuth"
selCol <- met.brewer("Isfahan1", n = 8, type = "continuous"); selCol <- selCol[c(2, 5, 7)]; colName <- "Isfahan1"

## Display:
library(scales)
show_col(selCol)

# -------------------------------------------- #
## Just visualize with small font size:

custom_coefplot(mod, plotSub = F, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol, FTS = 10)
custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol, FTS = 10)

# -------------------------------------------- #
## Plot with saving:

## Specify xLim manually?
min(coef(mod)[[1]][, 2:ncol(coef(mod)[[1]])])
max(coef(mod)[[1]][, 2:ncol(coef(mod)[[1]])])

## Without subs:
plotSub <- F
p <- custom_coefplot(mod, plotSub = plotSub, plotText = T, dropIntercept = T, 
                     # xLim = c(-2, 1),
                     revOrder = T, selCol = selCol)
plotHandle <- formula2handle(formula)
plotName <- paste0("coef_glmer_", plotHandle, "_", colName)
png(paste0(dirs$plotDir, plotName, ".png"), width = 900, height = 480)
print(p)
dev.off()

## With subs:
plotSub <- T
p <- custom_coefplot(mod, plotSub = plotSub, plotText = T, dropIntercept = T, 
                     # xLim = c(-2, 1),
                     revOrder = T, selCol = selCol)
plotHandle <- formula2handle(formula)
plotName <- paste0("coef_glmer_", plotHandle, "_", colName)
if (plotSub){plotName <- paste0(plotName, "_subs")}
png(paste0(dirs$plotDir, plotName, ".png"), width = 900, height = 480)
print(p)
dev.off()

# ---------------------------------------------------------------------------- #
#### Plot regressor inter-correlations: ####

corrplot_regressors(mod, perSub = F, savePNG = T) # intercorrelations overall
corrplot_regressors(mod, perSub = T, savePNG = T) # intercorrelations per subject averaged

# ---------------------------------------------------------------------------- #
#### Plot coefficient inter-correlations: ####

corrplot_coefficients(mod, savePNG = T) # plot inter-correlations between coefficients

# ---------------------------------------------------------------------------- #
#### Save coefficients per subject: ####

save_coefficients(mod) # save coefficients per subject
cat("Saved everything :-)\n")

# END
