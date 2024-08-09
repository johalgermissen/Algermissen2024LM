#!/usr/bin/env Rscript
# ============================================================================ #
## 06_mgngstakes_plots_final.R.
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
cat("Finished pre-processing :-)\n")

# ----------------------------------- #
## Select data, standardize variables:

modData <- select_standardize(data)
plotData <- modData
table(plotData$subject_f)
cat("Finished standardization :-)\n")

# ============================================================================ #
#### SELECTION OF PLOTS FOR PAPER: ####

## Global plot setting:
yLimRT <- c(0.4, 1.0)
plotHeight <- 480
plotWidth <- 480
plotWidthCoef <- 900

plotData <- modData

# ============================================================================ #
#### Figure01: #####

## Generate data:
hypData <- data.frame(condition_n <- 1:4)
hypData$congruency_f <- factor(c("cong", "cong", "incong", "incong"))
hypData$stakes_f <- factor(c("high", "low", "high", "low"))
hypData$Inv <- c(0.90, 0.85, 0.60, 0.65)
hypData$EVC <- c(0.85, 0.85, 0.70, 0.65)
hypData$EVCRT <- c(0.60, 0.60, 0.80, 0.72)

### Plotting settings:
LWD <- retrieve_plot_defaults("LWD") # 1.3
FTS <- retrieve_plot_defaults("FTS")
dodgeVal <- 0.6
colAlpha <- 1
yLim <- c(0.5, 1)
selCol <- met.brewer("Cassatt2", n = 10, type = "continuous"); selCol <- selCol[c(2, 5)]; colName <- "Cassatt2" # dark blue, blue, light blue

# ------------------------------------------------------- #
### Figure01D:
hypData$y <- hypData$Inv; yName <- "BiasInvigoration"

## Plot:
p <- ggplot(hypData, aes(x = congruency_f, fill = stakes_f, y = y)) + 
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = dodgeVal,
               lwd = LWD, color = "black") + 
  scale_fill_manual(values = selCol, limits = levels(hypData$stakes_f)) + 
  labs(x = "Congruency", y = "Accuracy", fill = "Stakes") +
  coord_cartesian(ylim = yLim) + 
  theme_classic() + 
  theme(axis.text = element_text(size = FTS),
        axis.title = element_text(size = FTS), 
        plot.title = element_text(size = FTS, hjust = 0.5), 
        legend.title = element_blank(), legend.position = "none",
        # legend.text = element_text(size = FTS),
        # legend.title = element_text(size = FTS),
        axis.line = element_line(colour = 'black')) # , linewidth = LWD)) # fixed font sizes
print(p)  
png(paste0(dirs$plot, "final/Figure01D.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure01E Left Panel (EVC-responses):
hypData$y <- hypData$EVC; yName <- "MotivationForControl"

## Plot:
p <- ggplot(hypData, aes(x = congruency_f, fill = stakes_f, y = y)) + 
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = dodgeVal,
               lwd = LWD, color = "black") + 
  scale_fill_manual(values = selCol, limits = levels(hypData$stakes_f)) + 
  labs(x = "Congruency", y = "Accuracy", fill = "Stakes") +
  coord_cartesian(ylim = yLim) + 
  theme_classic() + 
  theme(axis.text = element_text(size = FTS),
        axis.title = element_text(size = FTS), 
        plot.title = element_text(size = FTS, hjust = 0.5), 
        legend.title = element_blank(), legend.position = "none",
        # legend.text = element_text(size = FTS),
        # legend.title = element_text(size = FTS),
        axis.line = element_line(colour = 'black')) # , linewidth = LWD)) # fixed font sizes
print(p)  
png(paste0(dirs$plot, "final/Figure01ELeft.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure01E Right Panel (EVC-RTs):
hypData$y <- hypData$EVCRT; yName <- "MotivationForControl_RT"

## Plot:
p <- ggplot(hypData, aes(x = congruency_f, fill = stakes_f, y = y)) + 
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = dodgeVal,
               lwd = LWD, color = "black") + 
  scale_fill_manual(values = selCol, limits = levels(hypData$stakes_f)) + 
  labs(x = "Congruency", y = "RT (in sec.)", fill = "Stakes") +
  coord_cartesian(ylim = yLim) + 
  theme_classic() + 
  theme(axis.text = element_text(size = FTS),
        axis.title = element_text(size = FTS), 
        plot.title = element_text(size = FTS, hjust = 0.5), 
        legend.title = element_blank(), legend.position = "none",
        # legend.text = element_text(size = FTS),
        # legend.title = element_text(size = FTS),
        axis.line = element_line(colour = 'black')) # , linewidth = LWD)) # fixed font sizes
print(p)  
png(paste0(dirs$plot, "final/Figure01ERight.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ============================================================================ #
#### Figure02: #####

# ------------------------------------------------------- #
### Figure02A: Line plot response_n ~ reqAction_f * valence_f:

p <- custom_lineplot_gg(plotData, xVar = "cueRep_n", yVar = "response_n", zVar = "condition_f", subVar = "subject_f", 
                   selCol = c("#007174", "#FF654E", "#007174", "#FF654E"), selLineType = c(1, 1, 2, 2), breakVec = c(1, 5, 10, 15, 20), addLegend = F,
                   savePNG = T)
png(paste0(dirs$plot, "final/Figure02A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure02B: Bar plot response_n ~ reqAction_f * valence_f:

p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
png(paste0(dirs$plot, "final/Figure02B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure02C:

formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Demuth", n = nCoef, type = "continuous"); colName <- "Demuth"
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol)
png(paste0(dirs$plot, "final/Figure02C.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure02D:

p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f")
png(paste0(dirs$plot, "final/Figure02D.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure02E:

p <- custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "ACC_n", zVar = "stakes_f", subVar = "subject_f")
png(paste0(dirs$plot, "final/Figure02E.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure02F:

formula <- "ACC_n ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"
mod <- fit_lmem(formula)
selCol <- met.brewer("Isfahan1", n = 8, type = "continuous"); selCol <- selCol[c(2, 5, 7)]; colName <- "Isfahan1"
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol)
png(paste0(dirs$plot, "final/Figure02F.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# ============================================================================ #
#### Figure 3: #####

# ------------------------------------------------------- #
### Figure03A:

p <- customplot_density2(plotData, xVar = "RT_n", zVar = "stakes_f", addLegend = F) 
png(paste0(dirs$plot, "final/Figure03A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure03B:

p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)
png(paste0(dirs$plot, "final/Figure03B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure03C:

formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Demuth", n = nCoef, type = "continuous"); colName <- "Demuth"
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol)
png(paste0(dirs$plot, "final/Figure03C.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure03D:

p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
png(paste0(dirs$plot, "final/Figure03D.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
## Figure03E:

p <- custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "RTcleaned_n", zVar = "stakes_f", subVar = "subject_f", yLim = yLimRT)
png(paste0(dirs$plot, "final/Figure03E.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ------------------------------------------------------- #
### Figure03F:

formula <- "RTcleaned_z ~ reqCongruency_f * stakes_f + (reqCongruency_f * stakes_f|subject_f)"
mod <- fit_lmem(formula)
selCol <- met.brewer("Isfahan1", n = 8, type = "continuous"); selCol <- selCol[c(2, 5, 7)]; colName <- "Isfahan1"
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol)
png(paste0(dirs$plot, "final/Figure03F.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# END OF FILE.
