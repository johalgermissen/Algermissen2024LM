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

## Load custom functions:
source(paste0(dirs$codeDir, "functions/00_mgngstakes_functions_regression.R")) # Load functions

## Load packages:
library(ggplot2)
library(ggthemes)

# ============================================================================ #
#### Define data: ####

hypData <- data.frame(condition_n <- 1:4)
hypData$congruency_f <- factor(c("cong", "cong", "incong", "incong"))
hypData$stakes_f <- factor(c("high", "low", "high", "low"))
hypData$Inv <- c(0.90, 0.85, 0.60, 0.65)
hypData$EVC <- c(0.85, 0.85, 0.70, 0.65)
hypData$EVCRT <- c(0.60, 0.60, 0.80, 0.72)

# ============================================================================ #
#### Plot: ####

## Select color:
selCol <- c("orange", "purple")
selCol <- met.brewer("Troy", n = 8, type = "continuous"); selCol <- selCol[c(2, 4)]; colName <- "Troy" # dark blue, blue, light blue
selCol <- met.brewer("Cassatt1", n = 8, type = "continuous"); selCol <- selCol[c(1, 3)]; colName <- "Cassatt1" # dark blue, blue, light blue
selCol <- met.brewer("VanGogh1", n = 8, type = "continuous"); selCol <- selCol[c(2, 4)]; colName <- "VanGogh1" # dark blue, blue, light blue
selCol <- met.brewer("Cassatt2", n = 10, type = "continuous"); selCol <- selCol[c(2, 5)]; colName <- "Cassatt2" # dark blue, blue, light blue

## Display:
library(scales)
show_col(selCol)

## Plotting settings:
LWD <- retrieve_plot_defaults("LWD") # 1.3
FTS <- retrieve_plot_defaults("FTS")
dodgeVal <- 0.6
colAlpha <- 1
yLim <- c(0.5, 1)

## Select variable:
hypData$y <- hypData$Inv; yLab <- "Accuracy"; plotName <- "BiasInvigoration"
hypData$y <- hypData$EVC; yLab <- "Accuracy"; plotName <- "MotivationForControl"
hypData$y <- hypData$EVCRT; yLab <- "RT (in sec.)"; plotName <- "MotivationForControl_RT"

## Plot:
p <- ggplot(hypData, aes(x = congruency_f, fill = stakes_f, y = y)) + 
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = dodgeVal,
               lwd = LWD, color = "black") + 
  scale_fill_manual(values = selCol, limits = levels(hypData$stakes_f)) + 
  labs(x = "Congruency", y = yLab, fill = "Stakes") +
  coord_cartesian(ylim = yLim) + 
  theme_classic() + 
  theme(axis.text = element_text(size = FTS),
        axis.title = element_text(size = FTS), 
        plot.title = element_text(size = FTS, hjust = 0.5), 
        legend.title = element_blank(), legend.position = "none",
        axis.line = element_line(colour = 'black')) # , linewidth = LWD)) # fixed font sizes
print(p)  

## Save:
plotName <- paste0("hypothesis_", yName, "_", colName)
png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)
print(p)
dev.off()

# END OF FILE.
