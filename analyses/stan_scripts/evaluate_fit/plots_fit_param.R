# ==================================================================================== #
#' plots_fit_param.R
#' 1) Loop over models and extract WAIC/ LOOIC, plot with bar plot
#' 2) Plot group-level parameter densities
#' 3) Plot violin plot of subject-level mean parameters

rm(list = ls())

# install.packages("devtools")
# library(devtools)
# install_github("easyGgplot2", "kassambara")

# ==================================================================================== #
#### Load packages: ####

require(ggplot2)
# library(easyGgplot2) # http://www.sthda.com/english/wiki/ggplot2-violin-plot-easy-function-for-data-visualization-using-ggplot2-and-r-software
# require(ggpattern)
require(latex2exp)

require(plyr)
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

source(paste0(dirs$helpers, "custom_plots_stan.R"))

# ======================================================================================================= #
#### Loop over models to extract WAIC and LOOIC: ####

## Settings:
nMod <- 12
# suffix <- "_standard" # 
# suffix <- "_deltaIntSlope" #
suffix <- "_sepParam" #

## Initialize data frame:
evidence <- data.frame(iMod = seq(1, nMod),
                       modName = paste0("M", as.character(seq(1, nMod))), 
                       WAIC = rep(NA, nMod),
                       LOOIC = rep(NA, nMod))

## Select actual models

# allFiles <- list.files(dirs$results, pattern=paste0(suffix,"_fitted"),full=TRUE) # check existing files

modVec <- seq(1, nMod)

# modVec <- seq(1,8)
# modVec <- modVec[which(!(modVec %in% c(13,16,17)))] # exclude models bad models for full
# modVec <- modVec[which(!(modVec %in% c(15,16,18)))] # exclude models not fitted in _fixTresh

for (model2fit in modVec){ # model2fit <- 1
    model2fit_str <- str_pad(model2fit, width = 2, side = "left", pad = "0")
    cat(paste0("Load model ", model2fit_str, "\n"))
    stanFitsFile <- paste0(dirs$results, "M", model2fit_str, suffix, "_fitted.Rdata") # model results: f1 and f2
    load(stanFitsFile)
    # ---------------------------------------------------------------------- #
    ## WAIC:
    # cat(paste0(suffix, ": M", model2fit_str,": WAIC = ",round(loo::waic(extract_log_lik(f2))$waic,1),"\n"))
    evidence$WAIC[model2fit] <- round(loo::waic(extract_log_lik(f2))$waic, 1)
    # ---------------------------------------------------------------------- #
    ## LOOIC:
    # cat(paste0(suffix, ": M", model2fit_str,": LOOIC = ",round(loo(f2)$looic,1),"\n"))
    evidence$LOOIC[model2fit] <- round(loo(f2)$looic, 1)
  }
# nMod <- nrow(modVec) # update number of retrieved models
# evidence # check

# ======================================================================================================= #
#### Select models to plot: ####

## Copy over all models:
selEvidence <- evidence

## Exclude non-convergent models:
# selEvidence <- evidence[which(evidence$WAIC<3000),] # exclude non-converging models

# selEvidence <- evidence[1:8,] # select DDM, RL-DDM, stakes on bias, dwell time on drift

## For poster:
# selEvidence <- evidence[c(1,2,3,6,7),] # select DDM, RL-DDM, stakes on bias, dwell time on drift


# selEvidence <- evidence[c(1,2,4,6,12,13),] # select DDM, RL-DDM, stakes on bias, dwell time on drift

## Bias only:
# selEvidence <- evidence[c(1,2,3,5,18,19),] # select DDM, RL-DDM, stakes on bias, dwell time on drift

## Drift rate only:
# selEvidence <- evidence[c(1,2,4,6,8,9,10,11),] # select DDM, RL-DDM, stakes on bias, dwell time on drift

## Difference terms:
# selEvidence <- evidence[c(1,2,12,13,17,18,19),] # select DDM, RL-DDM, effect of differences
# selEvidence <- evidence[c(1,2,12,13,17),] # select DDM, RL-DDM, effect of differences on drift

# selEvidence <- evidence[c(1,2,4,6,8),] # select DDM, RL-DDM, effect on drift
# selEvidence <- evidence[c(1,2,8,9,10,11),] # select DDM, RL-DDM, single effects on drift
# selEvidence <- evidence[c(1,2,12,13,17),] # select DDM, RL-DDM, difference variables
# selEvidence <- evidence[c(1,2,7,8,11,17),] # select DDM, RL-DDM, "maximal models"
# selEvidence <- evidence[c(1,2,3,4,5,6,15),] # ,16),] # select DDM, RL-DDM, effect on single parameter

# ======================================================================================================= #
#### Define output name: ####

allModNames <- paste(selEvidence$modName, collapse = "")

## Rename models:
# selEvidence$modName <- paste0("M", as.character(seq(1,nrow(selEvidence))))

# ======================================================================================================= #
#### Delete outliers: ####

# maxVal <- 10000
# 
# ## Set non-converging models to negative difference value:
# selEvidence[which(selEvidence$WAIC>maxVal),"WAIC"] <- selEvidence$WAIC[1] + NA # set non-converging models to negative values
# selEvidence[which(selEvidence$LOOIC>maxVal),"LOOIC"] <- selEvidence$LOOIC[1] + NA # set non-converging models to negative values
# selEvidence

# ======================================================================================================= #
#### Select model metric: ####

## Option A: Select WAIC:
selEvidence$IC <- selEvidence$WAIC; ICname <- "WAIC"

## Option B: Select LOOIC:
selEvidence$IC <- selEvidence$LOOIC; ICname <- "LOOIC"

## Compute difference to minimal metric:
# selEvidence$difIC <- selEvidence$IC-max(selEvidence$IC,na.rm=T) # subtract minimum
selEvidence$difIC <- selEvidence$IC-selEvidence$IC[1] # subtract first model
# selEvidence$difIC[7] <- selEvidence$difIC[7] - 0.3
  
## Identify best model:
selEvidence$isMax <- factor(ifelse(selEvidence$difIC==min(selEvidence$difIC,na.rm=T),"best","other"))

# ======================================================================================================= #
#### Plot difference in IC to base model: ####

FTS <- 25

p <- ggplot(selEvidence, aes(x = reorder(modName, desc(iMod)), y = difIC, fill = isMax)) +
  geom_bar(stat="identity") + coord_flip()+ 
  scale_y_reverse(breaks = seq(0, -10000, -2000)) +
  scale_fill_manual(values=c("#FF8C00", "#FEE799")) +
  labs(title = "", y = paste0("Delta ", ICname), x = " ") +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.text = element_text(size = FTS), 
        axis.line = element_line(colour = "black"), # , linewidth = LWD), # fixed font sizes
        axis.title = element_text(size = FTS), 
        title = element_text(size = FTS+5), 
        legend.text = element_text(size = FTS))

print(p)
# ggsave(filename=paste(dirs$plots, ICname,"_", suffix, "_", paste(selEvidence$model, collapse = ""), ".png", sep = ""))
# ggsave(filename=paste0(dirs$plots, plotName, ".eps"), sep = "")

## Save:
plotName <- paste0(ICname, suffix, "_", allModNames)
png(paste0(dirs$plots, plotName, ".png"), width = 480*1.5, height = 480) # , type ="cairo")
print(p)
dev.off()

# ======================================================================================================= #
# ======================================================================================================= #
# ======================================================================================================= #
# ======================================================================================================= #
#### PARAMETERS: GROUP-LEVEL DENSITIES: ####

## Settings:
model2fit <- 12
# suffix <- "_standard"
# suffix <- "_deltaIntSlope"
suffix <- "_sepParam"

## Plot with ggplot:
source(paste0(dirs$helpers, "custom_plots_stan.R"))
plot_group_density(model2fit, suffix, pattern = "transf", isInteractive = F)
  
# ======================================================================================================= #
#### PARAMETERS: Violin plot SUBJECT-LEVEL MEANS: ####

## Settings:
# model2fit <- 12
# suffix <- "_standard"
# suffix <- "_deltaIntSlope"
# suffix <- "_sepParam"

## Plot:
source(paste0(dirs$helpers, "custom_plots_stan.R"))
plot_subject_means_density(model2fit, suffix, isInteractive = F)

# ======================================================================================================= #
#### PARAMETERS: Densities of all parameters (rasters with many plots): ####

## Settings:
model2fit <- 12
# suffix <- "_standard"
# suffix <- "_deltaIntSlope"

# source(paste0(dirs$helpers, "custom_plots_stan.R"))
# plot_stan_dens(model2fit = model2fit, suffix = suffix)

cat("Finished all plots :-) \n")

# ======================================================================================================= #
#### Difference between two parameters: ####

names(posterior)[grepl("transf", names(posterior))]

## var1 minus var2:
var1 <- "transf_mu_deltaWin"; var2 <- "transf_mu_deltaAvoid"
var1 <- "transf_mu_pi"; var2 <- "transf_mu_tau"
var1 <- "transf_mu_theta"; var2 <- "transf_mu_tau"
var1 <- "transf_mu_theta"; var2 <- "transf_mu_pi"

## Plot:
source(paste0(dirs$helpers, "custom_plots_stan.R"))
plot_density_diff(model2fit, suffix, var1, var2) # , LWD = 4, FTS = 70, WIDTH = 960, HEIGHT = 640)

# END OF FILE.