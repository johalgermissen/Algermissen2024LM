#### All packages and settings for analysis ####

# =========================================================================================== #
#### Packages: ####

cat("Start loading packages\n")

# ----------------------------------------- #
## Descriptives and plots:

require(car)
require(Rmisc) # for summarySEwithin
require(ggplot2)
require(lattice)
require(ltm)
require(pastecs) # for stat.desc
require(psych)
require(effsize) # for Cohen.d

require(stringr) # for str_pad etc.

# ----------------------------------------- #
## For reading Matlab files:

require(rmatio) # read.mat
require(data.table) # rbindlist

# ----------------------------------------- #
## RM-ANOVA:

require(ez)

# ----------------------------------------- #
## Tidyverse:

require(plyr)
require(dplyr)
require(magrittr)

# ----------------------------------------- #
## Linear mixed effects models:

require(lme4)
require(afex)
require(effects)
require(emmeans)
require(DescTools)

# ----------------------------------------- #
## Generalized additive mixed models:

require(mgcv)
require(itsadug)

# ----------------------------------------- #
## For corrplots:

library(corrplot) # for corrplots
library(synthesisr) # for line breaks in title

# ----------------------------------------- #
## Color bars:

require(sommer)
require(viridis)
require(RColorBrewer)
require(MetBrewer)
require(scales)

# ----------------------------------------- #
## For raincloud plots:

require(readr)
require(tidyr)
require(ggplot2)
require(Hmisc)
require(plyr)
require(RColorBrewer)
require(reshape2)
require(ggstatsplot) # for geom_flat_violin
require(gghalves) # for half plots
require(ggbeeswarm) # for ggbeeswarm
require(ggthemes)

# ----------------------------------------- #
## Facilitate detecting when model finished:

require(beepr)

# ============================================================================ #
#### General settings: ####

cat("Set scipen to 20\n")
options(scipen = 20)

cat("Set contrasts to sum-to-zero coding\n")
options(contrasts = c("contr.sum", "contr.poly"))

# ============================================================================ #
#### Set seed: ####

mySeed <- 70
cat(paste0("Set seed to ", mySeed, "\n"))
set.seed(mySeed)

# ============================================================================ #
#### sessionInfo: ####

# > sessionInfo()
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /opt/R/4.1.0/lib64/R/lib/libRblas.so
# LAPACK: /opt/R/4.1.0/lib64/R/lib/libRlapack.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] beepr_1.3          ggthemes_4.2.4     ggbeeswarm_0.7.2   gghalves_0.1.4     ggstatsplot_0.9.0  reshape2_1.4.4     Hmisc_5.1-0       
# [8] tidyr_1.3.0        readr_2.0.2        scales_1.2.1       MetBrewer_0.2.0    RColorBrewer_1.1-3 viridis_0.6.4      viridisLite_0.4.2 
# [15] sommer_4.3.2       crayon_1.5.2       synthesisr_0.3.0   corrplot_0.92      itsadug_2.4.1      plotfunctions_1.4  mgcv_1.9-0        
# [22] nlme_3.1-152       DescTools_0.99.49  emmeans_1.8.7      effects_4.2-2      afex_1.3-0         lme4_1.1-34        Matrix_1.6-0      
# [29] magrittr_2.0.3     dplyr_1.1.2        ez_4.4-0           data.table_1.14.8  rmatio_0.18.0      stringr_1.5.0      effsize_0.8.1     
# [36] psych_2.3.6        pastecs_1.3.21     ltm_1.2-0          polycor_0.8-1      msm_1.6.9          MASS_7.3-54        ggplot2_3.4.2     
# [43] Rmisc_1.5.1        plyr_1.8.8         lattice_0.21-8     car_3.1-2          carData_3.0-5     
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1           backports_1.4.1        splines_4.1.0          gmp_0.7-2              kSamples_1.2-9         TH.data_1.1-2         
# [7] SuppDists_1.1-9.7      digest_0.6.33          htmltools_0.5.6        lmerTest_3.1-3         fansi_1.0.4            memoise_2.0.1         
# [13] checkmate_2.2.0        paletteer_1.4.0        cluster_2.1.2          tzdb_0.4.0             sandwich_3.0-1         colorspace_2.1-0      
# [19] mitools_2.4            xfun_0.40              Exact_3.2              zeallot_0.1.0          survival_3.2-11        zoo_1.8-12            
# [25] glue_1.6.2             gtable_0.3.0           statsExpressions_1.2.0 Rmpfr_0.9-3            abind_1.4-5            mvtnorm_1.2-2         
# [31] DBI_1.1.3              PMCMRplus_1.9.2        Rcpp_1.0.11            performance_0.10.4     xtable_1.8-4           htmlTable_2.3.0       
# [37] foreign_0.8-81         proxy_0.4-27           Formula_1.2-5          survey_4.1-1           datawizard_0.8.0       htmlwidgets_1.5.4     
# [43] httr_1.4.6             ellipsis_0.3.2         reshape_0.8.9          pkgconfig_2.0.3        multcompView_0.1-9     nnet_7.3-16           
# [49] utf8_1.2.3             tidyselect_1.2.0       rlang_1.1.1            cachem_1.0.8           munsell_0.5.0          cellranger_1.1.0      
# [55] tools_4.1.0            cli_3.6.1              audio_0.1-10           generics_0.1.3         evaluate_0.21          fastmap_1.1.1         
# [61] BWStest_0.2.2          rematch2_2.1.2         knitr_1.36             admisc_0.33            purrr_1.0.2            rootSolve_1.8.2.3     
# [67] WRS2_1.1-3             correlation_0.7.1      compiler_4.1.0         rstudioapi_0.15.0      beeswarm_0.4.0         e1071_1.7-13          
# [73] tibble_3.2.1           stringi_1.7.12         parameters_0.21.1      nloptr_1.2.2.2         vctrs_0.6.3            stringdist_0.9.10     
# [79] mc2d_0.1-21            pillar_1.9.0           lifecycle_1.0.3        estimability_1.4.1     insight_0.19.3         lmom_2.9              
# [85] patchwork_1.1.1        R6_2.5.1               gridExtra_2.3          vipor_0.4.5            gld_2.6.6              codetools_0.2-18      
# [91] boot_1.3-28            withr_2.5.0            mnormt_2.1.1           multcomp_1.4-17        bayestestR_0.13.1      expm_0.999-6          
# [97] parallel_4.1.0         hms_1.1.1              grid_4.1.0             rpart_4.1-15           coda_0.19-4            class_7.3-19          
# [103] minqa_1.2.5            rmarkdown_2.11         numDeriv_2016.8-1.1    base64enc_0.1-3  

# > loadedNamespaces()
# [1] "readxl"           "backports"        "Hmisc"            "corrplot"         "plyr"             "ltm"              "splines"         
# [8] "polycor"          "gmp"              "kSamples"         "ggplot2"          "TH.data"          "SuppDists"        "digest"          
# [15] "htmltools"        "viridis"          "lmerTest"         "fansi"            "memoise"          "magrittr"         "checkmate"       
# [22] "paletteer"        "itsadug"          "cluster"          "tzdb"             "readr"            "methods"          "sandwich"        
# [29] "beepr"            "colorspace"       "mitools"          "xfun"             "dplyr"            "crayon"           "Exact"           
# [36] "lme4"             "zeallot"          "survival"         "zoo"              "glue"             "utils"            "gtable"          
# [43] "emmeans"          "statsExpressions" "car"              "Rmpfr"            "abind"            "scales"           "pastecs"         
# [50] "mvtnorm"          "DBI"              "graphics"         "ggthemes"         "PMCMRplus"        "synthesisr"       "Rcpp"            
# [57] "performance"      "viridisLite"      "xtable"           "htmlTable"        "foreign"          "proxy"            "effsize"         
# [64] "Formula"          "survey"           "MetBrewer"        "base"             "datawizard"       "htmlwidgets"      "httr"            
# [71] "RColorBrewer"     "ellipsis"         "reshape"          "pkgconfig"        "multcompView"     "nnet"             "utf8"            
# [78] "ez"               "tidyselect"       "rlang"            "reshape2"         "stats"            "cachem"           "munsell"         
# [85] "cellranger"       "tools"            "cli"              "audio"            "generics"         "evaluate"         "stringr"         
# [92] "fastmap"          "BWStest"          "rematch2"         "grDevices"        "knitr"            "sommer"           "admisc"          
# [99] "purrr"            "rootSolve"        "WRS2"             "rmatio"           "nlme"             "correlation"      "compiler"        
# [106] "rstudioapi"       "beeswarm"         "e1071"            "plotfunctions"    "tibble"           "afex"             "DescTools"       
# [113] "stringi"          "parameters"       "lattice"          "Matrix"           "psych"            "nloptr"           "vctrs"           
# [120] "effects"          "stringdist"       "msm"              "mc2d"             "pillar"           "lifecycle"        "estimability"    
# [127] "data.table"       "insight"          "lmom"             "patchwork"        "R6"               "Rmisc"            "gridExtra"       
# [134] "vipor"            "gld"              "codetools"        "boot"             "MASS"             "withr"            "mnormt"          
# [141] "multcomp"         "datasets"         "mgcv"             "bayestestR"       "expm"             "parallel"         "hms"             
# [148] "gghalves"         "grid"             "rpart"            "tidyr"            "coda"             "class"            "minqa"           
# [155] "rmarkdown"        "carData"          "numDeriv"         "base64enc"        "ggbeeswarm"       "ggstatsplot"  

# END
