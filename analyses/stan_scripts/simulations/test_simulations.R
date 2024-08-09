## test_simulations.R

# ------------------------------------------------------------------------- #
### Quick test version for debugging:

## Select model:
iMod <- 12
iMod_str <- str_pad(iMod, width = 2, side = "left", pad = "0")
cat(paste0("Load model ", iMod_str, "\n"))

## Load parameters:
load(paste0(dirs$results, "M", iMod_str, suffix, "_fitted.Rdata"))
fit <- f2
nSample <- dim(fit)[1]; nChain <- dim(fit)[2]; nParam <- dim(fit)[3]
posterior <- rstan::extract(fit, inc_warmup = FALSE, permuted = TRUE)

## Reshape into big matrix:
allParNames <- names(posterior) # names of all parameters
subParNames <- allParNames[grepl("sub_", allParNames, fixed = T)] # names of subject level parameters
nParamSub <- length(subParNames) # number subject level parameters
allSubParMat <- array(NA, c(nSub, nParamSub, nSample*nChain)) # subjects, parameters, iterations
for (iParam in 1:nParamSub){ # iParam <- 1 # for each parameter
  allSubParMat[, iParam, ] <- t(posterior[[subParNames[iParam]]]) # extract lists of all subjects for this parameter, transpose, save in matrix
}

## Select parameters:
iIter <- 1
par <- allSubParMat[, , selIter[iIter]]
colnames(par) <- subParNames # add column ways
cat(paste0("Simulate based on parameters ", paste0(colnames(par), collapse = ", "), "\n"))

# ------------------------------------------- #
## Start simulations:

source(paste0(dirs$simscripts, "RLDDM_modSims", suffix, ".R")) # model code for simulations
run_simulations <- eval(parse(text = paste0(simType, "_M", iMod_str)))
cat(paste0("Start simulations for M", iMod_str, " .... \n"))
tmp <- run_simulations(data, par) # takes 15 sec. for 35 people
cat(paste0("Finished simulations for M", iMod_str, "! :-)\n"))

# ------------------------------------------- #
### Inspect:

# head(tmp)

## Responses and accuracy:
table(tmp$sim_resp)
mean(tmp$sim_resp == 1)
mean(tmp$sim_resp == tmp$reqAction)
mean(tmp$sim_resp == tmp$resp)
mean(tmp$sim_resp[tmp$valence == 1]) - mean(tmp$sim_resp[tmp$valence == 0])
mean(tmp$sim_resp[tmp$stakes == 1]) - mean(tmp$sim_resp[tmp$stakes == 0])

## RTs:
densityplot(tmp$sim_rt)
mean(tmp$sim_rt > 1.3)
mean(tmp$sim_rt[tmp$valence == 0]) - mean(tmp$sim_rt[tmp$valence == 1])
mean(tmp$sim_rt[tmp$stakes == 1]) - mean(tmp$sim_rt[tmp$stakes == 0])

## Outcomes:
table(tmp$sim_outcome)

## Q-values:
if (all(is.na(tmp$QGo))){warning("QGo is all NA")}
densityplot(tmp$QGo)
densityplot(tmp$QNoGo)
# densityplot(tmp$QGo - tmp$QNoGo)
# plot(tmp$QGo - tmp$QNoGo)

## Separately for Go/NoGo cues:
selIdx <- which(tmp$reqAction == 1) 
densityplot(tmp$QGo[selIdx] - tmp$QNoGo[selIdx])
plot(tmp$QGo[selIdx] - tmp$QNoGo[selIdx])
mean((tmp$QGo[selIdx] - tmp$QNoGo[selIdx]) < 0)
stopifnot(mean((tmp$QGo[selIdx] - tmp$QNoGo[selIdx]) < 0) < 0.04)

selIdx <- which(tmp$reqAction == 0) 
densityplot(tmp$QGo[selIdx] - tmp$QNoGo[selIdx])
plot(tmp$QGo[selIdx] - tmp$QNoGo[selIdx])
mean((tmp$QGo[selIdx] - tmp$QNoGo[selIdx]) > 0)
stopifnot(mean((tmp$QGo[selIdx] - tmp$QNoGo[selIdx]) > 0) < 0.04)

## Separately for Win/Avoid cues:
selIdx <- which(tmp$valence == 1) 
densityplot(tmp$QGo[selIdx])
densityplot(tmp$QNoGo[selIdx])
sum(tmp$QNoGo[selIdx] < 0)
stopifnot(all(tmp$QNoGo[selIdx] > 0))

selIdx <- which(tmp$valence == 0) 
densityplot(tmp$QGo[selIdx])
densityplot(tmp$QNoGo[selIdx])
sum(tmp$QNoGo[selIdx] > 0)
stopifnot(all(tmp$QNoGo[selIdx] < 0))

# ------------------------------------------- #
## Selected subject:
iSub <- 1
selIdx <- which(tmp$subject == iSub) # select subjects
selIdx <- which(tmp$subject == iSub & tmp$reqAction == 1) # select one required action
selIdx <- which(tmp$subject == iSub & tmp$valence == 0) # select valence

densityplot(tmp$QGo[selIdx])
densityplot(tmp$QNoGo[selIdx])
densityplot(tmp$QGo[selIdx] - tmp$QNoGo[selIdx])
plot(tmp$QGo[selIdx] - tmp$QNoGo[selIdx])

# END OF FILE.