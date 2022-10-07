
rm(list=ls())
setwd("")


# Load packages
library(nimble)



######################
# Load TGMRF results #
######################

## Define the Scenario and Subscenario ##

Scenario <- 1  # 2,3
Subscenario <- 80  # 50,20


load(paste0("Results_SimuStudy2/Scenario", Scenario, "_TGMRF_GSC_cor", Subscenario, ".Rdata"))




#############################
# Summary and WAIC of TGMRF #
#############################

TGMRF.GSC <- vector("list", 0)
WAIC.GSC <- vector("list", 0)


n.sim <- 100


for (i in 1:n.sim){
  print(i)
  TGMRF.GSC[[i]] <- summary(Model.GSC[[i]]$samples)
  WAIC.GSC[[i]] <- Model.GSC[[i]]$WAIC
}



save(TGMRF.GSC, WAIC.GSC, 
     file=paste0("Results_SimuStudy2/Scenario", Scenario, "_summary_TGMRF_GSC_cor", Subscenario, ".Rdata"))
