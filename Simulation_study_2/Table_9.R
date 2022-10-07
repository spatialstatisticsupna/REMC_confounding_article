
rm(list=ls())
setwd("")


n.sim <- 100 


# Define the Scenario (1,2 or 3) and the Subscenario (cor=0.80,0.50 or 0.20)
Scenario <- 1  # 2,3
Subscenario <- 80  # 50,20




##########################################################
# Load models and save all the models in the same .Rdata #
##########################################################

load(paste0("Results_SimuStudy2/Scenario", Scenario, "_Null_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_RSR_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_summary_TGMRF_GSC_cor", Subscenario, ".Rdata"))

load(paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlus5_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlus10_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlus15_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlus20_cor", Subscenario, ".Rdata"))

load(paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlusP1_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlusTP1_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlusP2_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlusTP2_cor", Subscenario, ".Rdata"))


save.image(paste0("Results_SimuStudy2/Scenario", Scenario, "_all_models_cor", Subscenario, ".Rdata"))




###################################
# TABLE 9: type-S error for beta2 #
###################################

TypeS.error.inla <- function(Model){
  CI.L <- Model$summary.fixed$'0.025quant'[3]
  CI.U <- Model$summary.fixed$'0.975quant'[3]
  ifelse(0>=CI.L & 0<=CI.U, 0, 1)
}


TypeS.error.nimble <- function(Model){
  CI.L=Model$quantiles[grep('X2', rownames(Model$quantiles), value=TRUE), 1]                  
  CI.U=Model$quantiles[grep('X2', rownames(Model$quantiles), value=TRUE), 5]
  ifelse(0>=CI.L & 0<=CI.U, 0, 1)
}


TypeS.error.beta2 <- data.frame(Null=mean(unlist(lapply(Null, function(x) TypeS.error.inla(x))))*100,
                                Spat=mean(unlist(lapply(Spat, function(x) TypeS.error.inla(x))))*100,
                                RSR=mean(unlist(lapply(RSR, function(x) TypeS.error.inla(x))))*100,
                                TGMRF1=mean(unlist(lapply(TGMRF.GSC, function(x) TypeS.error.nimble(x))))*100,
                                SpatPlus5=mean(unlist(lapply(SpatPlus5, function(x) TypeS.error.inla(x))))*100,
                                SpatPlus10=mean(unlist(lapply(SpatPlus10, function(x) TypeS.error.inla(x))))*100,
                                SpatPlus15=mean(unlist(lapply(SpatPlus15, function(x) TypeS.error.inla(x))))*100,
                                SpatPlus20=mean(unlist(lapply(SpatPlus20, function(x) TypeS.error.inla(x))))*100,
                                SpatPlusP1=mean(unlist(lapply(SpatPlusP1, function(x) TypeS.error.inla(x))))*100,
                                SpatPlusTP1=mean(unlist(lapply(SpatPlusTP1, function(x) TypeS.error.inla(x))))*100,
                                SpatPlusP2=mean(unlist(lapply(SpatPlusP2, function(x) TypeS.error.inla(x))))*100,
                                SpatPlusTP2=mean(unlist(lapply(SpatPlusTP2, function(x) TypeS.error.inla(x))))*100)   



TypeS.error.beta2 <- as.data.frame(t(TypeS.error.beta2))
colnames(TypeS.error.beta2) <- paste0("cor0.", Subscenario)
print(TypeS.error.beta2)
