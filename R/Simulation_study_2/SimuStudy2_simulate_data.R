
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(MASS)
library(INLA)
library(rgdal)
library(splines)


# Define the Scenario and the Subscenario
Scenario <- 1
# Scenario <- 2
# Scenario <- 3

Subscenario <- 80
# Subscenario <- 50
# Subscenario <- 20



##################################################
# Load the data simulated for Simulation study 1 #
##################################################
load(paste0("../../Simulated_data/Simu1_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))
N <- nrow(Data)



#######################
# Simulate the counts #
#######################
if(Scenario==1){
  log.risk <- 0.2*Data$X1 + 0*Data$X2
}
if(Scenario==2){
  log.risk <- 0.2*Data$X1 + 0*Data$X2 + Data$spat
}
if(Scenario==3){
  log.risk <- 0.2*Data$X1 + 0*Data$X2 + Data$spat
}

lambda <- Data$E*exp(log.risk)


n.sim <- 100
simu.O <- vector("list", n.sim)

for(i in 1:n.sim) {
  set.seed(10+i)
  O <- rpois(N, lambda)
  simu.O[[i]] <- O
}



###########################
# Save the simulated data #
###########################
if(!file.exists("Simulated_data")) {dir.create("../../Simulated_data")}
save(simu.O, log.risk, Q.xi, Data, carto,
     file=paste0("../../Simulated_data/Simu2_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))


