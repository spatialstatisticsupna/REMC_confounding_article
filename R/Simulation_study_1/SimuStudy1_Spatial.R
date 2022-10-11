
rm(list=ls())
setwd("")


# Load packages
library(INLA)




#######################
# Load simulated data #
#######################

## Define the Scenario and Subscenario ##

Scenario <- 1  # 2,3
Subscenario <- 80  # 50,20


load(paste0("Simulated_data/Simu1_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))


Data$ID.area <- seq(1, 70, 1)




###############################
# Hyperpriors for INLA models #
###############################

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

compute.patterns <- FALSE  ## Set compute.patterns=FALSE if posterior patterns are not required
strategy <- "laplace"




#################
# SPATIAL MODEL #
#################

n.sim <- 100


f.Spat <- O ~ X1 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                     hyper=list(prec=list(prior=sdunif))) 


Spat <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))
  # Simulated counts
  Data$O <- simu.O[[i]]
  # Fit the model
  Spat[[i]] <- inla(f.Spat, family="poisson", data=Data, E=expected,
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy=strategy))
}




####################
# Save the results #
####################

save(Spat, file=paste0("Results_SimuStudy1/Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))



