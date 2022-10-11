
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


load(paste0("Simulated_data/Simu2_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))




###############################
# Hyperpriors for INLA models #
###############################

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

compute.patterns <- FALSE  ## Set compute.patterns=FALSE if posterior patterns are not required
strategy <- "laplace"




##############
# NULL MODEL #
##############

n.sim <- 100


f.Null <- O ~ 1 + X1 + X2


Null <- vector("list", n.sim)

for (i in 1:n.sim){
  print(paste0("i=", i))
  # Simulated counts
  Data$O <- simu.O[[i]]
  # Fit the model
  Null[[i]] <- inla(f.Null, family="poisson", data=Data, E=expected,
                        control.predictor=list(compute=TRUE, cdf=c(log(1))),
                        control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                        control.inla=list(strategy=strategy))
}



####################
# Save the results #
####################

save(Null, file=paste0("Results_SimuStudy2/Scenario", Scenario, "_Null_cor", Subscenario, ".Rdata"))


