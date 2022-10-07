
rm(list=ls())
setwd("")


# Load packages
library(TMGMRF)



#######################
# Load simulated data #
#######################

## Define the Scenario and Subscenario ##

Scenario <- 1  # 2,3
Subscenario <- 80  # 50,20


load(paste0("Simulated_data/Simu1_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))


Data$ID.area <- seq(1, 70, 1)
N <- dim(Data)[1]




#######################
# Configuration TGMRF #
#######################

family = "poisson"
type_data = "gamma-scale"
seed <- 1234


formula <- y ~ X1

nburnin <- 2000
thin <- 20
nsamp <- 10000




###############
# TGMRF MODEL #
###############

n.sim <- 100


Model.GSC <- vector("list", n.sim)

for (i in 1:n.sim){
  print(paste0("i=", i))
  # Simulated counts
  Data$O <- simu.O[[i]]
  df <- data.frame(y = Data$O, X1 = Data$X1)
  # Fit the model
  out <- tgmrf(
    formula = formula, data = df, E = Data$expected, n=N,
    W = Q.xi.TGMRF,
    family = family, model = type_data,
    nchains = 3,
    nsamp = nsamp, nburnin = nburnin, thin = thin,
    seed = seed
    )
  Model.GSC[[i]] <- out
}




####################
# Save the results #
####################

save(Model.GSC, file=paste0("Results_SimuStudy1/Scenario", Scenario, "_TGMRF_GSC_cor", Subscenario, ".Rdata"))

