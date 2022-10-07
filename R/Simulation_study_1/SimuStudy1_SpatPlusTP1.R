
rm(list=ls())
setwd("")


# Load packages
library(INLA)
library(sf)
library(mgcv)




#######################
# Load simulated data #
#######################

## Define the Scenario and Subscenario ##

Scenario <- 1  # 2,3
Subscenario <- 80  # 50,20


load(paste0("Simulated_data/Simu1_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))


Data$ID.area <- seq(1, 70, 1)
N <- dim(Data)[1]

coord <- coordinates(carto)
coord.df <-  as.data.frame(coord)
names(coord.df)<-c("x","y")




###############################
# Hyperpriors for INLA models #
###############################

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

compute.patterns <- FALSE  ## Set compute.patterns=FALSE if posterior patterns are not required
strategy <- "laplace"




##############################################
# Load spatial models to compute the weights #
##############################################

load(paste0("Results_SimuStudy1/Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))




###############################
# COVARIATE MODEL: TP splines #
###############################

n.sim <- 100

method<-"GCV.Cp" #"REML"


X1.Res <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))

  # Compute the weights
  W <- Spat[[i]]$summary.fitted.values$mode*Data$expected

  # Covariate model
  Covariate <- gam(Data$X1~1+s(x,y,k=17, fx=TRUE, m=2), weights=W, data=coord.df)$fitted.values
  res <- Data$X1-Covariate
  X1.Res[[i]] <- as.vector(scale(res))

}




########################
# SPATIAL+ MODEL: ICAR #
########################

f.SpatPlus <- O ~ 1 + X1.Res + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                 scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 


SpatPlusTP1 <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))

  Data$O <- simu.O[[i]]
  Data$X1.Res <- X1.Res[[i]]

  SpatPlusTP1[[i]] <- inla(f.SpatPlus, family="poisson", data=Data, E=expected,
                           control.predictor=list(compute=TRUE, cdf=c(log(1))),
                           control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                           control.inla=list(strategy=strategy))
}




####################
# Save the results #
####################

save(SpatPlusTP1, file=paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlusTP1_cor", Subscenario, ".Rdata"))

