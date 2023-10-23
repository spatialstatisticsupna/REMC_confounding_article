
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(INLA)




#######################################
# Load and prepare the simulated data #
#######################################

## Define the Scenario and Subscenario ##

Scenario <- 1
# Scenario <- 2
# Scenario <- 3

Subscenario <- 80
# Subscenario <- 50
# Subscenario <- 20


load(paste0("../../Simulated_data/Simu2_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))


Data$ID.area <- seq(1, 70, 1)
N <- nrow(Data)

ones.S <- rep(1, N)
p <- 2 # Number of covariates
Beta.df <- as.matrix(Data[,paste("X",1:p,sep="")])




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

load(paste0("Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))




#############
# RSR MODEL #
#############

n.sim <- 100


f.RSR <- O ~ -1 + Intercept + beta1 + beta2 + f(ID.area, model="besag", graph=Q.xi, rankdef=1, 
                                                constr=TRUE, hyper=list(prec=list(prior=sdunif)))

RSR <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))
  
  # Simulated counts
  Data$O <- simu.O[[i]]
  # Compute the weights
  Mod.Spat <- Spat[[i]]
  W <- diag(Mod.Spat$summary.fitted.values$mode*Data$expected)
  W.sqrt <- diag(sqrt(diag(W)))

  
  X <- cbind(ones.S, as.matrix(Beta.df))
  P <- W.sqrt%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W.sqrt
  Pc <- diag(N)-P

  eigen.Pc <- eigen(Pc)
  L <- eigen.Pc$vectors[,eigen.Pc$values>1e-12]

  M <- solve(W.sqrt)%*%L%*%t(L)%*%W.sqrt
  Z.area <- M%*%diag(N)

  Data.restricted <- list(O=Data$O, E=Data$expected,
                          X1=Data$X1, X2=Data$X2,
                          Intercept=c(1,rep(NA,N+2)),
                          beta1=c(NA,1,rep(NA,N+1)),
                          beta2=c(NA,NA,1,rep(NA,N)),
                          ID.area=c(rep(NA,1+1+1),1:N))

  # RSR model
  Apredictor <- cbind(ones.S, Beta.df, Z.area)

  RSR[[i]] <- inla(f.RSR, family="poisson", data=Data.restricted, E=E,
                   control.predictor=list(compute=TRUE, A=Apredictor),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.inla=list(strategy=strategy),
                   control.mode=list(theta=Mod.Spat$mode$theta, restart=FALSE),
                   control.fixed=list(prec=0))

}




####################
# Save the results #
####################

save(RSR, file=paste0("Scenario", Scenario, "_RSR_cor", Subscenario, ".Rdata"))


