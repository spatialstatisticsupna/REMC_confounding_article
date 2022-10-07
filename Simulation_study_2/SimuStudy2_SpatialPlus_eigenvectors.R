
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


Data$ID.area <- seq(1, 70, 1)
N <- dim(Data)[1]




###############################
# Hyperpriors for INLA models #
###############################

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

compute.patterns <- FALSE  ## Set compute.patterns=FALSE if posterior patterns are not required
strategy <- "laplace"




########################################################
# Eigenvectors and eigenvalues of the adjacency matrix #
########################################################

eigendecomp <- eigen(Q.xi)
eigendecomp$values


# Eigenvectors corresponding to 20 lowest non-null eigenvalues
eigen.vect <- eigendecomp$vectors[, (N-20):(N-1)]




##############################################
# Load spatial models to compute the weights #
##############################################

load(paste0("Results_SimuStudy2/Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))




#################################
# COVARIATE MODEL: eigenvectors #
#################################

n.sim <- 100

# SpatPlus5
f.X1.Cov5 <- X1.weights ~ -1 + W.intercept + V1+ V2 + V3 + V4 + V5
f.X2.Cov5 <- X2.weights ~ -1 + W.intercept + V1+ V2 + V3 + V4 + V5

# SpatPlus10
f.X1.Cov10 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10
f.X2.Cov10 <- X2.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10

# SpatPlus15
f.X1.Cov15 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        V11 + V12 + V13 + V14 + V15
f.X2.Cov15 <- X2.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                           V11 + V12 + V13 + V14 + V15

# SpatPlus20
f.X1.Cov20 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20
f.X2.Cov20 <- X2.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20


## Choose the covariate model ##
f.X1.Cov <- f.X1.Cov5
f.X2.Cov <- f.X2.Cov5


X1.Res <- vector("list", n.sim)
X2.Res <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))

  # Compute the weights
  W <- Spat[[i]]$summary.fitted.values$mode*Data$expected
  W.sqrt <- diag(sqrt(W))
  W.sqrt.inv <- diag(1/sqrt(W))
  
  # Compute the residuals
  W.eigen.vect <- W.sqrt%*%eigen.vect
  Data2 <- cbind(Data, W.eigen.vect)
  colnames(Data2) <- c(colnames(Data), paste0("V", 20:1))
  Data2$X1.weights <- W.sqrt%*%Data$X1
  Data2$X2.weights <- W.sqrt%*%Data$X2
  Data2$W.intercept <- W.sqrt%*%rep(1, N)
  
  Covariate.X1 <- inla(f.X1.Cov, family = "gaussian", data=Data2,
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                       control.predictor=list(compute=TRUE),
                       control.inla=list(strategy=strategy))
  Covariate.X2 <- inla(f.X2.Cov, family = "gaussian", data=Data2,
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                       control.predictor=list(compute=TRUE),
                       control.inla=list(strategy=strategy))

  
  X1.Res.eigen <- W.sqrt.inv%*%(Data2$X1.weights - Covariate.X1$summary.fitted.values[, 1])
  X1.Res[[i]] <- as.vector(scale(X1.Res.eigen))
  
  X2.Res.eigen <- W.sqrt.inv%*%(Data2$X2.weights - Covariate.X2$summary.fitted.values[, 1])
  X2.Res[[i]] <- as.vector(scale(X2.Res.eigen))
  
  Data2 <- NULL

}




########################
# SPATIAL+ MODEL: ICAR #
########################

f.SpatPlus <- O ~ 1 + X1.Res + X2.Res + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                          scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 



SpatPlus <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))
  # Simulated counts
  Data$O <- simu.O[[i]]
  Data$X1.Res <- X1.Res[[i]]
  Data$X2.Res <- X2.Res[[i]]

  # Fit the model
  SpatPlus[[i]] <- inla(f.SpatPlus, family="poisson", data=Data, E=expected,
                     control.predictor=list(compute=TRUE, cdf=c(log(1))),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     control.inla=list(strategy=strategy))
}




####################
# Save the results #
####################

SpatPlus5 <- SpatPlus
save(SpatPlus5, file=paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlus5_cor", Subscenario, ".Rdata"))


# SpatPlus10 <- SpatPlus
# save(SpatPlus10, file=paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlus10_cor", Subscenario, ".Rdata"))

# SpatPlus15 <- SpatPlus
# save(SpatPlus15, file=paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlus15_cor", Subscenario, ".Rdata"))

# SpatPlus20 <- SpatPlus
# save(SpatPlus20, file=paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlus20_cor", Subscenario, ".Rdata"))



