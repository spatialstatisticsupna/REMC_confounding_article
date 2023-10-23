
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(INLA)


# Define the number of eigenvectors for SpatPlus
model <- "SpatPlus5"
# model <- "SpatPlus10"
# model <- "SpatPlus15"
# model <- "SpatPlus20"




#######################
# Load simulated data #
#######################

## Define the Scenario and Subscenario ##

Scenario <- 1
# Scenario <- 2
# Scenario <- 3

Subscenario <- 80
# Subscenario <- 50
# Subscenario <- 20


load(paste0("../../Simulated_data/Simu1_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))


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

load(paste0("Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))




#################################
# COVARIATE MODEL: eigenvectors #
#################################

n.sim <- 100

# SpatPlus5
f.Cov5 <- X1.weights ~ -1 + W.intercept + V1+ V2 + V3 + V4 + V5

# SpatPlus10
f.Cov10 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10

# SpatPlus15
f.Cov15 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        V11 + V12 + V13 + V14 + V15

# SpatPlus20
f.Cov20 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20



## Choose the covariate model ##

if (model=="SpatPlus5"){
  f.Cov <- f.Cov5
}

if (model=="SpatPlus10"){
  f.Cov <- f.Cov10
}

if (model=="SpatPlus15"){
  f.Cov <- f.Cov15
}

if (model=="SpatPlus20"){
  f.Cov <- f.Cov20
}


X1.Res <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))

  # Compute the weights
  W <- Spat[[i]]$summary.fitted.values$mode*Data$expected
  W.sqrt <- diag(sqrt(W))
  W.sqrt.inv <- diag(1/sqrt(W))
  
  # Covariate model
  W.eigen.vect <- W.sqrt%*%eigen.vect
  Data2 <- cbind(Data, W.eigen.vect)
  colnames(Data2) <- c(colnames(Data), paste0("V", 20:1))
  Data2$X1.weights <- W.sqrt%*%Data$X1
  Data2$W.intercept <- W.sqrt%*%rep(1, N)
  
  Covariate <- inla(f.Cov, family = "gaussian", data=Data2,
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                control.predictor=list(compute=TRUE),
                control.inla=list(strategy=strategy))

  
  X1.Res.eigen <- W.sqrt.inv%*%(Data2$X1.weights - Covariate$summary.fitted.values[, 1])
  X1.Res[[i]] <- as.vector(scale(X1.Res.eigen))
  
  Data2 <- NULL

}




########################
# SPATIAL+ MODEL: ICAR #
########################

f.SpatPlus <- O ~ 1 + X1.Res + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                 scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 



SpatPlus <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))
  # Simulated counts
  Data$O <- simu.O[[i]]
  Data$X1.Res <- X1.Res[[i]]

  SpatPlus[[i]] <- inla(f.SpatPlus, family="poisson", data=Data, E=expected,
                     control.predictor=list(compute=TRUE, cdf=c(log(1))),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     control.inla=list(strategy=strategy))
}




####################
# Save the results #
####################

if (model=="SpatPlus5"){
  SpatPlus5 <- SpatPlus
  save(SpatPlus5, file=paste0("Scenario", Scenario, "_SpatPlus5_cor", Subscenario, ".Rdata"))
}

if (model=="SpatPlus10"){
  SpatPlus10 <- SpatPlus
  save(SpatPlus10, file=paste0("Scenario", Scenario, "_SpatPlus10_cor", Subscenario, ".Rdata"))
}

if (model=="SpatPlus15"){
  SpatPlus15 <- SpatPlus
  save(SpatPlus15, file=paste0("Scenario", Scenario, "_SpatPlus15_cor", Subscenario, ".Rdata"))
}

if (model=="SpatPlus20"){
  SpatPlus20 <- SpatPlus
  save(SpatPlus20, file=paste0("Scenario", Scenario, "_SpatPlus20_cor", Subscenario, ".Rdata"))
  
}








