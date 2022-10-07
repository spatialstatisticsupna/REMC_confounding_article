
rm(list=ls())
setwd("")


# Load packages
library(INLA)
library(sf)
library(splines)




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




##############################################
# Load spatial models to compute the weights #
##############################################

load(paste0("Results_SimuStudy2/Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))




##############################
# COVARIATE MODEL: P-splines #
##############################

# ------------------------------------ #
# Construct the spatial B-spline basis #
# ------------------------------------ #

p <- 3 	# Cubic B-splines
q <- 10	# Number of internal intervals


## Marginal basis for longitude ##	

x1 <- coordinates(carto)[,1]
x1 <- (x1-min(x1))/(max(x1)-min(x1))	## into the [0,1] interval

dist1 <- (max(x1)-min(x1))/q
x1l <- min(x1)-dist1*0.05
x1r <- max(x1)+dist1*0.05
dx1 <- (x1r-x1l)/q
knots1 <- seq(x1l-p*dx1, x1r+p*dx1, by=dx1)

B1 <- spline.des(knots1,x1,p+1)$design
k1 <- ncol(B1)



## Marginal basis for latitude ##

x2 <- coordinates(carto)[,2]		
x2 <- (x2-min(x2))/(max(x2)-min(x2))	## into the [0,1] interval

dist2 <- (max(x2)-min(x2))/q
x2l <- min(x2)-dist2*0.05
x2r <- max(x2)+dist2*0.05
dx2 <- (x2r-x2l)/q
knots2 <- seq(x2l-p*dx2, x2r+p*dx2, by=dx2)

B2 <- spline.des(knots2,x2,p+1)$design
k2 <- ncol(B2)



## Row-wise Kronecker product ##
Rten <- function(X1,X2){
  one1 <- matrix(1,1,ncol(X1))
  one2 <- matrix(1,1,ncol(X2))
  kronecker(X1,one2)*kronecker(one1,X2)
}

Bs <- Rten(B2,B1)
ks <- ncol(Bs)




# -------------------------- #
# Spatial structure matrices #
# -------------------------- #

order <- 2   # 1=RW1, 2=RW2

D1 <- diff(diag(k1),differences=order)
P1 <- t(D1)%*%D1

D2 <- diff(diag(k2),differences=order)
P2 <- t(D2)%*%D2

R1 <- kronecker(diag(k2),P1)
R2 <- kronecker(P2,diag(k1))

Cmat.s <- list(inla.as.sparse(R1),inla.as.sparse(R2))




# -------------------------- #
# COVARIATE MODEL: P-splines #
# -------------------------- #


n.sim <- 100


f.X1.Cov <- X1.weights ~ -1 + intercept + f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                                            hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))

f.X2.Cov <- X2.weights ~ -1 + intercept + f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                                            hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))


X1.Res <- vector("list", n.sim)
X2.Res <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))
  
  # Compute the weights
  W <- Spat[[i]]$summary.fitted.values$mode*Data$expected
  W.sqrt <- diag(sqrt(W))
  W.sqrt.inv <- diag(1/sqrt(W))
  
  # Compute the residuals
  Data$X1.weights <- W.sqrt%*%Data$X1
  Data$X2.weights <- W.sqrt%*%Data$X2
  Data$W.intercept <- W.sqrt%*%rep(1, N)
  Bs.W <- W.sqrt%*%Bs
  
  Data.splines.X1 <- list(X1.weights=Data$X1.weights, 
                          intercept=c(1,rep(NA,ks)),
                          ID.area=c(NA,1:ks))
  
  Data.splines.X2 <- list(X2.weights=Data$X2.weights, 
                          intercept=c(1,rep(NA,ks)),
                          ID.area=c(NA,1:ks))
  
  Apredictor <- cbind(Data$W.intercept, Bs.W)
  
  Covariate.X1 <- inla(f.X1.Cov, family="gaussian", data=Data.splines.X1,
                       control.predictor=list(compute=TRUE, A=Apredictor),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                       control.inla=list(strategy=strategy, verbose=TRUE))
  
  Covariate.X2 <- inla(f.X2.Cov, family="gaussian", data=Data.splines.X2,
                       control.predictor=list(compute=TRUE, A=Apredictor),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                       control.inla=list(strategy=strategy, verbose=TRUE))
  
  X1.Res.splines <- W.sqrt.inv%*%(Data$X1.weights - Covariate.X1$summary.fitted.values[1:70, 1])
  X1.Res[[i]] <- as.vector(scale(X1.Res.splines))
  
  X2.Res.splines <- W.sqrt.inv%*%(Data$X2.weights - Covariate.X2$summary.fitted.values[1:70, 1])
  X2.Res[[i]] <- as.vector(scale(X2.Res.splines))
}




#############################
# SPATIAL+ MODEL: P-splines #
#############################

f.SpatPlus <- O ~ -1 + intercept + beta1 + beta2 +
                  f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                    hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))


SpatPlusP2 <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))

  Data.splines <- list(O=simu.O[[i]], E=Data$expected,
                       X1.Res=X1.Res[[i]],
                       X2.Res=X2.Res[[i]],
                       intercept=c(1,rep(NA,ks+2)),
                       beta1=c(NA,1,rep(NA,ks+1)),
                       beta2=c(NA,NA,1,rep(NA,ks)),
                       ID.area=c(NA,NA,NA,1:ks))
  
  Apredictor <- cbind(rep(1,N), X1.Res[[i]], X2.Res[[i]], Bs)

  # Fit the model
  SpatPlusP2[[i]] <- inla(f.SpatPlus, family="poisson", data=Data.splines, E=E,
                     control.predictor=list(compute=TRUE, A=Apredictor),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     control.inla=list(strategy=strategy, verbose=TRUE))
}




####################
# Save the results #
####################

save(SpatPlusP2, file=paste0("Results_SimuStudy2/Scenario", Scenario, "_SpatPlusP2_cor", Subscenario, ".Rdata"))



