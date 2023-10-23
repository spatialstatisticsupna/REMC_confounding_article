
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(INLA)
library(sf)
library(splines)




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




##############################################
# Load spatial models to compute the weights #
##############################################

load(paste0("Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))




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


f.Cov <- X1.weights ~ -1 + intercept + f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                                         hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))


X1.Res <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))

  # Compute the weights
  W <- Spat[[i]]$summary.fitted.values$mode*Data$expected
  W.sqrt <- diag(sqrt(W))
  W.sqrt.inv <- diag(1/sqrt(W))

  # Covariate model
  Data$X1.weights <- W.sqrt%*%Data$X1
  Data$W.intercept <- W.sqrt%*%rep(1, N)
  Bs.W <- W.sqrt%*%Bs
  
  Data.splines <- list(X1.weights=Data$X1.weights, 
                       intercept=c(1,rep(NA,ks)),
                       ID.area=c(NA,1:ks))
  
  Apredictor <- cbind(Data$W.intercept, Bs.W)
  
  Covariate <- inla(f.Cov, family="gaussian", data=Data.splines,
                    control.predictor=list(compute=TRUE, A=Apredictor),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy=strategy, verbose=TRUE))
  
  X1.Res.splines <- W.sqrt.inv%*%(Data$X1.weights - Covariate$summary.fitted.values[1:70, 1])
  X1.Res[[i]] <- as.vector(scale(X1.Res.splines))
}




########################
# SPATIAL+ MODEL: ICAR #
########################

f.SpatPlus <- O ~ 1 + X1.Res + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                 scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 



SpatPlusP1 <- vector("list", n.sim)


for (i in 1:n.sim){
  print(paste0("i=", i))
  # Simulated counts
  Data$O <- simu.O[[i]]
  Data$X1.Res <- X1.Res[[i]]
  
  SpatPlusP1[[i]] <- inla(f.SpatPlus, family="poisson", data=Data, E=expected,
                        control.predictor=list(compute=TRUE, cdf=c(log(1))),
                        control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                        control.inla=list(strategy=strategy))
}



####################
# Save the results #
####################

save(SpatPlusP1, file=paste0("Scenario", Scenario, "_SpatPlusP1_cor", Subscenario, ".Rdata"))



