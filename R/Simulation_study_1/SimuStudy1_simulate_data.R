
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



#######################################################
# Load the dowry deaths data in Uttar Pradesh in 2001 #
#######################################################
load("../../Data/Dowry_death_2001.Rdata") 
N <- nrow(Data)



#############################
# Simulate the covariate X2 #
#############################
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}


# Simulate X2
if (Subscenario==80){
  rho <- 0.80
}
if (Subscenario==50){
  rho <- 0.50
}
if (Subscenario==20){
  rho <- 0.20
}

X2 <- complement(y=Data$X1, rho=rho)
cor(Data$X1, X2)


# Standardize the covariate X2
Data$X2 <- as.vector(scale(X2))



###############################################
# Simulate the additional spatial variability # 
###############################################
if(Scenario==1) {
  Data$spat <- rep(0,N)
}

if(Scenario==2) {
  var.u <- 0.20
  Sigma.u <- var.u*ginv(Q.xi)
  
  if (Subscenario==80){
    set.seed(1001)
  }
  if (Subscenario==50){
    set.seed(2005)
  }
  if (Subscenario==20){
    set.seed(3006)
  }
  spat <- mvrnorm(1,rep(0,N),Sigma.u)
  Data$spat <- spat
}

if(Scenario==3) {
  p <- 3 	# Cubic B-splines
  q <- 10	# Number of internal intervals
  
  # Marginal basis for longitude
  x1 <- coordinates(carto)[,1]
  x1 <- (x1-min(x1))/(max(x1)-min(x1))
  
  dist1 <- (max(x1)-min(x1))/q
  x1l <- min(x1)-dist1*0.05
  x1r <- max(x1)+dist1*0.05
  dx1 <- (x1r-x1l)/q
  knots1 <- seq(x1l-p*dx1, x1r+p*dx1, by=dx1)
  
  B1 <- spline.des(knots1,x1,p+1)$design
  k1 <- ncol(B1)
  
  # Marginal basis for latitude
  x2 <- coordinates(carto)[,2]
  x2 <- (x2-min(x2))/(max(x2)-min(x2))
  
  dist2 <- (max(x2)-min(x2))/q
  x2l <- min(x2)-dist2*0.05
  x2r <- max(x2)+dist2*0.05
  dx2 <- (x2r-x2l)/q
  knots2 <- seq(x2l-p*dx2, x2r+p*dx2, by=dx2)
  
  B2 <- spline.des(knots2,x2,p+1)$design
  k2 <- ncol(B2)
  
  # Row-wise Kronecker product
  Rten <- function(X1,X2){
    one1 <- matrix(1,1,ncol(X1))
    one2 <- matrix(1,1,ncol(X2))
    kronecker(X1,one2)*kronecker(one1,X2)
  }
  
  Bs <- Rten(B2,B1)
  ks <- ncol(Bs)

  # Spatial structure matrices
  order <- 2
  
  D1 <- diff(diag(k1),differences=order)
  P1 <- t(D1)%*%D1
  
  D2 <- diff(diag(k2),differences=order)
  P2 <- t(D2)%*%D2
  
  R1 <- kronecker(diag(k2),P1)
  R2 <- kronecker(P2,diag(k1))
  
  Cmat.s <- list(inla.as.sparse(R1),inla.as.sparse(R2))

  # Simulate the spatial random effects with B-splines
  lambda1 <- 1.22
  lambda2 <- 8.87
  
  Ps <- lambda1*R1 + lambda2*R2
  Sigma <- ginv(Ps)

  set.seed(1254) 
  spat <- Bs%*%mvrnorm(1,rep(0,169),Sigma)
  
  Data$spat <- spat
}  
  
cor(Data$X1, Data$spat)
  


#######################
# Simulate the counts #
#######################
log.risk <- 0.2*Data$X1 + 0.3*Data$X2 + Data$spat
lambda <-  Data$E*exp(log.risk)


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
     file=paste0("../../Simulated_data/Simu1_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))


