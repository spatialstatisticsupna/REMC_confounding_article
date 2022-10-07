

rm(list=ls())
setwd("")


# Load packages 
library(INLA)
library(splines)
library(sf)




#############################
# Load and prepare the data #
#############################

## Dowry death data 2001 ##
load("Data/Dowry_death_2001.Rdata")  # O=observed counts, E=expected counts, X1=standardized sex ratio
N <- nrow(Data)



## Slovenia stomach cancer data ##
# load("Data/Slovenia_stomach_cancer.Rdata")  # X=standardized socioeconomic indicator
# N <- nrow(Data)
# colnames(Data)[colnames(Data)=="X"] <- "X1"



## Scotland lip cancer data ##
# load("Data/Scotland_lip_cancer.Rdata")  # AFF=standardized AFF
# N <- nrow(Data)
# colnames(Data)[colnames(Data)=="AFF"] <- "X1"




###############################
# Hyperpriors for INLA models #
###############################

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

compute.patterns <- FALSE  ## Set compute.patterns=FALSE if posterior patterns are not required
strategy <- "laplace"




#############################
# Compute the weight matrix #
#############################

## Dowry death data 2001 ##
load("Results_Uttar/Spatial_SexRatio.Rdata")

## Slovenia stomach cancer data ##
# load("Results_Slovenia/Spatial_SE.Rdata")

## Scotland lip cancer data ##
# load("Results_Scotland/Spatial_AFF.Rdata")


weights <- Spat$summary.fitted.values$mode*Data$E
W.sqrt <- diag(sqrt(weights))





##############################
# COVARIATE MODEL: P-splines #
##############################

# ------------------------------------ #
# Construct the spatial B-spline basis #
# ------------------------------------ #

p <- 3 	# Cubic B-splines
q <- 10	# Number of internal intervals


## Dowry death data 2001: marginal basis for longitude ##
x1 <- coordinates(carto)[,1]

## Slovenia stomach cancer data: marginal basis for longitude ##
# x1 <- coord[,1]	

## Scotland lip cancer data: marginal basis for longitude ##
# centroids <- st_centroid(carto)
# coord <- st_coordinates(centroids)
# x1 <- coord[,1]	



x1 <- (x1-min(x1))/(max(x1)-min(x1))	

dist1 <- (max(x1)-min(x1))/q
x1l <- min(x1)-dist1*0.05
x1r <- max(x1)+dist1*0.05
dx1 <- (x1r-x1l)/q
knots1 <- seq(x1l-p*dx1, x1r+p*dx1, by=dx1)

B1 <- spline.des(knots1,x1,p+1)$design
k1 <- ncol(B1)



## Dowry death data 2001: marginal basis for latitude ##
x2 <- coordinates(carto)[,2]	

## Slovenia stomach cancer data: marginal basis for latitude ##
# x2 <- coord[,2]

## Scotland lip cancer data: marginal basis for latitude ##
# x2 <- coord[,2]



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

Data$X1.weights <- W.sqrt%*%Data$X1
Data$W.intercept <- W.sqrt%*%rep(1, N)
Bs.W <- W.sqrt%*%Bs


Data.splines <- list(X1.weights=Data$X1.weights, 
                     intercept=c(1,rep(NA,ks)),
                     ID.area=c(NA,1:ks))


f.Cov <- X1.weights ~ -1 + intercept + f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                                         hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))

Apredictor <- cbind(Data$W.intercept, Bs.W)


Covariate <- inla(f.Cov, family="gaussian", data=Data.splines,
                  control.predictor=list(compute=TRUE, A=Apredictor),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  control.inla=list(strategy=strategy, verbose=TRUE))


W.sqrt.inv <- diag(1/sqrt(weights))
Data$X1.Res <- W.sqrt.inv%*%(Data$X1.weights - Covariate$summary.fitted.values[1:N, 1])
Data$X1.Res <- as.vector(scale(Data$X1.Res))




#############################
# SPATIAL+ MODEL: P-splines #
#############################

f.SpatPlus <- O ~ -1 + intercept + beta1 +
                  f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                    hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))


Data.splines2 <- list(O=Data$O, E=Data$E,
                     X1.Res=Data$X1.Res,
                     intercept=c(1,rep(NA,ks+1)),
                     beta1=c(NA,1,rep(NA,ks)),
                     ID.area=c(NA,NA,1:ks))

Apredictor <- cbind(rep(1,N), Data$X1.Res, Bs)


SpatPlusP2 <- inla(f.SpatPlus, family="poisson", data=Data.splines2, E=E,
                   control.predictor=list(compute=TRUE, A=Apredictor),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.inla=list(strategy=strategy, verbose=TRUE))

save(SpatPlusP2, file="Results_Uttar/SpatPlusP2_SexRatio.Rdata")
# save(SpatPlusP2, file="Results_Slovenia/SpatPlusP2_SE.Rdata")
# save(SpatPlusP2, file="Results_Scotland/SpatPlusP2_AFF.Rdata")
