
rm(list=ls())
setwd("")


# Load packages
library(INLA)




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




########################################################
# Eigenvectors and eigenvalues of the adjacency matrix #
########################################################

eigendecomp <- eigen(Q.xi)
eigendecomp$values


## Dowry death data 2001 ##
# Eigenvectors corresponding to 20 lowest non-null eigenvalues
eigen.vect <- eigendecomp$vectors[, (N-20):(N-1)]


## Slovenia stomach cancer data ##
# Eigenvectors corresponding to 40 lowest non-null eigenvalues
# eigen.vect <- eigendecomp$vectors[, (N-40):(N-1)]  


## Scotland lip cancer data ##
# Eigenvectors corresponding to 15 lowest non-null eigenvalues
# eigen.vect <- eigendecomp$vectors[, (N-15):(N-1)]




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




#################################
# COVARIATE MODEL: eigenvectors #
#################################

W.eigen.vect <- W.sqrt%*%eigen.vect
Data2 <- cbind(Data, W.eigen.vect)


## Dowry death data 2001 ##
colnames(Data2) <- c(colnames(Data), paste0("V", 20:1))


## Slovenia stomach cancer data ##
# colnames(Data2) <- c(colnames(Data), paste0("V", 40:1))


## Scotland lip cancer data ##
# colnames(Data2) <- c(colnames(Data), paste0("V", 15:1))


Data2$X1.weights <- W.sqrt%*%Data$X1
Data2$W.intercept <- W.sqrt%*%rep(1, N)



## Choose the number of eigenvectors to include in the covariate model ##

# SpatPlus5
f.Cov5 <- X1.weights ~ -1 + W.intercept + V1+ V2 + V3 + V4 + V5

# SpatPlus10
# f.Cov10 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10

# SpatPlus15
# f.Cov15 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        # V11 + V12 + V13 + V14 + V15

# SpatPlus20
# f.Cov20 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        # V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20

# SpatPlus30
# f.Cov30 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        # V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + 
                        # V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30

# SpatPlus40
# f.Cov40 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        # V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + 
                        # V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + V31 + V32 + V33 + V34 + 
                        # V35 + V36 + V37 + V38 + V39 + V40



Covariate <- inla(f.Cov5, family = "gaussian", data=Data2, 
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  control.predictor=list(compute=TRUE),
                  control.inla=list(strategy="laplace"))


W.sqrt.inv <- diag(1/sqrt(weights))
Data$X1.Res <- W.sqrt.inv%*%(Data2$X1.weights - Covariate$summary.fitted.values[, 1])
Data$X1.Res <- as.vector(scale(Data$X1.Res))




########################
# SPATIAL+ MODEL: ICAR #
########################

f.SpatPlus <- O ~ 1 + X1.Res + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                 scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 


SpatPlus <- inla(f.SpatPlus, family="poisson", data=Data, E=E,
                 control.predictor=list(compute=TRUE, cdf=c(log(1))),
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 control.inla=list(strategy=strategy))


SpatPlus5 <- SpatPlus
save(SpatPlus5, file="Results_Uttar/SpatPlus5_SexRatio.Rdata")
# save(SpatPlus5, file="Results_Slovenia/SpatPlus5_SE.Rdata")
# save(SpatPlus5, file="Results_Scotland/SpatPlus5_AFF.Rdata")


# SpatPlus10 <- SpatPlus
# save(SpatPlus10, file="Results_Uttar/SpatPlus10_SexRatio.Rdata")
# save(SpatPlus10, file="Results_Slovenia/SpatPlus10_SE.Rdata")
# save(SpatPlus10, file="Results_Scotland/SpatPlus10_AFF.Rdata")


# SpatPlus15 <- SpatPlus
# save(SpatPlus15, file="Results_Uttar/SpatPlus15_SexRatio.Rdata")
# save(SpatPlus10, file="Results_Slovenia/SpatPlus15_SE.Rdata")
# save(SpatPlus15, file="Results_Scotland/SpatPlus15_AFF.Rdata")


# SpatPlus20 <- SpatPlus
# save(SpatPlus20, file="Results_Uttar/SpatPlus20_SexRatio.Rdata")
# save(SpatPlus20, file="Results_Slovenia/SpatPlus20_SE.Rdata")


# SpatPlus30 <- SpatPlus
# save(SpatPlus30, file="Results_Slovenia/SpatPlus30_SE.Rdata")


# SpatPlus40 <- SpatPlus
# save(SpatPlus40, file="Results_Slovenia/SpatPlus40_SE.Rdata")





