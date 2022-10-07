
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




##############
# NULL MODEL #
##############

f.Null <- O ~ X1

Null <- inla(f.Null, family="poisson", data=Data, E=E,
             control.predictor=list(compute=TRUE, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
             control.inla=list(strategy=strategy))


save(Null, file="Results_Uttar/Null_SexRatio.Rdata")
# save(Null, file="Results_Slovenia/Null_SE.Rdata")
# save(Null, file="Results_Scotland/Null_AFF.Rdata")




#################
# SPATIAL MODEL #
#################

f.Spat <- O ~ X1 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                     hyper=list(prec=list(prior=sdunif))) 

Spat <- inla(f.Spat, family="poisson", data=Data, E=E,
             control.predictor=list(compute=TRUE, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
             control.inla=list(strategy=strategy))


save(Spat, file="Results_Uttar/Spatial_SexRatio.Rdata")
# save(Spat, file="Results_Slovenia/Spatial_SE.Rdata")
# save(Spat, file="Results_Scotland/Spatial_AFF.Rdata")




#############
# RSR MODEL #
#############

ones.N <- rep(1, N)
Beta.df <- as.matrix(Data[, "X1"])


W <- diag(Spat$summary.fitted.values$mode*Data$E)
W.sqrt <- diag(sqrt(diag(W)))

X <- cbind(ones.N, as.matrix(Beta.df))
P <- W.sqrt%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W.sqrt
Pc <- diag(N)-P


eigen.Pc <- eigen(Pc)
L <- eigen.Pc$vectors[,eigen.Pc$values>1e-12]


eigen.Qs <- eigen(Q.xi)
Us <- eigen.Qs$vectors[,eigen.Qs$values>1e-12]
Ds <- diag(eigen.Qs$values[eigen.Qs$values>1e-12])


M <- solve(W.sqrt)%*%L%*%t(L)%*%W.sqrt
Z.area <- M%*%diag(N)

Data.RSR <- list(O=Data$O, E=Data$E,
                 X1=Data$X1,
                 Intercept=c(1,rep(NA,N+1)),
                 beta1=c(NA,1,rep(NA,N)),
                 ID.area=c(rep(NA,1+1),1:N))


f.RSR <- O ~ -1 + Intercept + beta1 + f(ID.area, model="besag", graph=Q.xi, 
                                        rankdef=1, constr=TRUE, hyper=list(prec=list(prior=sdunif)))


Apredictor <- cbind(rep(1,N), Beta.df, Z.area)

RSR <- inla(f.RSR, family="poisson", data=Data.RSR, E=E,
            control.predictor=list(compute=TRUE, A=Apredictor),
            control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
            control.inla=list(strategy=strategy),
            control.mode=list(theta=Spat$mode$theta, restart=FALSE),
            control.fixed=list(prec=0))


save(RSR, file="Results_Uttar/RSR_SexRatio.Rdata")
# save(RSR, file="Results_Slovenia/RSR_SE.Rdata")
# save(RSR, file="Results_Scotland/RSR_AFF.Rdata")

