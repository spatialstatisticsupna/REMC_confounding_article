
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(INLA)
library(splines)
library(sf)
library(mgcv)


# Define the dataset
dataset <- "Dowry"   
# dataset <- "Slovenia"
# dataset <- "Scotland"



#############################
# Load and prepare the data #
#############################

## Dowry death data 2001 ##
if(dataset=="Dowry") {
  load("../../Data/Dowry_death_2001.Rdata")  # O=observed counts, E=expected counts, X1=standardized sex ratio
  N <- nrow(Data)
}

## Slovenia stomach cancer data ##
if(dataset=="Slovenia") {
  load("../../Data/Slovenia_stomach_cancer.Rdata")  # X=standardized socioeconomic indicator
  N <- nrow(Data)
  colnames(Data)[colnames(Data)=="X"] <- "X1"
}

## Scotland lip cancer data ##
if(dataset=="Scotland") {
  load("../../Data/Scotland_lip_cancer.Rdata")  # AFF=standardized AFF
  N <- nrow(Data)
  colnames(Data)[colnames(Data)=="AFF"] <- "X1"
}



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




#################
# SPATIAL MODEL #
#################

f.Spat <- O ~ X1 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                     hyper=list(prec=list(prior=sdunif))) 

Spat <- inla(f.Spat, family="poisson", data=Data, E=E,
             control.predictor=list(compute=TRUE, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
             control.inla=list(strategy=strategy))




#############
# RSR MODEL #
#############

ones.N <- rep(1, N)
Beta.df <- as.matrix(Data[, "X1"])


W <- diag(Spat$summary.fitted.values$mode*Data$E)
W.sqrt <- diag(sqrt(diag(W)))
W.sqrt.inv <- diag(1/sqrt(diag(W)))

X <- cbind(ones.N, as.matrix(Beta.df))
P <- W.sqrt%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W.sqrt
Pc <- diag(N)-P


eigen.Pc <- eigen(Pc)
L <- eigen.Pc$vectors[,eigen.Pc$values>1e-12]


eigen.Qs <- eigen(Q.xi)
Us <- eigen.Qs$vectors[,eigen.Qs$values>1e-12]
Ds <- diag(eigen.Qs$values[eigen.Qs$values>1e-12])


M <- W.sqrt.inv%*%L%*%t(L)%*%W.sqrt
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




#####################################################
# SPATIAL+ MODEL: Covariate model with eigenvectors #
#####################################################

source("Covariate_model_eigenvectors.R")



###################
# SpatPlus5 MODEL #
###################

f.SpatPlus5 <- O ~ 1 + X1.Res5 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                   scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 


SpatPlus5 <- inla(f.SpatPlus5, family="poisson", data=Data, E=E,
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  control.inla=list(strategy=strategy))




####################
# SpatPlus10 MODEL #
####################

f.SpatPlus10 <- O ~ 1 + X1.Res10 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                   scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 


SpatPlus10 <- inla(f.SpatPlus10, family="poisson", data=Data, E=E,
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  control.inla=list(strategy=strategy))




####################
# SpatPlus15 MODEL #
####################

f.SpatPlus15 <- O ~ 1 + X1.Res15 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                     scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 


SpatPlus15 <- inla(f.SpatPlus15, family="poisson", data=Data, E=E,
                   control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.inla=list(strategy=strategy))




####################
# SpatPlus20 MODEL #
####################

if(dataset %in% c("Dowry", "Slovenia")) {
  f.SpatPlus20 <- O ~ 1 + X1.Res20 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                       scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 
  
  
  SpatPlus20 <- inla(f.SpatPlus20, family="poisson", data=Data, E=E,
                     control.predictor=list(compute=TRUE, cdf=c(log(1))),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     control.inla=list(strategy=strategy))
}





####################
# SpatPlus30 MODEL #
####################

if(dataset=="Slovenia") {
  f.SpatPlus30 <- O ~ 1 + X1.Res30 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                       scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 
  
  SpatPlus30 <- inla(f.SpatPlus30, family="poisson", data=Data, E=E,
                     control.predictor=list(compute=TRUE, cdf=c(log(1))),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     control.inla=list(strategy=strategy))
}



####################
# SpatPlus40 MODEL #
####################

if(dataset=="Slovenia") {
 f.SpatPlus40 <- O ~ 1 + X1.Res40 + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                      scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 

 SpatPlus40 <- inla(f.SpatPlus40, family="poisson", data=Data, E=E,
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy=strategy))
}



##################################################
# SPATIAL+ MODEL: Covariate model with P-splines #
##################################################

source("Covariate_model_Psplines.R")



####################
# SpatPlusP1 MODEL #
####################

f.SpatPlusP1 <- O ~ 1 + X1.ResP + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                    scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 


SpatPlusP1 <- inla(f.SpatPlusP1, family="poisson", data=Data, E=E,
                   control.predictor=list(compute=TRUE, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.inla=list(strategy=strategy))




####################
# SpatPlusP2 MODEL #
####################

f.SpatPlusP2 <- O ~ -1 + intercept + beta1 +
                    f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                      hyper=list(prec1=list(prior=sdunif), prec2=list(prior=sdunif)))


Data.splines2 <- list(O=Data$O, E=Data$E,
                      X1.ResP=Data$X1.ResP,
                      intercept=c(1,rep(NA,ks+1)),
                      beta1=c(NA,1,rep(NA,ks)),
                      ID.area=c(NA,NA,1:ks))

Apredictor <- cbind(rep(1,N), Data$X1.ResP, Bs)


SpatPlusP2 <- inla(f.SpatPlusP2, family="poisson", data=Data.splines2, E=E,
                   control.predictor=list(compute=TRUE, A=Apredictor),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.inla=list(strategy=strategy, verbose=TRUE))




###################################################
# SPATIAL+ MODEL: Covariate model with TP-splines #
###################################################

source("Covariate_model_TPsplines.R")



#####################
# SpatPlusTP1 MODEL #
#####################

f.SpatPlusTP1 <- O ~ 1 + X1.ResTP + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                      scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 


SpatPlusTP1 <- inla(f.SpatPlusTP1, family="poisson", data=Data, E=E,
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy=strategy))




#####################
# SpatPlusTP2 MODEL #
#####################

f.SpatPlusTP2 <- O ~ -1 + intercept + beta1 +
                     f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                       hyper=list(prec1=list(prior=sdunif), prec2=list(prior=sdunif)))


Data.splines2 <- list(O=Data$O, E=Data$E,
                      X1.ResTP=Data$X1.ResTP,
                      intercept=c(1,rep(NA,ks+1)),
                      beta1=c(NA,1,rep(NA,ks)),
                      ID.area=c(NA,NA,1:ks))

Apredictor <- cbind(rep(1,N), Data$X1.ResTP, Bs)


SpatPlusTP2 <- inla(f.SpatPlusTP2, family="poisson", data=Data.splines2, E=E,
                    control.predictor=list(compute=TRUE, A=Apredictor),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy=strategy, verbose=TRUE))




############################################################################
# TABLE 2 (Table 3, Table 4): Posterior mean, posterior standard deviation # 
#                             and 95% credible intervals of beta1          #
############################################################################

if(dataset=="Dowry") {
  Models <- list(Null, Spat, RSR, SpatPlus5, SpatPlus10, SpatPlus15, SpatPlus20,
                 SpatPlusP1, SpatPlusTP1, SpatPlusP2, SpatPlusTP2)

  names(Models) <- c("Null", "Spatial", "RSR", "SpatPlus5", "SpatPlus10",
                     "SpatPlus15", "SpatPlus20", "SpatPlusP1", "SpatPlusTP1",
                     "SpatPlusP2", "SpatPlusTP2")
}


if(dataset=="Slovenia") {
  Models <- list(Null, Spat, RSR, SpatPlus5, SpatPlus10, SpatPlus15, SpatPlus20,
                 SpatPlus30, SpatPlus40,
                 SpatPlusP1, SpatPlusTP1, SpatPlusP2, SpatPlusTP2)
  
  names(Models) <- c("Null", "Spatial", "RSR", "SpatPlus5", "SpatPlus10",
                     "SpatPlus15", "SpatPlus20", "SpatPlus30", "SpatPlus40",
                     "SpatPlusP1", "SpatPlusTP1", "SpatPlusP2", "SpatPlusTP2")
}


if(dataset=="Scotland") {
  Models <- list(Null, Spat, RSR, SpatPlus5, SpatPlus10, SpatPlus15,
                 SpatPlusP1, SpatPlusTP1, SpatPlusP2, SpatPlusTP2)
  
  names(Models) <- c("Null", "Spatial", "RSR", "SpatPlus5", "SpatPlus10",
                     "SpatPlus15",
                     "SpatPlusP1", "SpatPlusTP1", "SpatPlusP2", "SpatPlusTP2")
}


Tab.beta1 <- do.call(rbind,lapply(Models, function(x) x$summary.fixed[2, c(1:3, 5)]))
round(Tab.beta1,4)




####################
# SAVE THE RESULTS #
####################

if(dataset=="Dowry") {
  save(Models, Tab.beta1, file="Results_dowry_deaths.Rdata")
}

if(dataset=="Slovenia") {
  save(Models, Tab.beta1, file="Results_Slovenia_data.Rdata")
}

if(dataset=="Scotland") {
  save(Models, Tab.beta1, file="Results_Scotland_data.Rdata")
}







