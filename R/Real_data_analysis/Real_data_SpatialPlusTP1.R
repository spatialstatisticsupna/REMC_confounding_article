
rm(list=ls())
setwd("")

# Load packages 
library(INLA)
library(mgcv)
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




###############################
# COVARIATE MODEL: TP-splines #
###############################

method<-"GCV.Cp" #"REML"


## Dowry death data 2001 ##
coord <- coordinates(carto)
coord.df <-  as.data.frame(coord)


## Slovenia stomach cancer data ##
# coord.df <-  as.data.frame(coord)


## Scotland lip cancer data ##
# centroids <- st_centroid(carto)
# coord <- st_coordinates(centroids)
# coord.df <-  as.data.frame(coord[, c("X", "Y")])


names(coord.df)<-c("x","y")


Covariate <- gam(Data$X1~1+s(x,y,k=17, fx=TRUE, m=2), weights=weights, data=coord.df)$fitted.values
X1.Res <- Data$X1 - Covariate
Data$X1.Res <- as.vector(scale(X1.Res))




########################
# SPATIAL+ MODEL: ICAR #
########################

f.SpatPlus <- O ~ 1 + X1.Res + f(ID.area, model="besag", graph=Q.xi, constr=TRUE, 
                                 scale.model=FALSE, hyper=list(prec=list(prior=sdunif))) 


SpatPlusTP1 <- inla(f.SpatPlus, family="poisson", data=Data, E=E,
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy=strategy))


save(SpatPlusTP1, file="Results_Uttar/SpatPlusTP1_SexRatio.Rdata")
# save(SpatPlusTP1, file="Results_Slovenia/SpatPlusTP1_SE.Rdata")
# save(SpatPlusTP1, file="Results_Scotland/SpatPlusTP1_AFF.Rdata")
