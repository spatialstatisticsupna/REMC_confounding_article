

rm(list=ls())
setwd("")

# Packages
library(TMGMRF)
library(sf)
# library(rgdal)



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




##############
# TGMRF: GSC #
##############

family = "poisson"
type_data = "gamma-scale"  # "gamma-shape"=parameters in shape parameter
seed <- 1234


formula <- y ~ X1
df <- data.frame(y = Data$O, X1 = Data$X1)


nburnin <- 2000
thin <- 10
nsamp <- 10000


TGMRF.GSC <- tgmrf(
  formula = formula, data = df, E = Data$E, n=N,
  W = Q.xi.TGMRF,
  family = family, model = type_data,
  nchains = 3,
  nsamp = nsamp, nburnin = nburnin, thin = thin,
  seed = seed
)


save(TGMRF.GSC, file="Results_Uttar/TGMRF_GSC_SexRatio.Rdata")
# save(TGMRF.GSC, file="Results_Slovenia/TGMRF_GSC_SE.Rdata")
# save(TGMRF.GSC, file="Results_Scotland/TGMRF_GSC_AFF.Rdata")


