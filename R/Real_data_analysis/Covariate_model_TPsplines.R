

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


CovariateTP <- gam(Data$X1~1+s(x,y,k=17, fx=TRUE, m=2), weights=diag(W), data=coord.df)$fitted.values  # for Slovenia data k=30 is chosen
X1.ResTP <- Data$X1 - CovariateTP
Data$X1.ResTP <- as.vector(scale(X1.ResTP))




