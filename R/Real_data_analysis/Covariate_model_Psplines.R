

########################################
# Construct the spatial B-spline basis #
########################################

p <- 3 	# Cubic B-splines
q <- 10	# Number of internal intervals (for Slovenia data we have used 27 intervals)


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




##############################
# Spatial structure matrices #
##############################

order <- 2   # 1=RW1, 2=RW2

D1 <- diff(diag(k1),differences=order)
P1 <- t(D1)%*%D1

D2 <- diff(diag(k2),differences=order)
P2 <- t(D2)%*%D2

R1 <- kronecker(diag(k2),P1)
R2 <- kronecker(P2,diag(k1))

Cmat.s <- list(inla.as.sparse(R1),inla.as.sparse(R2))




##############################
# COVARIATE MODEL: P-splines #
##############################

Data$X1.weights <- W.sqrt%*%Data$X1
Data$W.intercept <- W.sqrt%*%rep(1, N)
Bs.W <- W.sqrt%*%Bs


Data.splines <- list(X1.weights=Data$X1.weights, 
                     intercept=c(1,rep(NA,ks)),
                     ID.area=c(NA,1:ks))


f.CovP <- X1.weights ~ -1 + intercept + f(ID.area, model="generic3", Cmatrix=Cmat.s, constr=TRUE, diagonal=1e-6,
                                         hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))

Apredictor <- cbind(Data$W.intercept, Bs.W)


CovariateP <- inla(f.CovP, family="gaussian", data=Data.splines,
                  control.predictor=list(compute=TRUE, A=Apredictor),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  control.inla=list(strategy=strategy, verbose=TRUE))


Data$X1.ResP <- W.sqrt.inv%*%(Data$X1.weights - CovariateP$summary.fitted.values[1:N, 1])
Data$X1.ResP <- as.vector(scale(Data$X1.ResP))



