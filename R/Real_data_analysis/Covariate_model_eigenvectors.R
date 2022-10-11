

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




#################################################
# Multiply X1 and eigenvectors with the weights #
#################################################

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




###################################
# COVARIATE MODEL: 5 eigenvectors #
###################################

f.Cov5 <- X1.weights ~ -1 + W.intercept + V1+ V2 + V3 + V4 + V5


Covariate5 <- inla(f.Cov5, family = "gaussian", data=Data2, 
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.predictor=list(compute=TRUE),
                   control.inla=list(strategy="laplace"))


Data$X1.Res5 <- W.sqrt.inv%*%(Data2$X1.weights - Covariate5$summary.fitted.values[, 1])
Data$X1.Res5 <- as.vector(scale(Data$X1.Res5))




####################################
# COVARIATE MODEL: 10 eigenvectors #
####################################

f.Cov10 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10


Covariate10 <- inla(f.Cov10, family = "gaussian", data=Data2, 
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.predictor=list(compute=TRUE),
                    control.inla=list(strategy="laplace"))


Data$X1.Res10 <- W.sqrt.inv%*%(Data2$X1.weights - Covariate10$summary.fitted.values[, 1])
Data$X1.Res10 <- as.vector(scale(Data$X1.Res10))




####################################
# COVARIATE MODEL: 15 eigenvectors #
####################################

f.Cov15 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        V11 + V12 + V13 + V14 + V15


Covariate15 <- inla(f.Cov15, family = "gaussian", data=Data2, 
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.predictor=list(compute=TRUE),
                    control.inla=list(strategy="laplace"))


Data$X1.Res15 <- W.sqrt.inv%*%(Data2$X1.weights - Covariate15$summary.fitted.values[, 1])
Data$X1.Res15 <- as.vector(scale(Data$X1.Res15))




####################################
# COVARIATE MODEL: 20 eigenvectors #
####################################

f.Cov20 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20


Covariate20 <- inla(f.Cov20, family = "gaussian", data=Data2, 
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.predictor=list(compute=TRUE),
                    control.inla=list(strategy="laplace"))


Data$X1.Res20 <- W.sqrt.inv%*%(Data2$X1.weights - Covariate20$summary.fitted.values[, 1])
Data$X1.Res20 <- as.vector(scale(Data$X1.Res20))




####################################
# COVARIATE MODEL: 30 eigenvectors #
####################################

# f.Cov30 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                          # V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + 
                          # V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30


# Covariate30 <- inla(f.Cov30, family = "gaussian", data=Data2, 
                    # control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    # control.predictor=list(compute=TRUE),
                    # control.inla=list(strategy="laplace"))


# Data$X1.Res30 <- W.sqrt.inv%*%(Data2$X1.weights - Covariate30$summary.fitted.values[, 1])
# Data$X1.Res30 <- as.vector(scale(Data$X1.Res30))




####################################
# COVARIATE MODEL: 40 eigenvectors #
####################################

# f.Cov40 <- X1.weights ~ -1 + W.intercept + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                        # V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + 
                        # V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + V31 + V32 + V33 + V34 + 
                        # V35 + V36 + V37 + V38 + V39 + V40


# Covariate40 <- inla(f.Cov40, family = "gaussian", data=Data2, 
                    # control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    # control.predictor=list(compute=TRUE),
                    # control.inla=list(strategy="laplace"))


# Data$X1.Res40 <- W.sqrt.inv%*%(Data2$X1.weights - Covariate40$summary.fitted.values[, 1])
# Data$X1.Res40 <- as.vector(scale(Data$X1.Res40))



