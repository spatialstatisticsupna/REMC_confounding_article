

rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


n.sim <- 100 


## Define the Scenario and Subscenario ##

Scenario <- 1
# Scenario <- 2
# Scenario <- 3

Subscenario <- 80
# Subscenario <- 50
# Subscenario <- 20




############################################
# Load models fitted to the simulated data #
############################################

load(paste0("Scenario", Scenario, "_Null_cor", Subscenario, ".Rdata"))
load(paste0("Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))
load(paste0("Scenario", Scenario, "_RSR_cor", Subscenario, ".Rdata"))

load(paste0("Scenario", Scenario, "_SpatPlus5_cor", Subscenario, ".Rdata"))
load(paste0("Scenario", Scenario, "_SpatPlus10_cor", Subscenario, ".Rdata"))
load(paste0("Scenario", Scenario, "_SpatPlus15_cor", Subscenario, ".Rdata"))
load(paste0("Scenario", Scenario, "_SpatPlus20_cor", Subscenario, ".Rdata"))

load(paste0("Scenario", Scenario, "_SpatPlusP1_cor", Subscenario, ".Rdata"))
load(paste0("Scenario", Scenario, "_SpatPlusTP1_cor", Subscenario, ".Rdata"))
load(paste0("Scenario", Scenario, "_SpatPlusP2_cor", Subscenario, ".Rdata"))
load(paste0("Scenario", Scenario, "_SpatPlusTP2_cor", Subscenario, ".Rdata"))




#######################
# Load simulated data #
#######################

load(paste0("Simu1_data_Scenario", Scenario, "_cor", Subscenario, ".Rdata"))




########################################################
# TABLE A.1: Average value of MARB and MRRMSE of beta1 #
########################################################

real.beta1 <- 0.2
real.beta1.vect <- rep(real.beta1, n.sim)


beta1.mean <- function(x){
  data.frame(alpha.mean=x$summary.fixed[2,1])
}



b1.Null <- lapply(Null, beta1.mean)
b1.Spat <- lapply(Spat, beta1.mean)
b1.RSR <- lapply(RSR, beta1.mean)
b1.SpatPlus5 <- lapply(SpatPlus5, beta1.mean)
b1.SpatPlus10 <- lapply(SpatPlus10, beta1.mean)
b1.SpatPlus15 <- lapply(SpatPlus15, beta1.mean)
b1.SpatPlus20 <- lapply(SpatPlus20, beta1.mean)
b1.SpatPlusP1 <- lapply(SpatPlusP1, beta1.mean)
b1.SpatPlusTP1 <- lapply(SpatPlusTP1, beta1.mean)
b1.SpatPlusP2 <- lapply(SpatPlusP2, beta1.mean)
b1.SpatPlusTP2 <- lapply(SpatPlusTP2, beta1.mean)



MARB.beta1 <- function(Model){
  data.frame(MARB=abs(Reduce("+",mapply(function(x,y){(x-y)/y}, 
                                        x=Model, y=real.beta1.vect, SIMPLIFY=FALSE)))/n.sim)
}


MRRMSE.beta1 <- function(Model){
  data.frame(sqrt(Reduce("+",mapply(function(x,y){((x-y)/y)^2}, 
                                    x=Model, y=real.beta1.vect, SIMPLIFY=FALSE))/n.sim))
}


Table.MARB.beta1 <- data.frame(Null=MARB.beta1(b1.Null), Spat=MARB.beta1(b1.Spat), RSR=MARB.beta1(b1.RSR),
                               SpatPlus5=MARB.beta1(b1.SpatPlus5), SpatPlus10=MARB.beta1(b1.SpatPlus10),
                               SpatPlus15=MARB.beta1(b1.SpatPlus15), SpatPlus20=MARB.beta1(b1.SpatPlus20),
                               SpatPlusP1=MARB.beta1(b1.SpatPlusP1), SpatPlusTP1=MARB.beta1(b1.SpatPlusTP1),
                               SpatPlusP2=MARB.beta1(b1.SpatPlusP2), SpatPlusTP2=MARB.beta1(b1.SpatPlusTP2))


Table.MRRMSE.beta1 <- data.frame(Null=MRRMSE.beta1(b1.Null), Spat=MRRMSE.beta1(b1.Spat), RSR=MRRMSE.beta1(b1.RSR),
                                 SpatPlus5=MRRMSE.beta1(b1.SpatPlus5), SpatPlus10=MRRMSE.beta1(b1.SpatPlus10),
                                 SpatPlus15=MRRMSE.beta1(b1.SpatPlus15), SpatPlus20=MRRMSE.beta1(b1.SpatPlus20),
                                 SpatPlusP1=MRRMSE.beta1(b1.SpatPlusP1), SpatPlusTP1=MRRMSE.beta1(b1.SpatPlusTP1),
                                 SpatPlusP2=MRRMSE.beta1(b1.SpatPlusP2), SpatPlusTP2=MRRMSE.beta1(b1.SpatPlusTP2))


colnames(Table.MARB.beta1) <- c("Null", "Spatial", "RSR",
                                "SpatPlus5", "SpatPlus10", "SpatPlus15", "SpatPlus20",
                                "SpatPlusP1", "SpatPlusTP1",
                                "SpatPlusP2", "SpatPlusTP2")

colnames(Table.MRRMSE.beta1) <- colnames(Table.MARB.beta1)



Table.Accuracy <- rbind(Table.MARB.beta1, Table.MRRMSE.beta1)
rownames(Table.Accuracy) <- c("MARB", "MRRMSE")

Table.Accuracy[, 1:11] <- round(Table.Accuracy[, 1:11], 4)
Table.Accuracy <- data.frame(t(Table.Accuracy))
print(Table.Accuracy)




############################################################
# TABLE A.2: length of the 95% credible intervals of beta1 # 
############################################################

Length.beta1 <- function(Model){
  data.frame(length=abs(Model$summary.fixed$'0.975quant'[2]-Model$summary.fixed$'0.025quant'[2]))
}


Length <- data.frame(Null=mean(unlist(lapply(Null, function(x) Length.beta1(x)))),
                     Spat=mean(unlist(lapply(Spat, function(x) Length.beta1(x)))),
                     RSR=mean(unlist(lapply(RSR, function(x) Length.beta1(x)))),
                     SpatPlus5=mean(unlist(lapply(SpatPlus5, function(x) Length.beta1(x)))),
                     SpatPlus10=mean(unlist(lapply(SpatPlus10, function(x) Length.beta1(x)))),
                     SpatPlus15=mean(unlist(lapply(SpatPlus15, function(x) Length.beta1(x)))),
                     SpatPlus20=mean(unlist(lapply(SpatPlus20, function(x) Length.beta1(x)))),
                     SpatPlusP1=mean(unlist(lapply(SpatPlusP1, function(x) Length.beta1(x)))),
                     SpatPlusTP1=mean(unlist(lapply(SpatPlusTP1, function(x) Length.beta1(x)))),
                     SpatPlusP2=mean(unlist(lapply(SpatPlusP2, function(x) Length.beta1(x)))),
                     SpatPlusTP2=mean(unlist(lapply(SpatPlusTP2, function(x) Length.beta1(x)))))  



Length <- as.data.frame(t(Length))
colnames(Length) <- paste0("cor0.", Subscenario)
Length[, 1] <- round(Length[, 1], 4)
print(Length)




############################################
# TABLE A.3: average value of DIC and WAIC # 
############################################

DIC.inla <- function(x){
  data.frame(mean.deviance=x$dic$mean.deviance, ## posterior mean deviance
             p.eff=x$dic$p.eff,                 ## effective number of parameters
             DIC=x$dic$dic,                     ## Deviance Information Criterion
             WAIC=x$waic$waic)                  ## Watanabe-Akaike information criterion
}



Null.DIC <- apply(do.call(rbind, lapply(Null, DIC.inla)), 2, mean)
Spat.DIC <- apply(do.call(rbind, lapply(Spat, DIC.inla)), 2, mean)
RSR.DIC <- apply(do.call(rbind, lapply(RSR, DIC.inla)), 2, mean)
SpatPlus5.DIC <- apply(do.call(rbind, lapply(SpatPlus5, DIC.inla)), 2, mean)
SpatPlus10.DIC <- apply(do.call(rbind, lapply(SpatPlus10, DIC.inla)), 2, mean)
SpatPlus15.DIC <- apply(do.call(rbind, lapply(SpatPlus15, DIC.inla)), 2, mean)
SpatPlus20.DIC <- apply(do.call(rbind, lapply(SpatPlus20, DIC.inla)), 2, mean)
SpatPlusP1.DIC <- apply(do.call(rbind, lapply(SpatPlusP1, DIC.inla)), 2, mean)
SpatPlusTP1.DIC <- apply(do.call(rbind, lapply(SpatPlusTP1, DIC.inla)), 2, mean)
SpatPlusP2.DIC <- apply(do.call(rbind, lapply(SpatPlusP2, DIC.inla)), 2, mean)
SpatPlusTP2.DIC <- apply(do.call(rbind, lapply(SpatPlusTP2, DIC.inla)), 2, mean)


Table.DIC <- rbind(Null.DIC, Spat.DIC, RSR.DIC, SpatPlus5.DIC, SpatPlus10.DIC,
                   SpatPlus15.DIC, SpatPlus20.DIC, SpatPlusP1.DIC, SpatPlusTP1.DIC,
                   SpatPlusP2.DIC, SpatPlusTP2.DIC)

Table.DIC  <- as.data.frame(Table.DIC)
Table.DIC[, 1:4] <- round(Table.DIC[, 1:4], 4)
rownames(Table.DIC) <- c("Null", "Spatial", "RSR",
                         "SpatPlus5", "SpatPlus10", "SpatPlus15", "SpatPlus20",
                         "SpatPlusP1", "SpatPlusTP1", "SpatPlusP2", "SpatPlusTP2")
print(Table.DIC)





#####################################################################
# TABLE A.4: Average value of MARB and MRRMSE of the relative risks #
#####################################################################

Real.risk <- exp(log.risk)

Real.risk.vect <- vector("list", n.sim)
for (i in 1:n.sim){
  Real.risk.vect[[i]] <- Real.risk
}


RR.Null <- lapply(Null, function(x) matrix(x$summary.fitted.values$'0.5quant',70,1))
RR.Spat <- lapply(Spat, function(x) matrix(x$summary.fitted.values$'0.5quant',70,1))
RR.RSR <- lapply(RSR, function(x) matrix(x$summary.fitted.values$'0.5quant'[1:70],70,1))

RR.SpatPlus5 <- lapply(SpatPlus5, function(x) matrix(x$summary.fitted.values$'0.5quant',70,1))
RR.SpatPlus10 <- lapply(SpatPlus10, function(x) matrix(x$summary.fitted.values$'0.5quant',70,1))
RR.SpatPlus15 <- lapply(SpatPlus15, function(x) matrix(x$summary.fitted.values$'0.5quant',70,1))
RR.SpatPlus20 <- lapply(SpatPlus20, function(x) matrix(x$summary.fitted.values$'0.5quant',70,1))

RR.SpatPlusP1 <- lapply(SpatPlusP1, function(x) matrix(x$summary.fitted.values$'0.5quant',70,1))
RR.SpatPlusTP1 <- lapply(SpatPlusTP1, function(x) matrix(x$summary.fitted.values$'0.5quant',70,1))
RR.SpatPlusP2 <- lapply(SpatPlusP2, function(x) matrix(x$summary.fitted.values$'0.5quant'[1:70],70,1))
RR.SpatPlusTP2 <- lapply(SpatPlusTP2, function(x) matrix(x$summary.fitted.values$'0.5quant'[1:70],70,1))




MARB.RR <- function(Model){
  data.frame(MARB=mean(abs(Reduce("+",mapply(function(x,y){(x-y)/y}, 
                                             x=Model, y=Real.risk.vect, SIMPLIFY=FALSE)))/n.sim))
}


MRRMSE.RR <- function(Model){
  data.frame(mean(sqrt(Reduce("+",mapply(function(x,y){((x-y)/y)^2}, 
                                   x=Model, y=Real.risk.vect, SIMPLIFY=FALSE))/n.sim)))
}


Table.MARB.RR <- data.frame(Null=MARB.RR(RR.Null), Spat=MARB.RR(RR.Spat), RSR=MARB.RR(RR.RSR),
                   SpatPlus5=MARB.RR(RR.SpatPlus5), SpatPlus10=MARB.RR(RR.SpatPlus10),
                   SpatPlus15=MARB.RR(RR.SpatPlus15), SpatPlus20=MARB.RR(RR.SpatPlus20),
                   SpatPlusP1=MARB.RR(RR.SpatPlusP1), SpatPlusTP1=MARB.RR(RR.SpatPlusTP1),
                   SpatPlusP2=MARB.RR(RR.SpatPlusP2), SpatPlusTP2=MARB.RR(RR.SpatPlusTP2))


Table.MRRMSE.RR <- data.frame(Null=MRRMSE.RR(RR.Null), Spat=MRRMSE.RR(RR.Spat), RSR=MRRMSE.RR(RR.RSR),
                              SpatPlus5=MRRMSE.RR(RR.SpatPlus5), SpatPlus10=MRRMSE.RR(RR.SpatPlus10),
                              SpatPlus15=MRRMSE.RR(RR.SpatPlus15), SpatPlus20=MRRMSE.RR(RR.SpatPlus20),
                              SpatPlusP1=MRRMSE.RR(RR.SpatPlusP1), SpatPlusTP1=MRRMSE.RR(RR.SpatPlusTP1),
                              SpatPlusP2=MRRMSE.RR(RR.SpatPlusP2), SpatPlusTP2=MRRMSE.RR(RR.SpatPlusTP2))


colnames(Table.MARB.RR) <- c("Null", "Spatial", "RSR", 
                            "SpatPlus5", "SpatPlus10", "SpatPlus15", "SpatPlus20",
                            "SpatPlusP1", "SpatPlusTP1", "SpatPlusP2", "SpatPlusTP2")

colnames(Table.MRRMSE.RR) <- colnames(Table.MARB.RR)

Table.Accuracy <- rbind(Table.MARB.RR, Table.MRRMSE.RR)
rownames(Table.Accuracy) <- c("MARB", "MRRMSE")

Table.Accuracy[, 1:11] <- round(Table.Accuracy[, 1:11], 4)
Table.Accuracy <- data.frame(t(Table.Accuracy))
print(Table.Accuracy)



