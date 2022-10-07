
rm(list=ls())
setwd("")


n.sim <- 100 


## Define the Scenario and Subscenario ##
Scenario <- 1  # 2,3
Subscenario <- 80  # 50,20




############################################
# Load models fitted to the simulated data #
############################################

load(paste0("Results_SimuStudy1/Scenario", Scenario, "_Null_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_Spatial_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_RSR_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_summary_TGMRF_GSC_cor", Subscenario, ".Rdata"))

load(paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlus5_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlus10_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlus15_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlus20_cor", Subscenario, ".Rdata"))

load(paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlusP1_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlusTP1_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlusP2_cor", Subscenario, ".Rdata"))
load(paste0("Results_SimuStudy1/Scenario", Scenario, "_SpatPlusTP2_cor", Subscenario, ".Rdata"))


save.image(paste0("Results_SimuStudy1/Scenario", Scenario, "_all_models_cor", Subscenario, ".Rdata"))




#############################################################
# TABLE 6: posterior means and standard deviations of beta1 #
#############################################################

beta1 <- function(x){
 data.frame(alpha.mean=x$summary.fixed[2,1], 
            alpha.sd=x$summary.fixed[2,2])
}

beta1.nimble <- function(x){
  data.frame(alpha.mean=x$statistics[grep('X1', rownames(x$statistics), value=TRUE), 1],
             alpha.sd=x$statistics[grep('X1', rownames(x$statistics), value=TRUE), 2])
}



Null.beta1 <- rbind(apply(do.call(rbind, lapply(Null, beta1)), 2, mean))
Spat.beta1 <- rbind(apply(do.call(rbind, lapply(Spat, beta1)), 2, mean))
RSR.beta1 <- rbind(apply(do.call(rbind, lapply(RSR, beta1)), 2, mean))
TGMRF.beta1 <- rbind(apply(do.call(rbind, lapply(TGMRF.GSC, beta1.nimble)), 2, mean))


SpatPlus5.beta1 <- rbind(apply(do.call(rbind, lapply(SpatPlus5, beta1)), 2, mean))
SpatPlus10.beta1 <- rbind(apply(do.call(rbind, lapply(SpatPlus10, beta1)), 2, mean))
SpatPlus15.beta1 <- rbind(apply(do.call(rbind, lapply(SpatPlus15, beta1)), 2, mean))
SpatPlus20.beta1 <- rbind(apply(do.call(rbind, lapply(SpatPlus20, beta1)), 2, mean))

SpatPlusP1.beta1 <- rbind(apply(do.call(rbind, lapply(SpatPlusP1, beta1)), 2, mean))
SpatPlusTP1.beta1 <- rbind(apply(do.call(rbind, lapply(SpatPlusTP1, beta1)), 2, mean))
SpatPlusP2.beta1 <- rbind(apply(do.call(rbind, lapply(SpatPlusP2, beta1)), 2, mean))
SpatPlusTP2.beta1 <- rbind(apply(do.call(rbind, lapply(SpatPlusTP2, beta1)), 2, mean))



Tab.beta1 <- rbind(Null.beta1, Spat.beta1, RSR.beta1, TGMRF.beta1,
                   SpatPlus5.beta1, SpatPlus10.beta1, SpatPlus15.beta1, 
                   SpatPlus20.beta1, SpatPlusP1.beta1, SpatPlusTP1.beta1,
                   SpatPlusP2.beta1, SpatPlusTP2.beta1)


Tab.beta1 <- as.data.frame(Tab.beta1)
rownames(Tab.beta1) <- c("Null", "Spatial", "RSR", "TGRMF.GSC", 
                         "SpatPlus5", "SpatPlus10", "SpatPlus15", "SpatPlus20",
                         "SpatPlusP1", "SpatPlusTP1",
                         "SpatPlusP2", "SpatPlusTP2")

colnames(Tab.beta1) <- c("mean", "sd")
Tab.beta1[, 1:2] <- round(Tab.beta1[, 1:2], 4)
print("Table beta1")
print(Tab.beta1)




###############################################################
# TABLE 7: simulated and estimated estandard errors for beta1 # 
###############################################################

beta1.mean <- function(x){
  data.frame(alpha.mean=x$summary.fixed[2,1])
}

beta1.mean.nimble <- function(x){
  data.frame(alpha.mean=x$statistics[grep('X1', rownames(x$statistics), value=TRUE), 1])
}


est.beta1 <- function(Model){
  data.frame(apply(do.call(rbind, lapply(Model, beta1))[2], 2, mean))
}


sim.beta1 <- function(Model){
  all.mean <- rep(rbind(apply(do.call(rbind, lapply(Model, beta1.mean)), 2, mean)), n.sim)
  data.frame(sqrt(Reduce("+",mapply(function(x,y){(x-y)^2}, 
                                    x=lapply(Model, beta1.mean), y=all.mean, SIMPLIFY=FALSE))/n.sim))
}


est.beta1.nimble <- function(Model){
  data.frame(apply(do.call(rbind, lapply(Model, beta1.nimble))[2], 2, mean))
}


sim.beta1.nimble <- function(Model){
  all.mean <- rep(rbind(apply(do.call(rbind, lapply(Model, beta1.mean.nimble)), 2, mean)), n.sim)
  data.frame(sqrt(Reduce("+",mapply(function(x,y){(x-y)^2}, 
                                    x=lapply(Model, beta1.mean.nimble), y=all.mean, SIMPLIFY=FALSE))/n.sim))
}



Table.est <- data.frame(Null=est.beta1(Null), Spat=est.beta1(Spat), RSR=est.beta1(RSR),
                        TGMRF=est.beta1.nimble(TGMRF.GSC), SpatPlus5=est.beta1(SpatPlus5),
                        SpatPlus10=est.beta1(SpatPlus10), SpatPlus15=est.beta1(SpatPlus15),
                        SpatPlus20=est.beta1(SpatPlus20), SpatPlusP1=est.beta1(SpatPlusP1),
                        SpatPlusTP1=est.beta1(SpatPlusTP1), SpatPlusP2=est.beta1(SpatPlusP2),
                        SpatPlusTP2=est.beta1(SpatPlusTP2))


Table.sim <- data.frame(Null=sim.beta1(Null), Spat=sim.beta1(Spat), RSR=sim.beta1(RSR),
                        TGMRF=sim.beta1.nimble(TGMRF.GSC), SpatPlus5=sim.beta1(SpatPlus5),
                        SpatPlus10=sim.beta1(SpatPlus10), SpatPlus15=sim.beta1(SpatPlus15),
                        SpatPlus20=sim.beta1(SpatPlus20), SpatPlusP1=sim.beta1(SpatPlusP1),
                        SpatPlusTP1=sim.beta1(SpatPlusTP1), SpatPlusP2=sim.beta1(SpatPlusP2),
                        SpatPlusTP2=sim.beta1(SpatPlusTP2))


colnames(Table.est) <- c("Null", "Spatial", "RSR", "TGRMF1", 
                         "SpatPlus5", "SpatPlus10", "SpatPlus15", "SpatPlus20",
                         "SpatPlusP1", "SpatPlusTP1", "SpatPlusP2", "SpatPlusTP2")

colnames(Table.sim) <- colnames(Table.est) 



Table.sim.est <- rbind(Table.est, Table.sim)
rownames(Table.sim.est) <- c("estimated", "simulated")

Table.sim.est[, 1:12] <- round(Table.sim.est[, 1:12], 4)
Table.sim.est <- data.frame(t(Table.sim.est))
print(Table.sim.est)





############################################################################
# TABLE 8: Empirical 95% coverage probabilities of the true value of beta1 # 
############################################################################

real.beta1 <- 0.2

Coverage.beta1 <- function(Model, true.value){
  CI.L <- Model$summary.fixed$'0.025quant'[2]
  CI.U <- Model$summary.fixed$'0.975quant'[2]
  ifelse(true.value>=CI.L & true.value<=CI.U, 1, 0)
}

Coverage.beta1.nimble <- function(Model, true.value){
  CI.L=Model$quantiles[grep('X1', rownames(Model$quantiles), value=TRUE), 1]                  
  CI.U=Model$quantiles[grep('X1', rownames(Model$quantiles), value=TRUE), 5]
  ifelse(true.value>=CI.L & true.value<=CI.U, 1, 0)
}


Coverage.beta1.95 <- data.frame(Null=mean(unlist(lapply(Null, function(x) Coverage.beta1(x,real.beta1))))*100,
                       Spat=mean(unlist(lapply(Spat, function(x) Coverage.beta1(x,real.beta1))))*100,
                       RSR=mean(unlist(lapply(RSR, function(x) Coverage.beta1(x,real.beta1))))*100,
                       TGMRF1=mean(unlist(lapply(TGMRF.GSC, function(x) Coverage.beta1.nimble(x,real.beta1))))*100,
                                
                       SpatPlus5=mean(unlist(lapply(SpatPlus5, function(x) Coverage.beta1(x,real.beta1))))*100,
                       SpatPlus10=mean(unlist(lapply(SpatPlus10, function(x) Coverage.beta1(x,real.beta1))))*100,
                       SpatPlus15=mean(unlist(lapply(SpatPlus15, function(x) Coverage.beta1(x,real.beta1))))*100,
                       SpatPlus20=mean(unlist(lapply(SpatPlus20, function(x) Coverage.beta1(x,real.beta1))))*100,
                       SpatPlusP1=mean(unlist(lapply(SpatPlusP1, function(x) Coverage.beta1(x,real.beta1))))*100,
                       SpatPlusTP1=mean(unlist(lapply(SpatPlusTP1, function(x) Coverage.beta1(x,real.beta1))))*100,
                       SpatPlusP2=mean(unlist(lapply(SpatPlusP2, function(x) Coverage.beta1(x,real.beta1))))*100,
                       SpatPlusTP2=mean(unlist(lapply(SpatPlusTP2, function(x) Coverage.beta1(x,real.beta1))))*100)   



Coverage.beta1.95 <- as.data.frame(t(Coverage.beta1.95))
colnames(Coverage.beta1.95) <- paste0("cor0.", Subscenario)
print(Coverage.beta1.95)






