

rm(list=ls())
setwd("")


n.sim <- 100 


# Load packages
library(ggplot2)
library(latex2exp)




####################
# Define functions #
####################

beta1 <- function(x){
  data.frame(alpha.mean=x$summary.fixed[2,1], 
             alpha.sd=x$summary.fixed[2,2])
}

beta1.nimble <- function(x){
  data.frame(alpha.mean=x$statistics[grep('X1', rownames(x$statistics), value=TRUE), 1],
             alpha.sd=x$statistics[grep('X1', rownames(x$statistics), value=TRUE), 2])
}




################################
# Load and prepare the results #
################################

## Subscenario 1 (cor=0.80) ##
Scenario <- 1  #2,3
Subscenario <- 80

load(paste0("Results_SimuStudy1/Scenario", Scenario, "_all_models_cor", Subscenario, ".Rdata"))

Null.cor080 <- do.call(rbind, lapply(Null, beta1))[, 1]
Spat.cor080 <- do.call(rbind, lapply(Spat, beta1))[, 1]
RSR.cor080 <- do.call(rbind, lapply(RSR, beta1))[, 1]
TGMRF.cor080 <- do.call(rbind, lapply(TGMRF.GSC, beta1.nimble))[, 1]


SpatPlus5.cor080 <- do.call(rbind, lapply(SpatPlus5, beta1))[, 1]
SpatPlus10.cor080 <- do.call(rbind, lapply(SpatPlus10, beta1))[, 1]
SpatPlus15.cor080 <- do.call(rbind, lapply(SpatPlus15, beta1))[, 1]
SpatPlus20.cor080 <- do.call(rbind, lapply(SpatPlus20, beta1))[, 1]


SpatPlusP1.cor080 <- do.call(rbind, lapply(SpatPlusP1, beta1))[, 1]
SpatPlusTP1.cor080 <- do.call(rbind, lapply(SpatPlusTP1, beta1))[, 1]
SpatPlusP2.cor080 <- do.call(rbind, lapply(SpatPlusP2, beta1))[, 1]
SpatPlusTP2.cor080 <- do.call(rbind, lapply(SpatPlusTP2, beta1))[, 1]



## Subscenario 2 (cor=0.50) ##
Scenario <- 1  #2,3
Subscenario <- 50

load(paste0("Results_SimuStudy1/Scenario", Scenario, "_all_models_cor", Subscenario, ".Rdata"))

Null.cor050 <- do.call(rbind, lapply(Null, beta1))[, 1]
Spat.cor050 <- do.call(rbind, lapply(Spat, beta1))[, 1]
RSR.cor050 <- do.call(rbind, lapply(RSR, beta1))[, 1]
TGMRF.cor050 <- do.call(rbind, lapply(TGMRF.GSC, beta1.nimble))[, 1]


SpatPlus5.cor050 <- do.call(rbind, lapply(SpatPlus5, beta1))[, 1]
SpatPlus10.cor050 <- do.call(rbind, lapply(SpatPlus10, beta1))[, 1]
SpatPlus15.cor050 <- do.call(rbind, lapply(SpatPlus15, beta1))[, 1]
SpatPlus20.cor050 <- do.call(rbind, lapply(SpatPlus20, beta1))[, 1]


SpatPlusP1.cor050 <- do.call(rbind, lapply(SpatPlusP1, beta1))[, 1]
SpatPlusTP1.cor050 <- do.call(rbind, lapply(SpatPlusTP1, beta1))[, 1]
SpatPlusP2.cor050 <- do.call(rbind, lapply(SpatPlusP2, beta1))[, 1]
SpatPlusTP2.cor050 <- do.call(rbind, lapply(SpatPlusTP2, beta1))[, 1]



## Subscenario 3 (cor=0.20) ##
Scenario <- 1  #2,3
Subscenario <- 20

load(paste0("Results_SimuStudy1/Scenario", Scenario, "_all_models_cor", Subscenario, ".Rdata"))

Null.cor020 <- do.call(rbind, lapply(Null, beta1))[, 1]
Spat.cor020 <- do.call(rbind, lapply(Spat, beta1))[, 1]
RSR.cor020 <- do.call(rbind, lapply(RSR, beta1))[, 1]
TGMRF.cor020 <- do.call(rbind, lapply(TGMRF.GSC, beta1.nimble))[, 1]


SpatPlus5.cor020 <- do.call(rbind, lapply(SpatPlus5, beta1))[, 1]
SpatPlus10.cor020 <- do.call(rbind, lapply(SpatPlus10, beta1))[, 1]
SpatPlus15.cor020 <- do.call(rbind, lapply(SpatPlus15, beta1))[, 1]
SpatPlus20.cor020 <- do.call(rbind, lapply(SpatPlus20, beta1))[, 1]


SpatPlusP1.cor020 <- do.call(rbind, lapply(SpatPlusP1, beta1))[, 1]
SpatPlusTP1.cor020 <- do.call(rbind, lapply(SpatPlusTP1, beta1))[, 1]
SpatPlusP2.cor020 <- do.call(rbind, lapply(SpatPlusP2, beta1))[, 1]
SpatPlusTP2.cor020 <- do.call(rbind, lapply(SpatPlusTP2, beta1))[, 1]



## Prepare a dataframe for the boxplots ##

model <- c(rep("Null", 100), rep("Spatial", 100), rep("RSR", 100), rep("TGMRF1", 100), 
           rep("SpatPlus5", 100), rep("SpatPlus10", 100), rep("SpatPlus15", 100), rep("SpatPlus20", 100),
           rep("SpatPlusP1", 100), rep("SpatPlusTP1", 100), rep("SpatPlusP2", 100), rep("SpatPlusTP2", 100))


data.cor080 <- data.frame(cor=rep("cor=0.80", 1200), model=model, 
                          beta1=c(Null.cor080, Spat.cor080, RSR.cor080, TGMRF.cor080,
                                  SpatPlus5.cor080, SpatPlus10.cor080, SpatPlus15.cor080, SpatPlus20.cor080,
                                  SpatPlusP1.cor080, SpatPlusTP1.cor080, SpatPlusP2.cor080, SpatPlusTP2.cor080))

data.cor050 <- data.frame(cor=rep("cor=0.50", 1200), model=model, 
                          beta1=c(Null.cor050, Spat.cor050, RSR.cor050, TGMRF.cor050,
                                  SpatPlus5.cor050, SpatPlus10.cor050, SpatPlus15.cor050, SpatPlus20.cor050,
                                  SpatPlusP1.cor050, SpatPlusTP1.cor050, SpatPlusP2.cor050, SpatPlusTP2.cor050))

data.cor020 <- data.frame(cor=rep("cor=0.20", 1200), model=model, 
                          beta1=c(Null.cor020, Spat.cor020, RSR.cor020, TGMRF.cor020,
                                  SpatPlus5.cor020, SpatPlus10.cor020, SpatPlus15.cor020, SpatPlus20.cor020,
                                  SpatPlusP1.cor020, SpatPlusTP1.cor020, SpatPlusP2.cor020, SpatPlusTP2.cor020))

data <- rbind(data.cor080, data.cor050, data.cor020)
str(data)
data$model <- factor(data$model, levels=c("Null", "Spatial", "RSR", "TGMRF1",
                                          "SpatPlus5", "SpatPlus10", "SpatPlus15", "SpatPlus20",
                                          "SpatPlusP1", "SpatPlusTP1", "SpatPlusP2", "SpatPlusTP2"))
data$cor <- factor(data$cor, levels=c('cor=0.80', 'cor=0.50', 'cor=0.20'))




#############################################################
# FIGURE 3 (Figures A.1, A.2 of the Supplementary material) #
#############################################################

# pdf(file="Figure3.pdf", width=14, height=12)
p <- ggplot(data, aes(x = model, y = beta1, color=cor)) + 
  geom_boxplot(position = position_dodge(0.9)) +
  geom_hline(yintercept=0.2, color="#e31a1c", linetype="dashed") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 0.9),
        text = element_text(size = 20),
        panel.spacing = unit(1.3, "lines")) + 
  xlab("") + ylab(TeX(r'($\beta_{1}$)'))  +
  scale_y_continuous(limits=c(-0.2, 0.85), breaks=c(-0.2, 0, 0.2, 0.4, 0.6, 0.8)) 
# Split in vertical direction
p + facet_grid(rows = vars(cor)) + scale_color_manual(values=c("#253494", "#225ea8", "#1d91c0")) 
# dev.off()







