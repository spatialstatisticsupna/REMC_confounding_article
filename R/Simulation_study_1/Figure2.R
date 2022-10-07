
rm(list=ls())
setwd("")


# Load packages
library(RColorBrewer)
library(rgdal)
library(tmap)




#############################
# Load and prepare the data #
#############################

## Subscenario 1 ##
load("Simulated_data/Simu1_data_Scenario2_cor080.Rdata")
Scenario2 <- data.2001
load("Simulated_data/Simu1_data_Scenario3_cor080.Rdata")
Scenario3 <- data.2001


## Subscenario 2 ##
# load("Simulated_data/Simu1_data_Scenario2_cor050.Rdata")
# Scenario2 <- data.2001
# load("Simulated_data/Simu1_data_Scenario3_cor050.Rdata")
# Scenario3 <- data.2001


## Subscenario 3 ##
# load("Simulated_data/Simu1_data_Scenario2_cor020.Rdata")
# Scenario2 <- data.2001
# load("Simulated_data/Simu1_data_Scenario3_cor020.Rdata")
# Scenario3 <- data.2001


rm(data.2001, log.risk, simu.O)
Data.plot <- data.frame(ID_area=1:70, X1=Scenario2$X1, X2=Scenario2$X2, 
                        ICAR=Scenario2$spat, Psplines=Scenario3$spat)


## Load the  cartography ##
load("Real_data_analysis/Data/Dowry_death_2001.Rdata")
rm(Data, Q.xi, Q.xi.TGMRF)




############
# FIGURE 2 #
############

# Merge cartography and the data
Carto <- merge(carto, Data.plot, by="ID_area")



paleta <- c("#fff7ec", brewer.pal(9,"YlOrRd"))
values <- c(-2.7, seq(-2, 2, 0.5), 2.6)



tmap_mode("plot")
plot.map <- tm_shape(Carto) +
  tm_polygons(col=c("X1", "X2", "ICAR", "Psplines"), palette=paleta, title="", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left", border.col="black") +
  tm_layout(main.title="", main.title.position="center", legend.text.size=1.5,
            panel.labels=c("X1","X2", "ICAR", "P-splines"), panel.label.bg.color="lightskyblue",
            legend.outside=T, legend.outside.position="right", legend.frame=F) +
  tm_borders(col="black") + 
  tm_facets(ncol=4, nrow=1)
print(plot.map)
tmap_save(plot.map, file="Simulation_study_1/Results_SimuStudy1/Figure2A.pdf", width=10, height=3)
# tmap_save(plot.map, file="Simulation_study_1/Results_SimuStudy1/Figure2B.pdf", width=10, height=3)
# tmap_save(plot.map, file="Simulation_study_1/Results_SimuStudy1/Figure2C.pdf", width=10, height=3)



