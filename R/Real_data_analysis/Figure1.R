
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(sp)
library(RColorBrewer)
library(tmap)




##############################
# Load dowry death data 2001 #
##############################

load("../../Data/Dowry_death_2001.Rdata")  # O=observed counts, E=expected counts, X1=standardized sex ratio




#################################
# Figure 1: Sex Ratio covariate #
#################################

carto$X1 <- Data$X1


paleta <- c("#fff7ec", brewer.pal(9,"YlOrRd"))
values <- c(-2.7, seq(-2, 2, 0.5), 2.6)


SexRatio.plot <- tm_shape(carto) +
  tm_polygons(col="X1", palette=paleta, border.alpha=1, title="",
              leyend.show=T, legend.reverse=T, style="fixed", breaks=values, interval.closure="left",
              border.col = "black") +
  tm_layout(legend.outside=T, legend.outside.position="right", legend.frame=F) +
  tm_borders(col="black") + tm_layout(frame=F)
print(SexRatio.plot)
tmap_save(SexRatio.plot, width=7, height=10, file="Figure1.pdf")




