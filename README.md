# Evaluating recent methods to overcome spatial confounding
This repository contains the R code to fit the models described in the paper entitled "Evaluating recent methods to overcome spatial confounding"


## Table of contents

- [Data](#Data)
- [Simulated data](#SimulatedData)
- [R code](#R-code)
- [References](#References)


# Data
### Dowry deaths data in Uttar Pradesh in 2001 (Vicente et al., 2020).

The [**Dowry_death_2001.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Dowry_death_2001.Rdata) file contains the following objects:
  - **_Data_**: contains the data set used. It is a dataframe with the following variables:
    - **_dist_**: districts of Uttar Pradesh
    - **_O_**: number of dowry deaths in each district in 2001
    - **_E_**: number of expected cases of each district in 2001
    - **_X1_**: sex ratio covariate (number of females per 1000 males)

  - **_carto_**: SpatialPolygonDataFrame object with the cartography of the 70 districts (year 2001) of Uttar Pradesh
  - **_Q.xi_**: adjacency matrix
  - **_Q.xi.TGMRF_**: adjacency matrix for _tgmrf()_ function of TMGMRF package
	

# R code
R code to fit the spatio-temporal models described in the paper has been included [here](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R).
Only models for the set of hyperprior distributions H1 are shown (to fit the models with H2 and H3 hyperprior distributions slight modifications are required in the code). 
- [icar_models](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/icar_models) and [bym_models](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/bym_models) folders contain the Rscripts with the spatio-temporal models fitted with ICAR and BYM spatial priors using R-INLA, Nimble 1 and Nimble 2. 
- [run](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/run) folder contains the Rscripts to run these models.
- [tables_figures_paper.R](https://github.com/ArantxaUrdangarin/Comparing-R-INLA-and-NIMBLE/blob/main/R/tables_figures_paper.R) contains the necessary functions to reproduce all the figures and tables of Spanish breast cancer mortality data analysis.

Computations were run using R-4.0.3, INLA version 21.02.23 and NIMBLE version 0.11.1.

# Acknowledgements
This work has been supported by Project PID2020-113125RB-I00/ MCIN/ AEI/ 10.13039/501100011033.

![image](https://github.com/spatialstatisticsupna/Comparing-R-INLA-and-NIMBLE/blob/main/micin-aei.jpg)
 
# References
[Urdangarin, A., Goicoa, T. and Ugarte, M.D. (2022). Space-time interactions in Bayesian disease mapping with recent tools: making things easier for practitioners. _Statistical Methods in Medical Research (online first)_. DOI: 10.1177/09622802221079351](https://journals.sagepub.com/doi/10.1177/09622802221079351)
	 

