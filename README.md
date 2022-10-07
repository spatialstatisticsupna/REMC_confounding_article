# Evaluating recent methods to overcome spatial confounding
This repository contains the R code to fit the models described in the paper entitled "Evaluating recent methods to overcome spatial confounding"


## Table of contents

- [Data](#Data)
- [Simulated data](#SimulatedData)
- [R code](#R-code)
- [References](#References)


# Data
### Dowry deaths data in Uttar Pradesh in 2001

The [**Dowry_death_2001.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Dowry_death_2001.Rdata) file contains the following objects:
  - **_Data_**: contains the data set used. It is a dataframe with the following variables:
    - **_dist_**: names of the districts of Uttar Pradesh
    - **_ID.area_**: numeric vector of district identifiers
    - **_O_**: number of dowry deaths in each district in 2001
    - **_E_**: number of expected cases of each district in 2001
    - **_X1_**: standardized sex ratio covariate (number of females per 1000 males)

  - **_carto_**: SpatialPolygonDataFrame object with the cartography of the 70 districts (year 2001) of Uttar Pradesh
  - **_Q.xi_**: spatial adjacency matrix
  - **_Q.xi.TGMRF_**: spatial adjacency matrix for _tgmrf()_ function of TMGMRF package


### Stomach cancer incidence data in Slovenia during the period 1995-2001

The [**Slovenia_stomach_cancer**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Slovenia_stomach_cancer.Rdata) file contains the following objects:
  - **_Data_**: contains the data set used. It is a dataframe with the following variables:
    - **_ID.area_**: numeric vector of area identifiers
    - **_O_**: number of stomach cancer cases in each area during 1995-2001
    - **_E_**: number of expected cases in each area during 1995-2001
    - **_X_**: standardized socioeconomic indicator
    
  - **_coord_**: a matrix that contains the coordinates of the 192 areas of Slovenia
  - **_Q.xi_**: spatial adjacency matrix
  - **_Q.xi.TGMRF_**: spatial adjacency matrix for _tgmrf()_ function of TMGMRF package
  
_**Note:** Slovenia data is obtained from the package RASCO of R [https://github.com/DouglasMesquita/RASCO](https://github.com/DouglasMesquita/RASCO)_


### Lip cancer incidence data in Scotland during 1975-1980

The [**Dowry_death_2001.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Dowry_death_2001.Rdata) file contains the following objects:
  - **_Data_**: contains the data set used. It is a dataframe with the following variables:
    - **_dist_**: names of the districts of Uttar Pradesh
    - **_ID.area_**: numeric vector of district identifiers
    - **_O_**: number of dowry deaths in each district in 2001
    - **_E_**: number of expected cases of each district in 2001
    - **_X1_**: sex ratio covariate (number of females per 1000 males)

  - **_carto_**: SpatialPolygonDataFrame object with the cartography of the 70 districts (year 2001) of Uttar Pradesh
  - **_Q.xi_**: spatial adjacency matrix
  - **_Q.xi.TGMRF_**: spatial adjacency matrix for _tgmrf()_ function of TMGMRF package





# Simulated data
[Simulated_data](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/Simulated_data) folder contains 18 .Rdata files (one file for each Scenario) used in Simulation Study 1 and Simulation Study 2. Each .Rdata file contains the same objects as [**Dowry_death_2001.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Dowry_death_2001.Rdata) but unobserved covariate **_X2_** is added to **_Data_**. Moreover the .Rdata contains the following objects:

- **_log.risk_**: a vector that contains the simulated log risks
- **_simu.O_**: a list with 100 simulated counts data sets


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

	 

