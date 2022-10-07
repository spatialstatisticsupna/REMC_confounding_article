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
    - **_ID.area_**: numeric vector of district identifiers
    - **_O_**: number of stomach cancer cases in each area during 1995-2001
    - **_E_**: number of expected cases in each area during 1995-2001
    - **_X_**: standardized socioeconomic indicator
    
  - **_coord_**: a matrix that contains the coordinates of the 192 areas of Slovenia
  - **_Q.xi_**: spatial adjacency matrix
  - **_Q.xi.TGMRF_**: spatial adjacency matrix for _tgmrf()_ function of TMGMRF package
  
_**Note:** Slovenia data is obtained from the package RASCO of R [https://github.com/DouglasMesquita/RASCO](https://github.com/DouglasMesquita/RASCO)_


### Lip cancer incidence data in Scotland during 1975-1980

The [**Scotland_lip_cancer.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Scotland_lip_cancer.Rdata) file contains the following objects:
  - **_Data_**: contains the data set used. It is a dataframe with the following variables:
    - **_ID.area_**: numeric vector of district identifiers
    - **_O_**: number of lip cancer cases in each area during 1975-1980
    - **_E_**: number of expected cases in each area during 1975-1980
    - **_AFF_**: standardized covariate indicating the proportion of the population engaged in agriculture, fishing, or forestry

  - **_carto_**: SpatialPolygonDataFrame object with the cartography of the 56 districts of Scotland
  - **_Q.xi_**: spatial adjacency matrix
  - **_Q.xi.TGMRF_**: spatial adjacency matrix for _tgmrf()_ function of TMGMRF package



# Simulated data
[Simulated_data](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/Simulated_data) folder contains a total of 18 .Rdata files (one file for each scenario and subscenario) used in Simulation Study 1 and Simulation Study 2. Each .Rdata file contains the same objects as [**Dowry_death_2001.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Dowry_death_2001.Rdata) (**_Data_**, **_carto_**, **_Q.xi_**, **_Q.xi.TGMRF_**) but simulated covariate **_X2_** is added to **_Data_**. Moreover, each .Rdata contains the following objects as well:

- **_log.risk_**: a vector that contains the simulated log risks
- **_simu.O_**: a list with 100 simulated counts data sets


# R code

R code to implement the procedures to alleviate spatial confounding described in the paper and reproduce the tables and figures of the paper has been included. 

- [R/Real_data_analysis](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/R/Real_data_analysis) folder contains the R code used in the real data analysis. The R code for dowry death data is presented, slight modifications included in the code as a comments should be done to fit the models to Slovenia and Scotland datasets.
 
- [R/Simulation_Study_1](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/R/Simulation_study_1) folder contains the R code used in Simulation Study 1. Before running the models, the options _Scenario_ (scenario 1, 2 or 3) and _Subscenario_ (cor=0.80,0.50 or 0.20) must be defined. It is important to run [Prepare_results_TGMRF.R](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/R/Simulation_study_1/Prepare_results_TGMRF.R) script before running the scripts to obtain the tables and figures. 

- [R/Simulation_Study_2](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/R/Simulation_study_2) folder contains the R code used in Simulation Study 2. Before running the models, the options _Scenario_ (scenario 1, 2 or 3) and _Subscenario_ (cor=0.80,0.50 or 0.20) must be defined. It is important to run [Prepare_results_TGMRF.R](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/R/Simulation_study_2/Prepare_results_TGMRF.R) script before running the scripts to obtain the tables and figures. 

Computations were run using R-4.0.3, INLA version 21.02.23 and NIMBLE version 0.11.1.

# Acknowledgements
This work has been supported by Project PID2020-113125RB-I00/ MCIN/ AEI/ 10.13039/501100011033.

![image](https://github.com/spatialstatisticsupna/Comparing-R-INLA-and-NIMBLE/blob/main/micin-aei.jpg)
 
# References

	 

