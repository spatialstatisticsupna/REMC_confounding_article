# Evaluating recent methods to overcome spatial confounding

This repository contains the R code to implement the methods described in the paper entitled "Evaluating recent methods to overcome spatial confounding" [(Urdangarin et al., 2023)](https://link.springer.com/article/10.1007/s13163-022-00449-8) as well as the R code to create the figures and tables presented in the paper.

## Table of contents

-   [Data](#Data)
-   [Simulated data](#Simulated-data)
-   [R code](#R-code)
-   [References](#References)

# Data

### Dowry deaths data in Uttar Pradesh in 2001 [(Vicente et al., 2020)](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssa.12545)

The [**Dowry_death_2001.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Dowry_death_2001.Rdata) file contains the following objects:

-   ***Data***: contains the data set used. It is a `data.frame` object with the following variables:

    -   ***dist***: names of the districts of Uttar Pradesh

    -   ***ID.area***: numeric identifiers of districts

    -   ***O***: number of dowry deaths in each district in 2001

    -   ***E***: number of expected cases of each district in 2001

    -   ***X1***: standardized sex ratio covariate (number of females per 1000 males)

-   ***carto***: `SpatialPolygonDataFrame` object with the cartography of the 70 districts (year 2001) of Uttar Pradesh

-   ***Q.xi***: spatial adjacency matrix

### Stomach cancer incidence data in Slovenia during the period 1995-2001 (Zadnik and Reich, 2006)

The [**Slovenia_stomach_cancer**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Slovenia_stomach_cancer.Rdata) file contains the following objects:

-   ***Data***: contains the data set used. It is a `data.frame` object with the following variables:

    -   ***ID.area***: numeric district identifiers

    -   ***O***: number of stomach cancer cases in each area during 1995-2001

    -   ***E***: number of expected cases in each area during 1995-2001

    -   ***X***: standardized socioeconomic indicator

-   ***coord***: a matrix that contains the coordinates of the 192 areas of Slovenia

-   ***Q.xi***: spatial adjacency matrix

*Slovenia data set is available from the package RASCO of R <https://github.com/DouglasMesquita/RASCO>. This dataset is also available from the [web page of James Hogdes](https://www.biostat.umn.edu/~hodges/RPLMBook/Datasets/Datasets.html)*.

### Lip cancer incidence data in Scotland during 1975-1980 [(Breslow and Clayton, 1993)](https://www.jstor.org/stable/2290687?origin=crossref#metadata_info_tab_contents)

The [**Scotland_lip_cancer.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Scotland_lip_cancer.Rdata) file contains the following objects:

-   ***Data***: contains the data set used. It is a dataframe with the following variables:

    -   ***ID.area***: numeric district identifiers

    -   ***O***: number of lip cancer cases in each area during 1975-1980

    -   ***E***: number of expected cases in each area during 1975-1980

    -   ***AFF***: standardized covariate indicating the proportion of the population engaged in agriculture, fishing, or forestry

-   ***carto***: SpatialPolygonDataFrame object with the cartography of the 56 districts of Scotland

-   ***Q.xi***: spatial adjacency matrix

# Simulated data

[Simulated_data](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/Simulated_data) folder contains a total of 18 .Rdata files (one file for each scenario and subscenario) used in Simulation Study 1 and Simulation Study 2. Each .Rdata file contains the same objects as [**Dowry_death_2001.Rdata**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/blob/main/Data/Dowry_death_2001.Rdata) (***Data***, ***carto***, ***Q.xi***) but a simulated covariate ***X2*** is added to ***Data***. Moreover, each .Rdata contains the following objects as well:

-   ***log.risk***: a vector that contains the simulated log risks
-   ***simu.O***: a list with 100 simulated counts data sets

The code used to simulate the data is available in [SimuStudy1_simulate_data.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_simulate_data.R) and [SimuStudy2_simulate_data.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy2_simulate_data.R).

# R code

R code to implement the procedures to alleviate spatial confounding described in the paper and to reproduce the tables and figures of the paper has been included.

-   [**R/Real_data_analysis**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/R/Real_data_analysis) folder contains the R code used in the real data analysis.

    The main file to fit the null, spatial, RSR and spatial+ models is [Null_Spatial_RSR_SpatPlus_models.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Real_data_analysis/Null_Spatial_RSR_SpatPlus_models.R). Before running the models, the `dataset` argument (one of either "Dowry", "Slovenia" or "Scotland") must be defined at the top of the code.

    -   [Figure1.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Real_data_analysis/Figure1.R): R script to reproduce Figure 1 of the paper.
    -   [Covariate_model_eigenvectors.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Real_data_analysis/Covariate_model_eigenvectors.R): R script to fit the covariate model based on the eigenvectors of the spatial precision matrix to remove the spatial dependence from the covariate before fitting the spatial+ model.
    -   [Covariate_model_Psplines.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Real_data_analysis/Covariate_model_Psplines.R): R script to fit the covariate model based on P-splines to remove the spatial dependence from the covariate before fitting the spatial+ model.
    -   [Covariate_model_TPsplines.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Real_data_analysis/Covariate_model_TPsplines.R): R script to fit the covariate model based on thin plate splines to remove the spatial dependence from the covariate before fitting the spatial+ model.

-   [**R/Simulation_Study_1**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/R/Simulation_study_1) folder contains the R code used in Simulation Study 1.

    Before running the models, the arguments `Scenario`(1, 2 or 3) and `Subscenario` (cor=80, 50 or 20) must be defined at the top of the code.

    -   [Figure2.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/Figure2.R): R code to reproduce Figure 2 of the paper.
    -   [SimuStudy1_Null.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_Null.R): R script to fit the null model to the 100 simulated datasets.
    -   [SimuStudy1_Spatial.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_Spatial.R): R script to fit the spatial model to the 100 simulated datasets.
    -   [SimuStudy1_RSR.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_RSR.R): R script to fit the RSR model to the 100 simulated datasets.
    -   [SimuStudy1_SpatialPlus_eigenvectors.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatialPlus_eigenvectors.R): R script to fit SpatPlus5, SpatPlus10, SpatPlus15 and SpatPlus20 models to the 100 simulated datasets.
    -   [SimuStudy1_SpatPlusP1.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatPlusP1.R): R script to fit SpatPlusP1 model to the 100 simulated datasets.
    -   [SimuStudy1_SpatPlusTP1.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatPlusTP1.R): R script to fit SpatPlusTP1 model to the 100 simulated datasets.
    -   [SimuStudy1_SpatPlusP2.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatPlusP2.R): R script to fit SpatPlusP2 model to the 100 simulated datasets.
    -   [SimuStudy1_SpatPlusTP2.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatPlusTP2.R): R script to fit SpatPlusTP2 model to the 100 simulated datasets.
    -   [Tables_6_7_8.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/Tables_6_7_8.R): R code to reproduce Table 6, 7 and 8 of the paper for each scenario and subscenario.
    -   [Tables_supplementary_A1_A2_A3_A4.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/Tables_supplementary_A1_A2_A3_A4.R): R code to reproduce Tables A1, A2 and A3 of the supplementary material for each scenario and subscenario.

-   [**R/Simulation_Study_2**](https://github.com/spatialstatisticsupna/Spatial_confounding_article/tree/main/R/Simulation_study_2) folder contains the R code used in Simulation Study 2.

    Before running the models, the arguments `Scenario`(1, 2 or 3) and `Subscenario` (cor=80, 50 or 20) must be defined at the top of the code.

    -   [SimuStudy2_Null.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/SimuStudy2_Null.R): R script to fit the null model to the 100 simulated datasets.
    -   [SimuStudy2_Spatial.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/SimuStudy2_Spatial.R): R script to fit the spatial model to the 100 simulated datasets.
    -   [SimuStudy2_RSR.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/SimuStudy2_RSR.R): R script to fit the RSR model to the 100 simulated datasets.
    -   [SimuStudy2_SpatialPlus_eigenvectors.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/SimuStudy2_SpatialPlus_eigenvectors.R): R script to fit SpatPlus5, SpatPlus10, SpatPlus15 and SpatPlus20 models to the 100 simulated datasets.
    -   [SimuStudy2_SpatPlusP1.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/SimuStudy2_SpatPlusP1.R): R script to fit SpatPlusP1 model to the 100 simulated datasets.
    -   [SimuStudy2_SpatPlusTP1.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/SimuStudy2_SpatPlusTP1.R): R script to fit SpatPlusTP1 model to the 100 simulated datasets.
    -   [SimuStudy2_SpatPlusP2.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/SimuStudy2_SpatPlusP2.R): R script to fit SpatPlusP2 model to the 100 simulated datasets.
    -   [SimuStudy2_SpatPlusTP2.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/SimuStudy2_SpatPlusTP2.R): R script to fit SpatPlusTP2 model to the 100 simulated datasets.
    -   [Table_9.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/Table_9.R): R code to reproduce Table 9 of the paper for each scenario and subscenario.
    -   [Figures_supplementary_A3_A4_A5.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_2/Figures_supplementary_A3_A4_A5.R): R code to reproduce Figures A3, A4 and A4 of the supplementary material for each scenario and subscenario.

Computations were run using R-4.0.4, INLA version 21.02.23, mgcv version 1.8-40.

# Acknowledgements

This work has been supported by Project PID2020-113125RB-I00/ MCIN/ AEI/ 10.13039/501100011033.

![image](https://github.com/spatialstatisticsupna/Comparing-R-INLA-and-NIMBLE/blob/main/micin-aei.jpg)

# References

[Urdangarin, A., Goicoa, T. and Ugarte, M.D. (2023). Evaluating recent methods to overcome spatial confounding. *Revista Matemática Complutense* **36**, 333-360. DOI: 10.1007/s13163-022-00449-8.](https://link.springer.com/article/10.1007/s13163-022-00449-8)
