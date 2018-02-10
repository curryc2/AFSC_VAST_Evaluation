# AFSC_VAST_Evaluation
Repository for evaluation of VAST package for single and multispecies spatial delta-generalized linear mixed models (delta-GLMM). The purpose of this analysis is to evaluate differences in indices and uncertainty estimates generated by VAST and design-based estimators, for a range of species in the Bering Sea and Gulf of Alaska.

## Analyses
### Individual analyses are contained within the [examples]("https://github.com/curryc2/AFSC_VAST_Evaluation/tree/develop/examples") folder, and inclue:

Name                        | Description
----------------------------|----------------------------------------------------
[Test_Parallel]("https://github.com/curryc2/AFSC_VAST_Evaluation/tree/develop/examples/Test_Parallel") | Example of how to create wrapper function to run VAST models in parallel across species.
[Test_Obs_Model]("https://github.com/curryc2/AFSC_VAST_Evaluation/tree/develop/examples/Test_Obs_Model") | Evaluation of differences in VAST model results for example species, with different observation models.
[Test_DesignBased_Estimator]("https://github.com/curryc2/AFSC_VAST_Evaluation/tree/develop/examples/Test_DesignBased_Estimator") | Simple comparison of design-based and VAST model-based indices across years and species.
[Test_Knot_Number]("https://github.com/curryc2/AFSC_VAST_Evaluation/tree/develop/examples/Test_Knot_Number") | Script to run comparison of VAST indices with differing levels of spatial complexity, as specified by the number of "knots".
[Test_Autoregressive]("https://github.com/curryc2/AFSC_VAST_Evaluation/tree/develop/examples/Test_Autoregressive") | Compares influence of different autoregressive structures for VAST model parameters over time.
[Test_DeltaModel]("https://github.com/curryc2/AFSC_VAST_Evaluation/tree/develop/examples/Test_DeltaModel") | Script to generate stratified Delta model biomass indices, and compare with design-based estimates.

***

## VAST Package Citation
#### Please cite if using the software
* Thorson, J.T., and Barnett, L.A.K. In press. Comparing estimates of abundance trends and distribution shifts using single- and multispecies models of fishes and biogenic habitat. ICES J. Mar. Sci. [URL](https://academic.oup.com/icesjms/article/74/5/1311/2907795)

***

## [VAST Package Website](https://github.com/James-Thorson/VAST)

***

## Data from NOAA: Resource Assessment and Conservation Engineering Division
Please visit the following website for more information: 
[RACE](https://www.afsc.noaa.gov/RACE/groundfish/bottom%20trawl%20surveys.php)