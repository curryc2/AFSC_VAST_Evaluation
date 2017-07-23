# AFSC_VAST_Evaluation functions
## Contains a series of functions necessary for evaluation.

***

### Helper Functions
Function Name                   | Details                                       | Instance
--------------------------------|-----------------------------------------------|----------------------
load-RACE-data.r                | Loads, merges, adds zero catch values for RACE bottom trawl survey data. | Multispecies
create-Data-Geostat.r           | Formats input for [VAST](https://github.com/James-Thorson/VAST) as Data_Geostat data frame, after calling load_RACE_data() function. | Multispecies 
create-VAST-input.r             | Creates the [VAST](https://github.com/James-Thorson/VAST) TMB input object: TmbData. TmbData contains VAST model control parametes, spatial grid configuration, and model input data from [RACE](https://www.afsc.noaa.gov/RACE/groundfish/bottom%20trawl%20surveys.php) surveys. | Multispecies
plot-VAST-output.r              | Wrapper function that generates [VAST](https://github.com/James-Thorson/VAST) diagnostic and prediction plots. | Single
cleanup-VAST-file.r             | Simple function to remove unnecessary VAST TMB model files after estimation. | Single
calc-design-based-index.R       | Calculates design-based biomass estimate from [RACE](https://www.afsc.noaa.gov/RACE/groundfish/bottom%20trawl%20surveys.php) surveys, along with Variance, SD, and CV. | Single
get-VAST-index.r                | Simple function to retreive [VAST](https://github.com/James-Thorson/VAST) index value and variance estimates. NOTE: This function only a simplified version of the [PlotIndex_Fn.R](https://github.com/nwfsc-assess/geostatistical_delta-GLMM/blob/master/R/PlotIndex_Fn.R) function in the [SpatialDeltaGLMM](https://github.com/nwfsc-assess/geostatistical_delta-GLMM) package. | Single
run-RE-model.R                  | Function to call random effects model for area apportionment. Builds input (.dat) file, copy and compiles ADMB-RE model, and executes.  | Single         