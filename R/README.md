# AFSC_VAST_Evaluation functions
## Contains a series of functions necessary for evaluation.

***

### Helper Functions
Function Name                   | Details
--------------------------------|-----------------------------------------------
load-RACE-data.r                | Loads, merges, adds zero catch values for RACE bottom trawl survey data.
create-Data-Geostat.r           | Formats input for [VAST](https://github.com/James-Thorson/VAST) as Data_Geostat data frame, after calling load_RACE_data() function. 
create-VAST-input.r             | Creates the [VAST](https://github.com/James-Thorson/VAST) TMB input object: TmbData. TmbData contains VAST model control parametes, spatial grid configuration, and model input data from [RACE](https://www.afsc.noaa.gov/RACE/groundfish/bottom%20trawl%20surveys.php) surveys.
plot-VAST-output.r              | Wrapper function that generates [VAST](https://github.com/James-Thorson/VAST) diagnostic and prediction plots.
cleanup-VAST-file.r             | Simple function to remove unnecessary VAST TMB model files after estimation.
calc-design-based-index.r       | Calculates design-based biomass estimate from [RACE](https://www.afsc.noaa.gov/RACE/groundfish/bottom%20trawl%20surveys.php) surveys, along with Variance, SD, and CV. 