# AFSC_VAST_Evaluation functions
=============
Contains a series of functions necessary for evaluation.


Helper functions
=============
1.  load-RACE-data.r - Loads, merges, adds zero catch values for RACE bottom trawl survey data.
2.  create-Data-Geostat.r - Formats input for [VAST](https://github.com/James-Thorson/VAST) as Data_Geostat data frame,
                              after calling load_RACE_data() function. 
3.  create-VAST-input.r - Creates the [VAST](https://github.com/James-Thorson/VAST) TMB input object: TmbData. 
                            TmbData contains VAST model control parametes, spatial grid configuration, and model input data from
                            [RACE](https://www.afsc.noaa.gov/RACE/groundfish/bottom%20trawl%20surveys.php) surveys.