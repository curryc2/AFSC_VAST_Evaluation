# Test_TrawlDuration_Covariate
#### R Script: Test-TrawlDuration-Covariate.R

***
Script to compare model-based [VAST](https://github.com/James-Thorson/VAST) biomass estimates when trawl duration (minutes) is included as a covariate. Given area swept effort accounts for differential trawl distrance, 

***

1.  Implement model with trawl duration as a "catchability" covariate:
    + Need to update or create new function to locate trawl duration.
    + Need to add trawl duration to input data frame.
2.  Run for the range of GOA and AI species and compare:
    + Scale and trend of indices
    + Precision (CV) in indices
    
*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)

***
#### Issues
1.  None.

***
#### Output Objects

Object                 | Description
-----------------------|-----------------------------------------------
Name                   | $object - Description
