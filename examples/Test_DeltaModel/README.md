# Test_DeltaModel
#### R Script: Test-DeltaModel.R

***
Script to generate stratified Delta model biomass indices, and compare with current design-based and [VAST](https://github.com/James-Thorson/VAST) model-based estimates.

***

1.  Implement method for running stratified delta-GLMM. Options:
    + Write own function
    + Use VAST with spatial and spatio-temporal components turned off and strata defined
    + Use the Bayesian West Coast delta-GLMM tool already built
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
