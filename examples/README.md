# AFSC_VAST_Evaluation: Example Scripts

*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)


***
## Run VAST Single Species.r
Example of calls to helper functions to:
1.  Read in data
2.  Construct VAST input
3.  Compile and build VAST (TMB) model
4.  Generate figures and model outputs

***
## Test_Parallel
Example of how to create wrapper function to run [VAST](https://github.com/James-Thorson/VAST) models in parallel across species.

***
## Test_Obs_Model
Evaluation of differences in [VAST](https://github.com/James-Thorson/VAST) model results for example species, with different observation models.

***
## Test_DesignBased_Estimator
Simple comparison of design-based and [VAST](https://github.com/James-Thorson/VAST) model-based indices across years and species. 

***
## Test_Knot_Number
Script to run comparison of [VAST](https://github.com/James-Thorson/VAST) indices with differing levels of spatial complexity, as specified by the number of "knots".
* Note: Directory for this example has two scripts "Test-Knot-Number.r" and "Test-Knot-Number-Updated.r"
    + Updated version should be used.
    + Test-Knot-Number.r -- creates a large file structure with much saved output. This was inefficient from a storage perspective.
    + Test-Knot-Number-Updated.r -- Saves VAST indices at .RData object, which can be readily compared to design-based estimates. Saves space...

***
## Test_Autoregressive
Compares influence of different autoregressive structures for [VAST](https://github.com/James-Thorson/VAST) model parameters over time.
    
***
## Test_DeltaModel
Script to generate stratified Delta model biomass indices, and compare with design-based estimates.

***

