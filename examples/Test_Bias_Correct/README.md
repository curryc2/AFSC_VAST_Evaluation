# Test_Bias_Correct
#### R Script: Test-Bias-Correct.R

***
Script to compare [VAST](https://github.com/James-Thorson/VAST) model-based estimates with and without bias correction. To determine if the scale of indices are more sensitive to: a) spatial complexity (number of knots), or b) bias correction.

***

1.  Loop through bias.correction =  TRUE/FALSE
2.  Loop through trial knot numbers: 100, 500, 1000.
3.  Plot comparing both to design-based estimates. 
    
*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)

***
#### Issues
1.  None.

***
#### Output Objects

Object                 | Description
-----------------------|-----------------------------------------------
Name                   | $object - Description
