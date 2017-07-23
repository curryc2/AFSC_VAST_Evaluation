# Test_DesignBased_Estimator
#### R Script: Test-DesignBased-Estimator.r

***
##### Purpose: Compare index standardization methods. 
  1.  Generate output object with design-based estimates for each species x survey combination.
    + Output file with design-based estimtes: "Design Based Estimates.xlsx"
  2.  Compare model-based [VAST](https://github.com/James-Thorson/VAST) with design-based indices.
    + Bias Correction is used and 250 knots assumed (knot assumption tested in alternative [example](https://github.com/curryc2/AFSC_VAST_Evaluation/tree/develop/examples/Test_Knot_Number)).

*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)
