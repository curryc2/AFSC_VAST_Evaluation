# Test_Knot_Number
Script to run comparison of VAST indices with differing knot choices.
(1) Calcualte design-based index.(NOTE: 2001 will be excluded from comparison due to incomplete survey)
(2) Loop through trial knot numbers.
(3) Build folder structure.
(4) Call [VAST](https://github.com/James-Thorson/VAST) in parallel across species, for the sake of efficency.
(5) No plots, only Save.RData. Then use cleanup function call to reduce storage req.
(6) Final plotting section to make multipannel comparing [VAST](https://github.com/James-Thorson/VAST) model-based estimates with differing numbers of knots, with design based estimates, across species. 
(7) Boxplot (ggplot2) comparison of CV's across time, knot number, and species. 

*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)



