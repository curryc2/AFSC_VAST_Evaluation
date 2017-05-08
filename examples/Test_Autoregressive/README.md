# Test_Autoregressive
#### R Script: Test-Autoregressive.r

***
Script to compare influence of different autoregressive structures for [VAST](https://github.com/James-Thorson/VAST) model parameters over time.
+ All other examples assume independence of 
***

1.  Use baseline knot number = 500. 
2.  Loop over RhoConfig specifications.
    + Intercept (Beta) RhoConfig[1:2] = (0) Independent among years, (2) Random walk, (4) Autoregressive
    + Spatio-temporal RE (Epsilon) RhoConfig[3:4] = (0) Independent among years, (2) Random walk, (4) Autoregressive
3.  Call [VAST](https://github.com/James-Thorson/VAST) in parallel across species, for the sake of efficency.
    + Only Save.RData. Then use cleanup function call to reduce storage requirements.
4.  Loop through and delete unnecessary folder structure, while building overall list object.
5.  Save output list object, to be loaded for summary and figure generation.


*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)



