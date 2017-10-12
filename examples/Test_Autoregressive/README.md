# Test_Autoregressive
#### R Script: Test-Autoregressive.r

***
Compares influence of different autoregressive structures for [VAST](https://github.com/James-Thorson/VAST) model parameters over time.
***

1.  Use baseline knot numbers = 100, 500. 
2.  Loop over RhoConfig specifications.
    + Intercept (Beta) RhoConfig[1:2] = (0) Independent among years, (2) Random walk, (4) Autoregressive
    + Spatio-temporal RE (Epsilon) RhoConfig[3:4] = (0) Independent among years, (2) Random walk, (4) Autoregressive
3.  Call [VAST](https://github.com/James-Thorson/VAST) in parallel across species, for the sake of efficency.
    + Only Save.RData. Then use cleanup function call to reduce storage requirements.
4.  Loop through and delete unnecessary folder structure, while building overall list object.
5.  Saves index values and uncertainty, optimization results for retreiving AIC, and Report object.

*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)

* Note: We have remove EBS_SHELF Arrowtooth from this example, given error. 

***
#### Issues
1.  Error prevented exection for EBS_SHELF Arrowtooth Flounder.
    + Rho specs: int_RW-FE stRE_IaY Rho specs.
    + 500, but NOT 100 knots. 

***
#### Output Objects

Object                 | Description
-----------------------|-----------------------------------------------
vast_est.output.RData  | $vast_est - Index Values and Uncertainty
                       | $Opt - VAST optimizaton results
                       | $Report - VAST report object
vast_specs.csv         | Specifications for the different trial models, describing knot numbers and rho specifications on intercepts and spatio-temporal RE


