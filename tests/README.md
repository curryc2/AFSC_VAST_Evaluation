# AFSC_VAST_Evaluation: Testing Scripts

This directory contains specific tests and explorations assicated with updates to TMB and VAST software. 

***

Name                            | Description | Result
--------------------------------|-------------|------------------------
Test_bias.correct_Efficency     | Testing two alternative methods for reducing RAM load for estimation with epsilon-estimator for bias correction and high (>300) knots. | Use of bias correction control options `bias.correct.control=list(nsplit=200, split=NULL, sd=FALSE)` allows estimation of 1,000 knots, but is time consuming (~7 hrs). Specifying model with bias correction only for the index (normal space) is effective and significantly faster (~1.5 hrs).

***

*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)