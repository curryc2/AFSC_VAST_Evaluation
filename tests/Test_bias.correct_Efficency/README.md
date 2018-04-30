# Test-bias.correct-Efficency.r

Trial implementation of several alternatives for reducing RAM use for bias.correction with large numbers of knots.

##  Methods

Option Number | Description                     | Code                            | Result
--------------|---------------------------------|---------------------------------|-------------------------
1 | Original standard implementation. | `bias.correct.control` not used.  | Estimation fails on GOA species with `bias.correct=TRUE` and knots > ~ 300 on machine with 32 GB RAM.
2 | Use `bias.correct.control` to split bias correction estimation into many (`nsplit`) smaller processes. | `bias.correct.control=list(nsplit=200, split=NULL, sd=FALSE)` | Can estimate model with 1,000 knots and bias.correct=TRUE **without** running out of memory and failing. However estimation is **slow** ~7 hours.
3 | Estimate model without bias correction, **then** conduct bias correction for specified list of model parameters. | `SD = sdreport( obj=Obj, par.fixed=Opt$par, hessian.fixed=h, bias.correct=TRUE, bias.correct.control=list(sd=FALSE, split=Which, nsplit=NULL) )` | Estimation successful and was **significantly** faster at ~1.5 hours.
4 | Use updated specification in `bias.correct.control` that allows for specification of specific parameters on which to use epsilon-estimatior, through the `split` argument. | `bias.correct.control=list(sd=FALSE, split="Index_cyl", nsplit=100)`  | **Faster** than **Option 2**, taking ~ 3 hours.

* Current specifications are not considered final. 

***

*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)



