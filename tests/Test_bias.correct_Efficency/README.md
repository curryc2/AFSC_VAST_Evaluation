# Test-bias.correct-Efficency.r

Trial implementation of several alternatives for reducing RAM use for bias.correction with large numbers of knots.

##  Methods

Description                     | Code                            | Result
--------------------------------|---------------------------------|-------------------------
Use `bias.correct.control` to split bias correction estimation into many (`nsplit`) smaller processes. | `bias.correct.control=list(nsplit=200, split=NULL, sd=FALSE)` | Can estimate model with 1,000 knots and bias.correct=TRUE **without** running out of memory and failing. However estimation is **slow** ~7 hours.
Estimate model without bias correction, **then** conduct bias correction for specified list of model parameters. | `SD = sdreport( obj=Obj, par.fixed=Opt$par, hessian.fixed=h, bias.correct=TRUE, bias.correct.control=list(sd=FALSE, split=Which, nsplit=NULL) )` | Estimation successful and was **significantly** faster at ~1.5 hours.

* Current specifications are not considered final. 

***

*	Note: Wrapper functions rely on calls to existing functions developed by J. Thorson [VAST](https://github.com/James-Thorson/VAST)



