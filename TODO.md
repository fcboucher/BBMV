**To do, really**

In *lnL_BBMV* and *lnL_FPK* add check of match btw tip names and trait names

In *fin.mle_FPK* change the default starting point for optim: very small values of dCoeff work much better (init.optim=c(0,rep(0,model$ncoeff)))

Develop function for joint inference of model on multiple independent trees.

**Less pressure...**

Generate error message when stat. distrib. does not converge to 0 at the bounds in FPK --> check 99% HPD and suggest to use lnL_BBMV instead with bounds further apart.

Develop/check functions for incorporating standard error in trait data?

Generate helpful message if optimization failed, suggesting to either reduce Npts or change the optimization routine? Explained in the tutorial already.

rj-MCMC algorithm 'à la Bayou' for multiple FPK processes in the tree.

Function that enables users to DRAW any kind of potential?

Logo?
