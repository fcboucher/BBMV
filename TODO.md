**To do, really**

Remove unnecessary paramaters in MCMC help functions: 'trait,bounds' not used in 'log_prior_5pars_root_bounds', 'trait' not used in 'proposal_5pars_root_bounds'

Update Tutorial to show usage of functions with measurement error

Further test functions for joint inference of model on multiple independent trees with measurement error incorporated in it, then push to Github.

Develop uncertainty function for the parameters of the multiclade fit

Develop MCMC algorithm for the multiclade fit


**Less pressure...**

Generate error message when stat. distrib. does not converge to 0 at the bounds in FPK --> check 99% HPD and suggest to use lnL_BBMV instead with bounds further apart.

Generate helpful message if optimization failed, suggesting to either reduce Npts or change the optimization routine? Explained in the tutorial already.

rj-MCMC algorithm 'à la Bayou' for multiple FPK processes in the tree.

Function that enables users to DRAW any kind of potential?

Logo?
