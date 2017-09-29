**To do, really**

Further test and then upload functions for joint inference of model on multiple independent trees.

Develop uncertainty function for the parameters of the multiclade fit

Develop MCMC algorithm for the multiclade fit

Develop/check functions for incorporating standard error in trait data --> do it directly by modifying the 'VectorPos_bounds' function?


**Less pressure...**

Generate error message when stat. distrib. does not converge to 0 at the bounds in FPK --> check 99% HPD and suggest to use lnL_BBMV instead with bounds further apart.

Generate helpful message if optimization failed, suggesting to either reduce Npts or change the optimization routine? Explained in the tutorial already.

rj-MCMC algorithm 'à la Bayou' for multiple FPK processes in the tree.

Function that enables users to DRAW any kind of potential?

Logo?
