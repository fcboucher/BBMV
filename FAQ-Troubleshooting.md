# FAQ/Troubleshooting

### How do I fit the bounded Brownian motion model?
The bounded Brownian motion model (Boucher & Démery 2016 Syst. Biol.) is a special case of the *FPK* model in which there are bounds on the trait values and the macroevolutionary landscape is flat. Its likelihood can be calculated using the *BBMV* package as follows:
```r
ll_BBMV0=lnL_BBMV(tree,TRAITb,Npts=50,bounds=bounds,a=0,b=0,c=0)
```
However, this implies that you know the bounds that act on trait values. Estimating them is possible using functions of the [*BBM* repository](https://github.com/fcboucher/BBM), but the optimization procedure used is less reliable than that of the *BBMV* package. Although we lack an analytical proof of it, it seems that the minimum and maximum of observed trait values at the tips of the phylogeny are the ML estimators of the bounds of the trait interval (see Boucher & Démery 2016 Syst. Biol. for details). If you don't have an a priori of which bounds to expect, you could thus use the *BBMV* package to fit the bounded Brownian motion model as follows:
```r
ll_BBMV0=lnL_BBMV(tree,TRAITb,Npts=50,bounds=c(min(TRAITb),max(TRAITb)),a=0,b=0,c=0)
```

### ML estimation
Numerical errors can occur when trying to fit the *FPK* model to empirical data. Here are a few suggestions that can help solving some issues:
- changing the optimization method used in the *find.mle_FPK* function (e.g. from *Nelder-Mead* to *L-BFGS-B*)
- changing *Npts* to an odd number, or reducing it
- changing the initial parameters used to start the optimization, especially reducing the value of the first parameter (the diffusion coefficient, log(sigsq/2)). This can be done as follows:
```r
fit4b=find.mle_FPK(model=ll_FPK4,init.optim = c(-12,0,1,-5))
```
- using three different sets of initial parameters used to start the optimization, which can be done by setting the *safe* argument to TRUE:
```r
fit4b=find.mle_FPK(model=ll_FPK4,safe=TRUE)
```

In cases where you see an error which looks like *error in solve.default ... reciprocal condition number =*, you can try decreasing the tolerance used when calling the *solve* function, which we use for inverting the diffusion matrix. This can be done by manually editing the *prep_mat_exp* function in the *utils_BBMV.r* script.

**Pay attention to the likelihoods of nested models**

When you are comparing models with different shapes of the evolutionary potential, remember to check that complex models should always have a higher likelihood than simpler models (which are nested within the complex ones). Optimization of the model is difficult, and this kind of situations unfortunately happens...

### MCMC estimation

In order to ensure that the Gibbs sampler mixes well you need to fine-tune the proposal function (the function that determines how you can move to one parameter value to another in the next step of the MCMC chain). This is done by specifying a vector of values as the *proposal_sensitivity* argument and can be quite tricky. One solution to get an idea of the values this parameter should take is to make an test MCMC run while sampling the prior only. This is rather quick since the likelihood needs not be calculated (which takes the most computing time) and can be done by setting *prior.only=F*. At least by doing that you would be able to rule out values of the proposal sensitvity that are too small for the prior to be fully explored in a reasonable number of steps of the MCMC algorithm.