# Tutorial for the R package *BBMV*

The purpose of the **BBMV** package is to fit highly flexible models for continuous traits evolving on phylogenies. The package implements two main models: the *FPK* model and the *BBM+V* model.

Under the *FPK* model, a continuous trait evolves according to random diffusion (Brownian motion) and is also subject to an 'evolutionary potential' that creates a force that pulls the trait towards specific regions of the trait interval. In theory, this force can be of any conceivable shape but for the present implementation we have chosen a parametric shape for the potential of the shape: *V(x)=ax<sup>4</sup>+bx<sup>2</sup>+cx*. 

The *BBM+V* model is a special case of *FPK* in which the trait is subject to an 'evolutionary potential' but also evolves between two reflective bounds.

This parametrization is rather flexible since it allows for no force, directional trends, attraction towards a trait value within the interval, attraction towards the two bounds, or attraction towards several distinct trait values within the interval. In this tutorial we will see how to estimate the model parameters using maximum-likelihood and MCMC integration.

We first need to load the only R package on which **BBMV** depends: **ape**
```r
library(ape)
```
Then we need to source all functions in the **BBMV** package, which can be done via:
```r
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/get.landscape.BBMV_no_bounds.r',chdir=F)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/ML_functions.r',chdir=F)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/utils_BBMV.r',chdir=F)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/Simulate\ BBM+V.r',chdir = F)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/charac_time.r',chdir=F)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/Uncertainty_BBMV.r',chdir=F)
```

## Simulation function
For this tutorial we will simulate data and then infer parameters of the BBM+V model on this simulated dataset. We need the R package **geiger** to simulate phylogenetic trees.
```r
library(geiger)
tree=sim.bdtree(stop='taxa',n=50)
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
```
Here we have simulated a tree with only 50 tips. This is rather small but will allow for functions to run quickly. We have rescaled the total tree depth to 100 arbitrary time units for an easier interpretation of model parameters. 

Next we will simulate a macroevolutionary landscape with two peaks of equal heights, defined over the interval [-1.5,1.5]. We first create an evolutionary potential with two wells:
```r
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x))
V6=10*(x^4-0.5*(x^2)+0.*x)
```

From this evolutionary potential (*V6*), we calculate the normalized macroevolutionary landscape, which is the stationary distribution of the FPK model. It has two peaks of equal height:

```r
step_size=(max(bounds)-min(bounds))/(100-1)
V6_norm=exp(-V6)/sum(exp(-V6)*step_size)
par(mfrow=c(1,1))
plot(V6_norm)
```

Now we will use the function *Sim_BBMV* to simulate a continuous trait evolving on the tree we just simulated. To do so we need to provide a rate of evolution (*sigma*), bounds on the trait interval (which might not influence the whole process because they are located far away), a value for the trait at the root of the tree (*x0*), and the evolutionary potential that we just simulated (*V6*): 

```r
TRAIT= Sim_BBMV(tree,x0=0.5,V=V6,sigma=1,bounds=bounds)
hist(TRAIT,breaks=20)
```

## Maximum-likelihood estimation
Maximum-likelihood (ML) estimation of the parameters of the *FPK* model is done in two steps: (i) creating the likelihood function and (ii) finding its maximum.

The likelihood function is created using the function *lnL_FPK*, which takes the phylogenetic tree and the vector of trait values at the tips of the tree as main arguments. In addition, we need to specify how finely we want to discretize the trait interval: our implementation of the *FPK* process indeed works by divinding the continuous trait intervals into a regular grid of points ranging from the lower to the upper bound. The finer the discretization the better the accuracy in the calculation of the likelihood, but the longer it takes. Here we will only take 25 points to discretize the interval so that the test is quick, but more (at least 50) should be used when analyzing data seriously. For this example we will use the *Nelder-Mead* optimization routine, which seems to perform better than others in the tests we have made. Finally, we need to specify the shape of the potential which we want to fit. The most complex form has three parameters (see above) but we can fit simpler shapes by fixing unnecessary parameters to 0.

We'll start with the most general shape of the potential, which can accomodate up to two peaks in the macroevolutionary landscape: 

```r
ll_FPK4=lnL_FPK(tree,TRAIT,Npts=25,a=NULL,b=NULL,c=NULL) # the full model
```
Then we can actually fit simpler versions of the FPK model. Here is the Ornstein-Uhlenbeck model, with *V(x)=b.x^2+c.x*:
```r
ll_FPK2=lnL_FPK(tree,TRAIT,Npts=25,a=0,b=NULL,c=NULL)
```

And here is Brownian motion, with *V(x)=0*:
```r
ll_FPK0=lnL_FPK(tree,TRAIT,Npts=25,a=0,b=0,c=0)
```

Once these likelihood functions are created, we need to use the *find.mle_FPK* function to estimate their maxima. This function takes one single argument: a likelihood function created by *lnL_FPK*. Once each model is fit, we can plot the macroevolutionary landscape estimated using the function *get.landscape.BBMV* and compare it with the macroevolutionary landscape that we simulated:

```r
fit4=find.mle_FPK(model=ll_FPK4)
get.landscape.BBMV(fit=fit4)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit2=find.mle_FPK(model=ll_FPK2)
get.landscape.BBMV(fit=fit2) # this shape of the landscape cannot have 2 peaks
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit0=find.mle_FPK(model=ll_FPK0) 
get.landscape.BBMV(fit=fit0) # this one is forced to be flat
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))
```

Now let's compare the goodness of fit of these three models using AIC:
```r
fit4$aic 
fit2$aic
fit0$aic
```
*fit4* should be the best model by far since no other form of the FPK model can accomodate two peaks.

We can measure the time it takes to reach stationarity in the FPK model using the function *charac_time*:
```r
charac_time(fit=fit4)
charac_time(fit=fit2)
charac_time(fit=fit0)
```

And we should compare it with tree depth to see how far are we from stationarity: 
```r
max(branching.times(tree))
```

In this case we've reached stationarity since a while. Next, we will have a look at the probability distribution at the root of the tree. Be careful with the y-scale of the plot: this might actually be quite flat!
```r
plot(fit4$root,type='l') # 
```

We can also estimate the uncertainty around maximum-likelihood parameter estimates using the function *Uncertainty_BBMV*, which produces graphs of the likelihood of the model as a function of the value of each parameter. In this function, the parameter *effort_uncertainty* determines how many values of each parameter will be evaluated. The function returns confidence interval that contain the 95% highest probability density around parameter estimates while fixing other parameters to their maximum likelihood estimate. We first look at the MLEs of parameters estimated, which are contained in the *$par* argument of a fitted FPK model:
```r
fit4$par
```
And from that we can choose the *scope* of the uncertainty search for each parameter so that they include your MLEs:
```r
Uncertainty_BBMV(fit=fit4,tree,trait=TRAIT,Npts=25,effort_uncertainty= 100,scope_a=c(-1,10),scope_b=c(-5,5),scope_c=c(-2,2))
```
One particularly interesting result from this will be to see if the confidence intervals for each parameter of the potential (*a*, *b*, and *c*) contain 0.

Finally, we can fit the OU and BM models using the package **geiger** to see if likelihoods match with those calculated using BBMV. They should be quite close but remember that the FPK model uses an approximation of the likelihood, hence the implementation in **geiger** is more accurate:
```r
OU=fitContinuous(phy=tree,dat=TRAIT,model="OU")
OU$opt$lnL ; fit2$lnL # 
BM=fitContinuous(phy=tree,dat=TRAIT,model="BM")
BM$opt$lnL; fit0$lnL 
```

In the **BBMV** package we can also fit a special case of the FPK model in which there are actual (reflective) bounds on the trait interval: the *BBM+V* model. Here we will use a simulated dataset on which there is a trend towards one of the bounds of the trait interval:
```r
par(mfrow=c(1,1))
Vb=3*x
Vb_norm=exp(-Vb)/sum(exp(-Vb)*step_size)
plot(Vb_norm)
TRAITb= Sim_BBMV(tree,x0=0,V=Vb,sigma=2,bounds=bounds)
hist(TRAITb,breaks=20)
```

Then we create four different likelihood functions using a variant of the *lnL_FPK* function, called *lnL_BBMV*. The four scenarios correspond to either 3, 2, 1 or 0 polynomial terms in the evolutionary potential, *V*: 
```r
ll_BBMV4=lnL_BBMV(tree,TRAITb,Npts=25,bounds=bounds,a=NULL,b=NULL,c=NULL)
ll_BBMV2=lnL_BBMV(tree,TRAITb,Npts=25,bounds=bounds,a=0,b=NULL,c=NULL)
ll_BBMV1=lnL_BBMV(tree,TRAITb,Npts=25,bounds=bounds,a=0,b=0,c=NULL)
ll_BBMV0=lnL_BBMV(tree,TRAITb,Npts=25,bounds=bounds,a=0,b=0,c=0) # this is the BBM model
```

Next we fit these four models and plot the macroevolutionary landscapes:
```r
fit4b=find.mle_FPK(model=ll_BBMV4)
get.landscape.BBMV(fit=fit4b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit2b=find.mle_FPK(model=ll_BBMV2)
get.landscape.BBMV(fit=fit2b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit1b=find.mle_FPK(model=ll_BBMV1)
get.landscape.BBMV(fit=fit1b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit0b=find.mle_FPK(model=ll_BBMV0)
get.landscape.BBMV(fit=fit0b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))
```
We can also compare the AIC of these four models:
```r
fit4b$aic
fit2b$aic
fit1b$aic 
fit0b$aic
```
*fit1b* this should have the lowest AIC since it is the model we simulated... at least if we would have had a large enough tree. More complex models (fit4b and fit2b) will also do good job since they can accommodate this trend, but they have more parameters.




## Markov Chain Monte Carlo estimation
We can also estimate parameters of the full model using an MCMC chain with the Metropolis Hastings algorithm and a simple Gibbs sampler using the function *MH_MCMC_FPK*. Here we will do only a few generations so that computation time is not too long but for analysing real datasets you should monitor convergence of the MCMC chain (see below).

Parameters of the MCMC functions are the following:
- tree and TRAIT: the phylogenetic tree and data vector
- Nsteps: the number of generations in the MCMC chain
- record_every: the interval used for sampling the MCMC chain
- plot_every: the interval at which the chain is plotted (if plot=TRUE).
- Npts: the number of point on the grid
- pars_init: the initial parameters for starting the algorithm, c(log(sig2/2),a,b,c,x0). Be careful since x0 is actually the point on the grid (between 1 and Npts), not the actual root value
- prob_update: the relative frequencies of update of the different parameters of the model
- verbose: if TRUE, will print some generations of the chain to the screen
- plot: if TRUE, the chain is plotted from time to time
- save_to: file to which the chain is saved (can be useful in case the chain crashes)
- save_every: sets how often the chain is saved
- type_priors: the type of priors used, can be either normal (preferred) or uniform for log(sig2/2), a, b and c, ; and can only be discrete uniform for x0
- shape_priors: list that gives the shape for each prior. (mean,sd) for normal priors and (min,max) for continuous uniform priors. The shape is not specified for the root prior, since it is fixed to be discrete uniform on the grid.
- proposal_type: the type of proposal function, only uniform is available
- proposal_sensitivity: the width of the uniform proposal. The entire value for x0 gives how many steps at a time can be travelled on the trait grid (better to keep it to 1)

Let's load and use the MCMC function on the simulated dataset with two peaks:

```r
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/MCMC_function_BBMV.r',chdir=F)
MH_MCMC_FPK(tree,trait=TRAIT,bounds=c(-1.5,1.5),Nsteps=200000,record_every=100,plot_every=100,Npts=20,pars_init=c(0,-4,-4,0,1),prob_update=c(0.2,0.25,0.25,0.25,0.05),verbose=TRUE,plot=TRUE,save_to='~/Desktop/MCMC_FPK_test.Rdata',save_every=100,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F)
```

#########
EDIT FROM HERE

Now we can measure the effective sample size of the chain using the package *coda*. This value should be above, say, 100 for a chain to have converged and in addition we should run several chains and check that they have converged to the same posterior distribution. We can also plot the posterior distribution of the model parameters. We remove the 50 first samples as burnin just for fun, but we should probably be running the chain for much longer and discard way more samples:
```r
library(coda)
apply(MCMC[-c(1:50),2:11],2,effectiveSize)

par(mfrow=c(2,4))
hist(log(MCMC[-c(1:50),2]/2),breaks=100,main='log(sigsq/2)',ylab=NULL)
hist(MCMC[-c(1:50),3],breaks=100,main='a (x^4 term)',ylab=NULL)
hist(MCMC[-c(1:50),4],breaks=100,main='b (x^2 term)',ylab=NULL)
hist(MCMC[-c(1:50),5],breaks=100,main='c (x term)',ylab=NULL)
hist(MCMC[-c(1:50),6],breaks=100,main='root',ylab=NULL)
hist(MCMC[-c(1:50),9],breaks=100,main='lnprior',ylab=NULL)
hist(MCMC[-c(1:50),10],breaks=100,main='lnlik',ylab=NULL)
hist(MCMC[-c(1:50),11],breaks=100,main='quasi-lnpost',ylab=NULL)
```
Finally, we can also estimate a simpler version of the model using MCMC. Here we will run a chain with the potential forced to be linear (i.e. what we simulated). We do this by fixing the intial values of a and b to 0 and setting their probabilities of update to zero:
```r
MCMC_trend= MH_MCMC_V_ax4bx2cx_root_bounds(tree,trait=TRAIT,Nsteps=20000,record_every=100,plot_every=500,Npts_int=20,pars_init=c(-8,0,0,0,5,min(TRAIT),max(TRAIT)),prob_update=c(0.05,0.,0.,0.15,0.15,0.05,0.05),verbose=TRUE,plot=TRUE,save_to='testMCMC_linear.Rdata',save_every=1000,type_priors=c(rep('Normal',4),rep('Uniform',3)),shape_priors=list(c(0,2),c(0,2),c(0,2),c(0,2),NA,30,30),proposal_type='Uniform',proposal_sensitivity=c(1,0.5,0.5,0.5,1,1,1),prior.only=F)

# sample size and plots
apply(MCMC_trend[-c(1:50),c(2,5:11)],2,effectiveSize)
par(mfrow=c(2,4))
hist(log(MCMC_trend[-c(1:50),2]/2),breaks=100,main='log(sigsq/2)',ylab=NULL)
hist(MCMC[-c(1:50),5],breaks=100,main='c (x term)',ylab=NULL)
hist(MCMC_trend[-c(1:50),6],breaks=100,main='root',ylab=NULL)
hist(MCMC_trend[-c(1:50),7],breaks=100,main='bmin',ylab=NULL)
hist(MCMC_trend[-c(1:50),8],breaks=100,main='bmax',ylab=NULL)
hist(MCMC[-c(1:50),9],breaks=100,main='lnprior',ylab=NULL)
hist(MCMC_trend[-c(1:50),10],breaks=100,main='lnlik',ylab=NULL)
hist(MCMC_trend[-c(1:50),11],breaks=100,main='quasi-lnpost',ylab=NULL)
```


### Troubleshooting ML estimation
Numerical errors can occur when trying to fit the *BBM+V* model to empirical data. Here are a few suggestions that can help solving some issues:
- changing the optimization method used in the *fit_BBMV* function (e.g. from *Nelder-Mead* to *L-BFGS-B*)
- changing *Npts* to an odd number, or reducing it
- in cases where you see an error which looks like *error in solve.default ... reciprocal condition number =*, you can try decreasing the tolerance used when calling the *solve* function. This can be done by manually editing the *prep_mat_exp* function in the *BBM+V_functions_MLoptim.R* script
- if you encounter errors when using *fit_BBMV* with fixed bounds, then you can try changing the intial parameters used to start the optimization. This can be done as follows:
```r
BBM_x_init=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='linear',bounds=c(-5,5),init.optim = c(log(BBM_x$par$sigsq/2),0))
```

Also, when you are comparing models with different shapes of the evolutionary potential, remember to check that complex models should always have a higher likelihood than simpler models (which are nested). Optimization of the model is difficult, and this kind of situations unfortunately happen...
