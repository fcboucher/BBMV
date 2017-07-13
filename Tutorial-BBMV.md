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
The function to perform maximum-likelihood (ML) estimation of model parameters is *fit_BBMV*. It takes the phylogenetic tree and the vector of trait values at the tips of the tree as main arguments. In addition, we need to specify how finely we want to discretize the trait interval: our implementation of the *BBM+V* process indeed works by divinding the continuous trait intervals into a regular grid of points ranging from the lower to the upper bound. The finer the discretization the better the accuracy in the calculation of the likelihood, but the longer it takes. Here we will only take 20 points to discretize the interval so that the test is quick, but more (at least 50) should be used when analyzing data seriously. For this example we will use the *Nelder-Mead* optimization routine, which seems to perform better than others in the tests we have made. Finally, we need to specify the shape of the potential which we want to fit. The most complex form has three parameters (see above) but we can fit simpler shapes.

We'll start with a flat potential, i.e. there is no force acting on the trait and the trait only evolves according to bounded Brownian Motion (*BBM*):

```r
BBM=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='flat')
BBM$par
```
The *$par* element of the model fit gives the ML values of the parameters. Here we have the evolutionary rate (sigsq), the value of the trait at the root and the positions of the bounds of the trait interval.
Now we can fit increasingly complex models starting with a linear potential, adding a quadratic term, and finally fitting the full model with a *x<sup>4</sup>* term:
```r
BBM_x=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='linear')
BBM_x$par

BBM_x2x=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='quadratic')
BBM_x2x$par

BBM_full=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='full')
BBM_full$par
```
We can also compare our models with classic models of evolution like Brownian Motion (without bounds) and the Ornstein-Uhlenbeck process since likelihoods are comparable with the ones calculated in the **geiger** package:
```r
BM=fitContinuous(phy=tree,dat=TRAIT,model='BM') # Brownian motion with no bounds
OU=fitContinuous(phy=tree,dat=TRAIT,model='OU') # Ornstein-Uhlenbeck process with a single optimum
```
Yes, calculations in **geiger** are much faster than in **BBMV**. This is (mostly) because both BM and OU produce trait distribution that are multivariate normal, which simplifies calculations a lot. Unfortunately, this is not the case for the *BBM+V* model (trait distributions can anyway not be normal since the trait interval is bounded).

Now we will compare the fit of all of these models by looking at their AICs corrected for small sample sizes:
```r
BBM$aicc
BBM_x$aicc
BBM_x2x$aicc
BBM_full$aicc
BM$opt$aicc
OU$opt$aicc
```
*BBM_x* should be the model with the lowest AICc since this is the model we used for simulating the data. However, since we have a rather small dataset (20 tips) and since *BBM+V* is highly stochastic it might not always be the case. If you're not convinced, try running an example with 100 tips instead of 20.

If we want, we can also fix bounds that we think make sense: this can be sensible in some applications, for example if our trait is a probability or is naturally bounded. This can be done by specifying bounds when calling the *fit_BBMV* function: 
```r
BBM_x_fixed=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='linear',bounds=c(-5,5))
BBM_x_fixed$aicc
```

The **BBMV** package has a function for plotting what we call the 'macroevolutionary landscape' estimated by the model. The macroevolutionary landscape is simply the opposite of the evolutionary potential *landscape(x)=-V(x)*. Here, we again need to specify the number of points used to discretize the trait interval but this is just for plotting purposes:
```r
get.landscape.BBMV(model=BBM_x,Npts=100)
```
As for adaptive landscapes in quantitative genetics, we see peaks towards which trait values are attracted and valleys from which traits are repulsed. The plot shows you the macroevolutionary landscape over the whole trait interval, from the lower to the upper bound. We can also plot the macroevolutionary landscapes inferred by the four different versions of *BBMV* we have fitted (notice the flat landscape imposed in the first model):
```r
get.multiple.landscapes.BBMV(models=list(BBM,BBM_x, BBM_x2x, BBM_full),Npts=100,ylim=c(0,0.06))
```
An important measure in the *BBM+V* process is the time it takes for the process to reach stationarity. This is quite similar to the measure of the phylogenetic half-life for an OU process and we label it the *characteristic time*. Comparing this value to the total tree depth (100 in this example) gives us an idea of how far we are from stationarity:
```r
charac_time(Npts=20,BBM)
charac_time(Npts=20, BBM_x)
charac_time(Npts=20, BBM_x2x)
charac_time(Npts=20, BBM_full)
```

The *$ACE* element of the model fit gives the probability distribution of ancestral trait values. It is a list, with one table for each node, and elements of this list are numbered like nodes in the *phylo* object (the tree we used). Each of these ACE tables has two columns: the first ones gives the position of each point on the trait grid that was used for calculations and the second one the associated probability that the trait has this value. Here we can look at the probability distribution at the root of the tree:
```r
plot(BBM_x$ACE[[21]],type='l')
```
As shown in our [manuscript](https://github.com/fcboucher/BBMV/blob/master/Boucher_et_al_main_text.pdf), these ancestral character estimations will often be uninformative, especially when the process has reached stationarity. However, this might be useful in some cases.  

Finally, we can estimate the uncertainty around maximum-likelihood parameter estimates. This is done using the function *Uncertainty_BBMV*, which takes has input a model fitted using *fit_BBMV*, the phylogenetic tree, and the trait vector. The parameter 'effort_uncertainty' determines how many values of each parameter will be evaluated. The function produces graphs of the likelihood of the model as a function of the value of each parameter (the ML estimate is shown with a red line) and returns confidence intervals that contain the 95% highest probability density around parameter estimates while fixing other parameters to their maximum likelihood estimate. One particularly interesting result from this will be to see if the confidence intervals for each parameter of the potential (*a*, *b*, and *c*) contain 0.
```r
Uncertainty_BBMV(BBM,tree,trait= TRAIT,Npts=20,effort_uncertainty= 100)
Uncertainty_BBMV(BBM_x,tree,trait= TRAIT,Npts=20,effort_uncertainty= 100)
Uncertainty_BBMV(BBM_x2x,tree,trait=TRAIT,Npts=20,effort_uncertainty= 100)
Uncertainty_BBMV(BBM_full,tree,trait= TRAIT,Npts=20,effort_uncertainty= 100)
```

## Markov Chain Monte Carlo estimation
We can also estimate parameters of the full model using an MCMC chain with the Metropolis Hastings algorithm and a simple Gibbs sampler. This is done through the *MH_MCMC_V_ax4bx2cx_root_bounds* function. For explanations on each parameter the function takes as input please have a look at the [manual of the **BBMV** package](https://github.com/fcboucher/BBMV/blob/master/BBMV-manual.pdf).

Here we will run a quick example with only 20,000 generations and default parameters for the priors and proposal functions. In verbose mode, we get the state of the chain printed to the screen every at every sampled generation. If you allow plots, you will also see the trace of the chain:

```r
MCMC= MH_MCMC_V_ax4bx2cx_root_bounds(tree,trait=TRAIT,Nsteps=20000,record_every=100,plot_every=500,Npts_int=20,pars_init=c(-8,0,0,0,5,min(TRAIT),max(TRAIT)),prob_update=c(0.05,0.3,0.3,0.15,0.15,0.05,0.05),verbose=TRUE,plot=TRUE,save_to='testMCMC.Rdata',save_every=1000,type_priors=c(rep('Normal',4),rep('Uniform',3)),shape_priors=list(c(0,2),c(0,2),c(0,2),c(0,2),NA,30,30),proposal_type='Uniform',proposal_sensitivity=c(1,0.5,0.5,0.5,1,1,1),prior.only=F)
```
Now we can measure the effective sample size of the chain using the package *coda*. This value should be above, say, 100 for a chain to have converged and in addition we should run several chains and check that they have converged to the same posterior distribution. We can also plot the posterior distribution of the model parameters. We remove the 50 first samples as burnin just for fun, but we should probably be running the chain for much longer and discard way more samples:
```r
library(coda)
apply(MCMC[-c(1:50),2:11],2,effectiveSize)

par(mfrow=c(2,5))
hist(log(MCMC[-c(1:50),2]/2),breaks=100,main='log(sigsq/2)',ylab=NULL)
hist(MCMC[-c(1:50),3],breaks=100,main='a (x^4 term)',ylab=NULL)
hist(MCMC[-c(1:50),4],breaks=100,main='b (x^2 term)',ylab=NULL)
hist(MCMC[-c(1:50),5],breaks=100,main='c (x term)',ylab=NULL)
hist(MCMC[-c(1:50),6],breaks=100,main='root',ylab=NULL)
hist(MCMC[-c(1:50),7],breaks=100,main='bmin',ylab=NULL)
hist(MCMC[-c(1:50),8],breaks=100,main='bmax',ylab=NULL)
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
