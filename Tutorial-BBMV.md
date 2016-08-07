# Tutorial for the R package *BBMV*

The purpose of **BBMV** is to fit highly flexible models for continuous traits evolving on phylogenies. Under the BBM+V model, a continuous trait evolves between two reflective bounds according to random diffusion (Brownian motion). In addition, the trait is   subject to an 'evolutionary potential', which creates a force that pulls the trait towards specific regions of the trait interval. In theory, this force can be of any conceivable shape but for the present implementation we have chosen a parametric shape for the potential of the shape:
```r
V(x)=ax^4+bx^2+cx 
```

This parametrization is rather flexible since it allows for no force, directional trends, attraction towards a trait value within the interval, attraction towards the two bounds, or attraction towards several distinct trait values within the interval. In this tutorial we will see how to estimate the model parameters using maximum-likelihood and MCMC integration.

We first need to load the only R package on which **BBMV** depends: *ape*
```r
library(ape)
```
Then we need to source all functions in the **BBMV** package, which can be done via:
```r
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/BBM+V_functions_MLoptim.R', chdir = TRUE)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/charac_time.R', chdir = TRUE)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/MCMC functions BBM+V.R', chdir = TRUE)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/plot.landscape.BBMV.R', chdir = TRUE)
source('SET_PATH_TO_THIS_FILE_ON_YOUR_COMPUTER/Simulate BBM+V.R', chdir = TRUE)
```

### Simulation function
For this tutorial we will simulate data and then infer parameters of the BBM+V model on this simulated dataset. We need the R package *geiger* to simulate phylogenetic trees.
```r
library(geiger)
tree=sim.bdtree(stop='taxa',n=20)
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
```
Here we have simulated a tree with only 20 tips. This is rather small but will allow for functions to run quickly. We have rescaled the total tree depth to 100 arbitrary time units for an easier interpretation of model parameters. 

Next we will use the function *Sim_BBMV* to simulate a continuous trait evolving on the tree we just simulated. To do so we need to provide a rate of evolution (*sigma*), bounds on the trait interval, a value for the trait at the root of the tree (*x0*), and an evolutionary potential (*V*). Here we simulate a potential linearly increasing towards *higher* values of the trait. This will create a force towards *smaller* values of the trait, which we can see as the distribution of values at the tips of the tree is strongly left skewed.
```r
TRAIT= Sim_BBMV(tree,x0=0,V=seq(from=0,to=5,length.out=50),sigma=10,bounds=c(-5, 5))
hist(TRAIT,breaks=20)
```
### Maximum-likelihood estimation
The function to perform maximum-likelihood (ML) estimation of model parameters is *fit_BBMV*. It takes the phylogenetic tree and the vector of trait values at the tips of the tree as main arguments. In addition, we need to specify how finely we want to discretize the trait interval: the BBM+V process indeed works by divinding the continuous trait intervals into a regular grid of points ranging from the lower to the upper bound. The finer the discretization the better the accuracy in the calculation of the likelihood, but the longer it takes. Here we will only take 20 points to discretize the interval so that the test is quick, but more (at least 50) should be used when analyzing data seriously. For this example we will use the *Nelder-Mead* optimization routine, which seems to perform better than others in the tests we have made. Finally, we need to specify the shape of the potential which we want to fit. The most complex form has three parameters (see above) but we can fit simpler shapes.

We'll start with a flat potential, i.e. there is no force acting on the trait and the trait only evolves according to bounded Brownian Motion (BBM):

```r
BBM=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='flat')
BBM$par
```
The *$par* element of the model fit gives the ML values of the parameters. Here we have the evolutionary rate (sigsq), the value of the trait at the root and the positions of the bounds of the trait interval.
Now we can fit increasingly complex models starting with a linear potential, adding a quadratic term, and finally fitting the full model with a *x^4* term:
```r
BBM_x=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='linear')
BBM_x$par

BBM_x2x=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='quadratic')
BBM_x2x$par

BBM_full=fit_BBMV(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T,V_shape='full')
BBM_full$par
```

```r
# Fit other classic models of evolution implemented in package {geiger}
BM=fitContinuous(phy=tree,dat=TRAIT,model='BM') # Brownian motion with no bounds
OU=fitContinuous(phy=tree,dat=TRAIT,model='OU') # Ornstein-Uhlenbeck process with a single optimum

# AIC comparison of all the models fitted
BBM$aicc
BBM_x$aicc # best model... normally (this is the one we simulated)
BBM_x2x$aicc
BBM_full$aicc
BM$opt$aicc
OU$opt$aicc

# Now plot the adaptive landscape estimated by the best model
plot.landscape.BBMV(model=BBM_x,Npts=100)

# Plot landscapes estimated by all 4 versions of BBM+V fitted...
plot.multiple.landscapes.BBMV(models=list(BBM,BBM_x, BBM_x2x, BBM_full),Npts=100,ylim=c(0,0.06))

# measures times to reach stationarity
charac_time(Npts=20,BBM)
charac_time(Npts=20, BBM_x)
charac_time(Npts=20, BBM_x2x)
charac_time(Npts=20, BBM_full)
# compare them with tree depth: how far are we from equilibrium?
max(branching.times(tree))


###############################################
##### Markov Chain Monte Carlo estimation #####
###############################################

# Estimate parameters of the full model using an MCMC chain with the Metropolis Hastings algorithm and a simple Gibbs sampler
# You need to specify the file to which the chain is saved ('save_to' parameter)
# Here we will do only a few generations so that computation time is not too long but for analysing real datasets you should monitor MCMC convergence (see below)

# Parameters of the MCMC functions are the following:
# tree and TRAIT: the phylogenetic tree and data vector
# Nsteps: the number of generations in the MCMC chain
# record_every: the interval used for sampling the MCMC chain
# plot_every: the interval at which the chain is plotted (if plot=TRUE).
# Npts_int: the number of point on the grid between min(trait) and max(trait)
# pars_init: the initial parameters for starting the algorithm, c(log(sig2/2),a,b,c,x0,Bmin,Bmax). Be careful since x0 is actually the point on the grid (between 1 and Npts_int), not the actual root value
# prob_update: the relative frequencies of update of the different parameters of the model
# verbose: if TRUE, will print some generations of the chain to the screen
# plot: if TRUE, the chain is plotted from time to time
# save_to: file to which the chain is saved (can be useful in case the chain crashes)
# save_every: sets how often the chain is saved
# type_priors: the type of priors used, can be either normal (preferred) or uniform for log(sig2/2), a, b and c, ; and can only be discrete uniform for bounds and x0
# shape_priors: list that gives the shape for each prior. (mean,sd) for normal priors and (min,max) for continuous uniform priors. The shape is not specified for the root prior, since it is fixed to be discrete uniform on the grid. Values for the priors on the bounds (discrete uniform) give the maximum number of points that can be added on the trait grid outside of the observed trait interval
# proposal_type: the type of proposal function, only uniform is available
# proposal_sensitivity: the width of the uniform proposal. The entire value for x0, Bmin, and Bmax give how many steps at a time can be travelled on the trait grid (better to keep it to 1)

MCMC= MH_MCMC_V_ax4bx2cx_root_bounds(tree,trait=TRAIT,Nsteps=20000,record_every=100,plot_every=500,Npts_int=20,pars_init=c(-8,0,0,0,5,min(TRAIT),max(TRAIT)),prob_update=c(0.05,0.3,0.3,0.15,0.15,0.05,0.05),verbose=TRUE,plot=TRUE,save_to='~/Desktop/testMCMC1.Rdata',save_every=1000,type_priors=c(rep('Normal',4),rep('Uniform',3)),shape_priors=list(c(0,2),c(0,2),c(0,2),c(0,2),NA,30,30),proposal_type='Uniform',proposal_sensitivity=c(1,0.5,0.5,0.5,1,1,1),prior.only=F)


# Explore MCMC outputs
library(coda)
apply(MCMC[-c(1:50),2:11],2,effectiveSize) # Effective Sample Size for sampling of parameters, ideally we should aim for something >100. Here we have removed the 50 first samples as burnin.

# plot posterior distributions of parameters
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

# Now we will run a chain with the potential forced to be linear (i.e. what we simulated)
# We do this by fixing the intial values of a and b to 0 and setting their probabilities of update to zero: they will never be updated
MCMC_trend= MH_MCMC_V_ax4bx2cx_root_bounds(tree,trait=TRAIT,Nsteps=20000,record_every=100,plot_every=500,Npts_int=20,pars_init=c(-8,0,0,0,5,min(TRAIT),max(TRAIT)),prob_update=c(0.05,0.,0.,0.15,0.15,0.05,0.05),verbose=TRUE,plot=TRUE,save_to='~/Desktop/testMCMC1.Rdata',save_every=1000,type_priors=c(rep('Normal',4),rep('Uniform',3)),shape_priors=list(c(0,2),c(0,2),c(0,2),c(0,2),NA,30,30),proposal_type='Uniform',proposal_sensitivity=c(1,0.5,0.5,0.5,1,1,1),prior.only=F)

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
