# Tutorial for the R package *BBMV*

The purpose of the **BBMV** package is to fit highly flexible models for continuous traits evolving on phylogenies. The package implements two main models: the *FPK* model and the *BBMV* model.

Under the *FPK* model, a continuous trait evolves according to random diffusion (Brownian motion) and is also subject to an 'evolutionary potential' that creates a force that pulls the trait towards specific regions of the trait interval. In theory, this force can be of any conceivable shape but for the present implementation we have chosen a parametric shape for the potential of the shape: *V(x)=ax<sup>4</sup>+bx<sup>2</sup>+cx*. 

The *BBMV* model is a special case of *FPK* in which the trait is subject to an 'evolutionary potential' but also evolves between two reflective bounds. This parametrization is rather flexible since it allows for no force, directional trends, attraction towards a trait value within the interval, attraction towards the two bounds, or attraction towards several distinct trait values within the interval. In this tutorial we will see how to estimate the model parameters using maximum-likelihood and MCMC.

Our implementation of the *FPK* and *BBMV* models does **not** require that the phylogenetic tree you provide be ultrametric: you can use it with fossil data for example.

The **BBMV** package can be installed from CRAN:
```r
install.packages('BBMV')
library(BBMV)
```

## Simulation function
For this tutorial we will simulate data and then infer parameters of the *BBMV* model on this simulated dataset. We need the R package **geiger** to simulate phylogenetic trees.
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

Now we will use the function *Sim_FPK* to simulate a continuous trait evolving on the tree we just simulated. To do so we need to provide a rate of evolution (*sigma*), bounds on the trait interval (which might not influence the whole process because they are located far away), a value for the trait at the root of the tree (*x0*), and the evolutionary potential that we just simulated (*V6*): 

```r
TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=1,bounds=bounds)
hist(TRAIT,breaks=20)
```

## Maximum-likelihood estimation

### FPK model

Maximum-likelihood (ML) estimation of the parameters of the *FPK* model is done in two steps: (i) creating the likelihood function and (ii) finding its maximum.

The likelihood function is created using the function *lnL_FPK*, which takes the phylogenetic tree and the vector of trait values at the tips of the tree as main arguments. In addition, we need to specify how finely we want to discretize the trait interval: our implementation of the *FPK* process indeed works by divinding the continuous trait intervals into a regular grid of points ranging from the lower to the upper bound. The finer the discretization the better the accuracy in the calculation of the likelihood, but the longer it takes. Here we will only take 50 points to discretize the interval so that the test is quick, but more should be used if computational ressources allow it. For this example we will use the *Nelder-Mead* optimization routine, which seems to perform better than others in the tests we have made. Finally, we need to specify the shape of the potential which we want to fit. The most complex form has three parameters (see above) but we can fit simpler shapes by fixing unnecessary parameters to 0.

We'll start with the most general shape of the potential, which can accommodate up to two peaks in the macroevolutionary landscape: 

```r
ll_FPK4=lnL_FPK(tree,TRAIT,Npts=50,a=NULL,b=NULL,c=NULL) # the full model
```
Then we can actually fit simpler versions of the FPK model. Here is the Ornstein-Uhlenbeck model, with *V(x)=b.x^2+c.x*:
```r
ll_FPK2=lnL_FPK(tree,TRAIT,Npts=50,a=0,b=NULL,c=NULL)
```

And here is Brownian motion, with *V(x)=0*:
```r
ll_FPK0=lnL_FPK(tree,TRAIT,Npts=50,a=0,b=0,c=0)
```

Once these likelihood functions are created, we need to use the *find.mle_FPK* function to estimate their maxima. This function takes one single argument: a likelihood function created by *lnL_FPK*. Once each model is fitted, we can plot the macroevolutionary landscape estimated using the function *get.landscape.BBMV* and compare it with the macroevolutionary landscape that we simulated:

```r
fit4=find.mle_FPK(model=ll_FPK4)
get.landscape.FPK(fit=fit4)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit2=find.mle_FPK(model=ll_FPK2)
get.landscape.FPK(fit=fit2) # this shape of the landscape cannot have 2 peaks
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit0=find.mle_FPK(model=ll_FPK0) 
get.landscape.FPK(fit=fit0) # this one is forced to be flat
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))
```

Now let's compare the goodness of fit of these three models using AIC:
```r
fit4$aic 
fit2$aic
fit0$aic
```
*fit4* should be the best model by far since no other form of the FPK model can accommodate two peaks.

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

It actually appears that the process reaches stationarity on most branches of the tree since the characteristic time is smaller than the median branch length:
```r
summary(tree$edge.length)
```


In this case we've reached stationarity since a while. Next, we will have a look at the probability distribution at the root of the tree. Be careful with the y-scale of the plot: this might actually be quite flat!
```r
plot(fit4$root,type='l') 
```

We can also estimate the uncertainty around maximum-likelihood parameter estimates using the function *Uncertainty_FPK*, which produces graphs of the likelihood of the model as a function of the value of each parameter. In this function, the parameter *effort_uncertainty* determines how many values of each parameter will be evaluated. The function returns confidence interval that contain the 95% highest probability density around parameter estimates while fixing other parameters to their maximum likelihood estimate. We first look at the MLEs of parameters estimated, which are contained in the *$par* argument of a fitted FPK model:
```r
fit4$par
```
And from that we can choose the *scope* of the uncertainty search for each parameter so that they include your MLEs:
```r
Uncertainty_FPK(fit=fit4,tree,trait=TRAIT,Npts=50,effort_uncertainty= 100,scope_a=c(-1,10),scope_b=c(-5,5),scope_c=c(-2,2))
```
One particularly interesting result from this will be to see if the confidence intervals for each parameter of the potential (*a*, *b*, and *c*) contain 0.

Finally, we can fit the OU model using the package **geiger** to see if its likelihood matches the one calculated using BBMV. They should be quite close but remember that the FPK model uses an approximation of the likelihood, hence the implementation in **geiger** is more accurate:
```r
OU=fitContinuous(phy=tree,dat=TRAIT,model="OU")
OU$opt$lnL ; fit2$lnL 
```

We can also look at the AIC of the BM model, fitted using **geiger**:
```r
BM=fitContinuous(phy=tree,dat=TRAIT,model="BM")
BM$opt$aic
```

Ancestral character estimations for internal nodes of the tree can be obtained as follows:
```r
ACE_nodes=ACE_FPK(fit4,specific.point=NULL)
plot(ACE_nodes[[90]],type='l')
```

The function *ACE_FPK* returns a list with one element per internal node in the tree. Each element is a table giving trait values on the trait grid in the first column and their associated density in the second column. Here we have plotted the probability density of the trait at node '90' but we can also ask for an ACE at any point in the tree. Below we will ask for an ACE in the middle of the first branch (as ordered in the 'phylo' object), this is down by passing a vector with (i) the parent node, (ii) the child node, and (iii) the time from the begining of the branch to the 'specific.point' argument:
```r
ACE_1st_branch=ACE_FPK(fit4,specific.point=c(fit4$tree$edge[1,1],fit4$tree$edge[1,2],fit4$tree$edge.length[1]/2))
plot(ACE_1st_branch,type='l')
```

### BBMV model

In the **BBMV** package we can also fit a special case of the FPK model in which there are actual (reflective) bounds on the trait interval: the *BBMV* model. Here we will use a simulated dataset on which there is a trend towards one of the bounds of the trait interval:
```r
par(mfrow=c(1,1))
Vb=3*x
Vb_norm=exp(-Vb)/sum(exp(-Vb)*step_size)
plot(Vb_norm)
TRAITb= Sim_FPK(tree,x0=0,V=Vb,sigma=2,bounds=bounds)
hist(TRAITb,breaks=20)
```

Then we create four different likelihood functions using a variant of the *lnL_FPK* function, called *lnL_BBMV*. The four scenarios correspond to either 3, 2, 1 or 0 polynomial terms in the evolutionary potential, *V*: 
```r
ll_BBMV4=lnL_BBMV(tree,TRAITb,Npts=50,bounds=bounds,a=NULL,b=NULL,c=NULL)
ll_BBMV2=lnL_BBMV(tree,TRAITb,Npts=50,bounds=bounds,a=0,b=NULL,c=NULL)
ll_BBMV1=lnL_BBMV(tree,TRAITb,Npts=50,bounds=bounds,a=0,b=0,c=NULL)
ll_BBMV0=lnL_BBMV(tree,TRAITb,Npts=50,bounds=bounds,a=0,b=0,c=0) # this is the BBM model
```

Next we fit these four models and plot the macroevolutionary landscapes:
```r
fit4b=find.mle_FPK(model=ll_BBMV4)
get.landscape.FPK(fit=fit4b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit2b=find.mle_FPK(model=ll_BBMV2)
get.landscape.FPK(fit=fit2b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit1b=find.mle_FPK(model=ll_BBMV1)
get.landscape.FPK(fit=fit1b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit0b=find.mle_FPK(model=ll_BBMV0)
get.landscape.FPK(fit=fit0b)
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

### Measurement error in trait data

The FPK and BBMV models can be fitted while accounting for measurement error in the value of the trait at the tips of the tree. This error is not given as a standard-deviation or variance, but rather a vector of multiple measurements of the trait value is given for each tip. The trait object passed as argument to either *lnL_FPK* or *lnL_BBMV* thus becomes a list instead of a vector.

If you only have an estimate of the standard-error and the mean trait value, then you can provide random draws of the trait for each species, as done below. Here we introduce random measurement error in the values of the TRAIT simulated earlier under the FPK model.

```r
TRAIT2=list()
for (i in 1:length(TRAIT)){
  TRAIT2[[i]]=rnorm(n=runif(n=1,min=2,max=10),mean=TRAIT[i],sd=runif(n=1,min=0,max=0.03))
}
names(TRAIT2)=names(TRAIT)
```

We can have a look at how trait measurments look for the first 4 species: different tips have different number of measures and different levels of error.

```r
TRAIT2[c(1:4)] 
```

Then we can use the same functions as before to fit the model to this dataset. Here we only demonstrate the inclusion of measurement error when fitting the FPK model, but the same works for the BBMV model.

```r
ll_FPK4_with_ME=lnL_FPK(tree,TRAIT2,Npts=50,a=NULL,b=NULL,c=NULL) # the full model
fit4_with_ME=find.mle_FPK(model=ll_FPK4_with_ME)
```

Now we can compare the macroevolutionary landscapes estimated using both methods: they should be rather similar unless measurement error is large relative to the spread of the trait.

```r
par(mfrow=c(1,2))
get.landscape.FPK(fit=fit4)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))
get.landscape.FPK(fit=fit4_with_ME)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))
```
Uncertainty in parameters can also be estimated when fitting the model with measurement error:

```r
fit4_with_ME$par
Uncertainty_FPK(fit=fit4_with_ME,tree,trait=TRAIT2,Npts=50,effort_uncertainty= 100,scope_a=c(0,100),scope_b=c(-15,5),scope_c=c(-5,5))
```

## Markov Chain Monte Carlo estimation
We can also estimate parameters of the full model using an MCMC chain with the Metropolis Hastings algorithm using the function *MH_MCMC_FPK*. Here we will do only a few generations so that computation time is not too long but for analysing real datasets you should monitor convergence of the MCMC chain (see below).

Parameters of the MCMC functions are the following:
- tree and TRAIT: the phylogenetic tree and data vector
- bounds: the bounds on the trait interval
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

Let's load and use the MCMC function on the simulated dataset with two peaks (FPK model). Here we fix the bounds of the trait interval to be the same as when we fitted the FPK model: they are far away from the observed trait interval and thus do not influence the process, we are thus fitting the FPK model. If we would like to fit the BBMV model we could place them closer to the observed trait, e.g. bounds=c(min(TRAIT),max(TRAIT)).

```r
MCMC=MH_MCMC_FPK(tree,trait=TRAIT,bounds=fit4$par_fixed$bounds,Nsteps=10000,record_every=100,plot_every=100,Npts=50,pars_init=c(0,-4,-4,0,1),prob_update=c(0.2,0.25,0.25,0.25,0.05),verbose=TRUE,plot=TRUE,save_to='~/MCMC_FPK_test.Rdata',save_every=100,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F)
```

Now we can measure the effective sample size of the chain using the package **coda**. This value should be above, say, 100 for a chain to have converged and in addition we should run several chains and check that they have converged to the same posterior distribution. We can also plot the posterior distribution of model parameters. We remove the 50 first samples as burnin just for the example, but we should probably be running the chain for much longer and discard way more samples:
```r
library(coda)
apply(MCMC[-c(1:50),2:9],2,effectiveSize)
par(mfrow=c(2,3))
hist(log(MCMC[-c(1:50),2]/2),breaks=20,main='log(sigsq/2)',ylab=NULL)
hist(MCMC[-c(1:50),6],breaks=20,main='root',ylab=NULL)
plot(1,1)
hist(MCMC[-c(1:50),3],breaks=20,main='a (x^4 term)',ylab=NULL)
hist(MCMC[-c(1:50),4],breaks=20,main='b (x^2 term)',ylab=NULL)
hist(MCMC[-c(1:50),5],breaks=20,main='c (x term)',ylab=NULL)
```

In order to visualize the posterior of our MCMC analysis we can plot credible intervals (CI) of macroevolutionary landscapes estimated in the MCMC run. The function *get.landscape.FPK.MCMC* plots the median value of the macroevolutionary landscape across the posterior in a solid lines and draws a polygon that stretches between two quantiles, defined by probs.CI. Here we set a burnin fraction of 50% since the MCMC run was extremely short and we look at the 95% CI. 

```r
get.landscape.FPK.MCMC(chain=MCMC,bounds=fit4$par_fixed$bounds,Npts=100,burnin=0.5,probs.CI=c(0.025,0.975),COLOR_MEDIAN='red',COLOR_FILL='red',transparency=0.3,main='Macroevolutionary landscapes MCMC',ylab='N.exp(-V)',xlab='Trait',xlim=NULL,ylim=NULL)
```

We can add the macroevolutionary landscape fitted using maximum-likelihood to this plot:
```r
add.ML.landscape.FPK(fit=fit4,Npts=100,COLOR=1,LTY='dashed')
```

We can plot the distribution of the prior vs. posterior for each parameter of the model to see how much information came from the data. Even though the chain is far from having converged we might still get strong information from the likelihood (i.e., the data)... or we might still be strongly influenced by the initial point of the MCMC chain.
```r
par(mfrow=c(2,2))
posterior_vs_prior(chain=MCMC,param='a',Npts=100,burnin=0.2,type_prior='Normal',shape_prior=c(0,10))
posterior_vs_prior(chain=MCMC,param='b',Npts=100,burnin=0.2,type_prior='Normal',shape_prior=c(0,10))
posterior_vs_prior(chain=MCMC,param='c',Npts=100,burnin=0.2,type_prior='Normal',shape_prior=c(0,10))
posterior_vs_prior(chain=MCMC,param='sigsq',Npts=100,burnin=0.2,type_prior='Normal',shape_prior=c(0,10))
```

Finally, we can also estimate a simpler version of the model using MCMC. Here we will run a chain with the potential forced to be quadratic (i.e. an OU model). We do this by fixing the intial values of *a* to 0 and setting its probability of update to zero:
```r
MCMC_OU=MH_MCMC_FPK(tree,trait=TRAIT,bounds=fit4$par_fixed$bounds,Nsteps=10000,record_every=100,plot_every=100,Npts=50,pars_init=c(0,0,-4,0,1),prob_update=c(0.25,0,0.35,0.35,0.05),verbose=TRUE,plot=TRUE,save_to='~/MCMC_FPK_test.Rdata',save_every=100,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F)

apply(MCMC_OU[-c(1:50),2:9],2,effectiveSize)
par(mfrow=c(2,3))
hist(log(MCMC_OU[-c(1:50),2]/2),breaks=20,main='log(sigsq/2)',ylab=NULL)
hist(MCMC_OU[-c(1:50),6],breaks=20,main='root',ylab=NULL)
plot(1,1)
hist(MCMC_OU[-c(1:50),3],breaks=20,main='a (x^4 term)',ylab=NULL)
hist(MCMC_OU[-c(1:50),4],breaks=20,main='b (x^2 term)',ylab=NULL)
hist(MCMC_OU[-c(1:50),5],breaks=20,main='c (x term)',ylab=NULL)
```

MCMC estimation of the FPK model can also be done while incorporating measurement error in tip data. This is provided exactly as for ML estimation.

## Fit the *FPK* model on multiple clades at once
It is also possible to fit the *FPK* model on multiple clades together in order to statistically test whether they share a similar macroevolutionary landscape. This is an alternative to existing methods that infer heterogeneity in macroevolutionary dynamics across large phylogenies (e.g. packages **l1ou** or **bayou** for OU models). The advantage of this procedure is that it does not require a backbone tree connecting the different clades one wants to compare. Drawbacks are that clades are considered independent and that only *these* clades can be compared (*i.e.*, the algorithm does not explore the possibility that any subsclade has different dynamics). This method could also help in cases where one has data at hand for multiple small clades in which traits are believed to have evolved under similar dynamics: pooling clades together might improve estimation of the macroevolutionary landscape.

We first create a potential that we will use to simulate trait evolution: it has two peaks of very unequal heights.

```r
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x))
a=8 ; b=-4 ; c=1
V6=a*x^4+b*(x^2)+c*x
step_size=(max(bounds)-min(bounds))/(100-1)
V6_norm=exp(-V6)/sum(exp(-V6)*step_size)
par(mfrow=c(1,1))
plot(V6_norm,type='l')
```

Now we simulate a tree and a continuous trait for 3 independent clades. The trait evolves in the same macroevolutionary landscape for the 3 clades, but with different evolutionary rates (parameter *sigma* in the *Sim_FPK* function).

```r
tree=sim.bdtree(stop='taxa',n=25)
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=1,bounds=bounds)
tree1=tree ; TRAIT1=TRAIT

tree=sim.bdtree(stop='taxa',n=25)
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=0.5,bounds=bounds) 
tree2=tree ; TRAIT2=TRAIT

tree=sim.bdtree(stop='taxa',n=25)
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=0.1,bounds=bounds) 
tree3=tree ; TRAIT3=TRAIT
rm(tree) ; rm(TRAIT)
```

The phylogenies and trait vectors for each clade are put in a list. This is the format required for the multiclade functions.

```r
TREES=list(tree1,tree2,tree3)
TRAITS=list(TRAIT1,TRAIT2,TRAIT3)
```

We can have a look at the trait distribution in each clade: from left to right the evolutionary rate decreases and thus the macroevolutionary landscape is less well explored

```r
par(mfrow=c(1,3))
hist(TRAITS[[1]],breaks=20)
hist(TRAITS[[2]],breaks=20)
hist(TRAITS[[3]],breaks=20)
```

 Now we will fit three different scenarios of the FPK model with the most complex form of the potential, V(x)=a.x^4+b.x^2+c.x.

 1) All clades share the same macroevolutionary landscape and evolutionary rate :

```r 
testFPK4=lnl_FPK_multiclades_same_V_same_sig2(trees=TREES,traits=TRAITS,a=NULL,b=NULL,c=NULL,Npts=50)
fitFPK4=find.mle_FPK_multiple_clades_same_V_same_sig2(model=testFPK4,method='Nelder-Mead',init.optim=NULL)
```

2) All clades share the same macroevolutionary landscape but they have different evolutionary rates:

```r
testbFPK4=lnl_FPK_multiclades_same_V_different_sig2(trees=TREES,traits=TRAITS,a=NULL,b=NULL,c=NULL,Npts=50)
fitbFPK4=find.mle_FPK_multiple_clades_same_V_different_sig2(model=testbFPK4,method='Nelder-Mead',init.optim=NULL)
```

3) All clades have their own dynamics, both the macroevolutionary landscape and the evolutionary rate vary:
```r
fitmFPK4=fit_FPK_multiple_clades_different_V_different_sig2(trees=TREES,traits=TRAITS,a=NULL,b=NULL,c=NULL,Npts=50)
```

We can compare the fits of these three scenarios using AIC:
```r
fitFPK4$aic
fitbFPK4$aic
fitmFPK4$aic
```

The second model should have the lowest AIC (this is the one we simulated under). We can have a look at parameter estimates:

```r
fitbFPK4$par
```

In order to do an ACE, calculate the characteristic time of the FPK process or even plot the macroevolutionary landscape, we first need to separate results for each clade. The function 'reformat_multiclade_results' separates the results for each clade and returns a list with results for each clade, as done by 'find.mle_FPK'. 

```r
fits=reformat_multiclade_results(fitbFPK4)
```

And now we can do various things on clade #1.

```r
ace_tree1=ACE_FPK(fits$fit_clade_1)
par(mfrow=c(1,1))
plot(ace_tree1[[length(fits$fit_clade_1$tree$tip.label)+1]],type='l')
charac_time(fit=fits$fit_clade_1)
get.landscape.FPK(fit=fits$fit_clade_1)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))
```



