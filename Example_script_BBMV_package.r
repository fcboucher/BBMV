rm(list=ls())
#Below is an example of the use of the BBMV package with simulated data:
# The BBMV package only depends on the package 'ape', which we need to load
library(ape)

# You first need to source all the functions of the package that we will need
source('/Users/florianboucher/Documents/Flo_BACKUPS/Travail/BBM\ plus\ potentiel/BBMV_Github/R/get.landscape.BBMV_no_bounds.r',chdir=F)
source('/Users/florianboucher/Documents/Flo_BACKUPS/Travail/BBM\ plus\ potentiel/BBMV_Github/R/ACE_FPK.r',chdir=F)
source('/Users/florianboucher/Documents/Flo_BACKUPS/Travail/BBM\ plus\ potentiel/BBMV_Github/R/ML_functions.r',chdir=F)
source('/Users/florianboucher/Documents/Flo_BACKUPS/Travail/BBM\ plus\ potentiel/BBMV_Github/R/utils_BBMV.r',chdir=F)
source('/Users/florianboucher/Documents/Flo_BACKUPS/Travail/BBM\ plus\ potentiel/BBMV_Github/R/Simulate\ BBM+V.r',chdir = F)
source('/Users/florianboucher/Documents/Flo_BACKUPS/Travail/BBM\ plus\ potentiel/BBMV_Github/R/charac_time.r',chdir=F)
source('/Users/florianboucher/Documents/Flo_BACKUPS/Travail/BBM\ plus\ potentiel/BBMV_Github/R/Uncertainty_BBMV.r',chdir=F)

# Simulate data: tree + continuous trait
library(geiger) # we will use geiger for simulating the tree
tree=sim.bdtree(stop='taxa',n=50) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree to a total depth of 100

# Here we define a macroevolutionary landscape with two peaks of equal heights, defined over the interval [-1.5,1.5]
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x)) # the bounds we use for simulating: they are just here for technical purposes but are not reached (see the plot of V6_norm below) and do not influence the process
V6=10*(x^4-0.5*(x^2)+0.*x) # this is the evolutionary potential: it has two wells

# We calculate the normalized macroevolutionary landscape and plot it
step_size=(max(bounds)-min(bounds))/(100-1)
V6_norm=exp(-V6)/sum(exp(-V6)*step_size) # the step size on the grid
par(mfrow=c(1,1))
plot(V6_norm) # two peaks of equal height

TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=0.5,bounds=bounds) # TRAIT simulated on the tree, evolving on the macroevolutionary landscape defined above: for that you need to source the function 'Sim_BBMV.R'
hist(TRAIT,breaks=20) # the distribution of the trait at the tips of the tree: it should reflect the landscape simulated... more or less

#####################################
#####################################
# Maximum-likelihood inference
#####################################
#####################################

#####################################
############ FPK model ##############
#####################################

# We will fit the FPK model to the data simulated using maximum-likelihood: this is done in two steps: (i) creating the likelihood function and (ii) finding its maximum
# The evolutionary potential can have three parameters when it is of the form V(x)=a.x^4+b.x^2+c.x, which allows for two peaks on the macroevolutionary landscape, but you can fix some parameters (to zero usually) to fit simpler versions of the model
# We use only 25 points for discretizing the trait interval to make it faster, but more points should be used on empirical datasets

# 1) We first need to prepare likelihood functions
ll_FPK4=lnL_FPK(tree,TRAIT,Npts=25,a=NULL,b=NULL,c=NULL) # the full model
ll_FPK2=lnL_FPK(tree,TRAIT,Npts=25,a=0,b=NULL,c=NULL) # we fix the x^4 coefficient to 0: this is actually the OU model
ll_FPK0=lnL_FPK(tree,TRAIT,Npts=25,a=0,b=0,c=0) # we fix all coefficients to 0:this is actually the BM model

# 2) fit these models and plot the landscapes
fit4=find.mle_FPK(model=ll_FPK4)
get.landscape.FPK(fit=fit4)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit2=find.mle_FPK(model=ll_FPK2)
get.landscape.FPK(fit=fit2) # this shape of the landscape cannot have 2 peaks
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit0=find.mle_FPK(model=ll_FPK0) 
get.landscape.FPK(fit=fit0) # this one is forced to be flat
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

# Compare model fits using AIC
fit4$aic # this should be the best model by far: no other form of the FPK model can accomodate two peaks
fit2$aic
fit0$aic

# measure the time it takes to reach stationarity in the FPK model
charac_time(fit=fit4)
charac_time(fit=fit2)
charac_time(fit=fit0)
# compare it with tree depth and typical branch length: how far are we from equilibrium? --> We've reached it!
max(branching.times(tree))
summary(tree$edge.length)

# We can have a look at the probability distribution at the root of the tree
plot(fit4$root,type='l') # be careful with the y-scale of the plot: this might actually be quite flat!

# We can also estimate the uncertainty around maximum-likelihood parameter estimates
# The function 'Uncertainty_BBMV' will produce graphs of the likelihood of the model as a function of the value of each parameter
# 'effort_uncertainty' determines how many values of each parameter will be evaluated
# The function returns confidence interval that contain the 95% highest probability density around parameter estimates while fixing other parameters to their maximum likelihood estimate
fit4$par # the MLEs of parameters
# You can specify the scope of the uncertainty search for each parameter so that they include your MLEs (fit4$par)
Uncertainty_FPK(fit=fit4,tree,trait=TRAIT,Npts=25,effort_uncertainty= 100,scope_a=c(-1,20),scope_b=c(-15,5),scope_c=c(-5,5))

# We can fit the OU and BM models using the package geiger to see if likelihoods match with those calculated using BBMV
OU=fitContinuous(phy=tree,dat=TRAIT,model="OU")
OU$opt$lnL ; fit2$lnL # they should be quite close: remember that the FPK model uses an approximation of the likelihood
BM=fitContinuous(phy=tree,dat=TRAIT,model="BM")
BM$opt$lnL

# Ancestral character estimations for internal nodes of the tree can be obtained as follows:
ACE_nodes=ACE_FPK(fit4,specific.point=NULL)
par(mfrow=c(1,1))
plot(ACE_nodes[[90]],type='l') # plot the probability density of the trait at node '90'

# We can also ask for an ACE at any point in the tree. Here we ask for an ACE in the middle of the first branch (as ordered in the 'phylo' object), this is down by specifying (i) the parent node, (ii) the child node, and (iii) the time from the begining of the branch in the 'specific.point' argument:
ACE_1st_branch=ACE_FPK(fit4,specific.point=c(fit4$tree$edge[1,1],fit4$tree$edge[1,2],fit4$tree$edge.length[1]/2))
plot(ACE_1st_branch,type='l')

#####################################
############ BBMV model #############
#####################################

# In the BBMV package, we can also fit a special case of the FPK model in which there are actual (reflective) bounds on the trait interval. We call it the BBMV model, and we need to specify bounds when we create the likelihood function.
# We will use a simulated dataset on which there is a trend towards one of the bounds of the trait interval
par(mfrow=c(1,1))
Vb=3*x # this is the evolutionary potential: it has a linear trend
Vb_norm=exp(-Vb)/sum(exp(-Vb)*step_size) # the step size on the grid
plot(Vb_norm) # a peak on the left bound

TRAITb= Sim_FPK(tree,x0=0,V=Vb,sigma=2,bounds=bounds)
hist(TRAITb,breaks=20) # the distribution of the trait at the tips of the tree: it should reflect the landscape simulated... more or less

# Create four different likelihood functions
ll_BBMV4=lnL_BBMV(tree,TRAITb,Npts=25,bounds=bounds,a=NULL,b=NULL,c=NULL)
ll_BBMV2=lnL_BBMV(tree,TRAITb,Npts=25,bounds=bounds,a=0,b=NULL,c=NULL)
ll_BBMV1=lnL_BBMV(tree,TRAITb,Npts=25,bounds=bounds,a=0,b=0,c=NULL)
ll_BBMV0=lnL_BBMV(tree,TRAITb,Npts=25,bounds=bounds,a=0,b=0,c=0) # this is the BBM model

# fit the four models and plot the macroevolutionary landscapes
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

# compare the fits using AIC
fit4b$aic
fit2b$aic
fit1b$aic # this should be the lowest... at least if we would have had a large enough tree. More complex models (fit4b and fit2b) will also do good job since they can accommodate this trend, but they have more parameters
fit0b$aic

###############################################
####### Measurement error in trait data #######
###############################################
# The FPK and BBMV models can be fit while accounting for measurement error in the value of the trait at the tips of the tree
# This error is not given as a standard-deviation or variance, but rather a vector of multiple measurements of the trait value is given for each tip. The trait object thus becomes a list instead of a vector.
# If you only have an estimate of the standard-error and the mean trait value, then you can provide random draws of the trait for each species, as done below. Here we introduce random measurement error in the values of the TRAIT simulated earlier under the FPK model
TRAIT2=list()
for (i in 1:length(TRAIT)){
  TRAIT2[[i]]=rnorm(n=runif(n=1,min=2,max=10),mean=TRAIT[i],sd=runif(n=1,min=0,max=0.03))
}
names(TRAIT2)=names(TRAIT)
TRAIT2[c(1:4)] # trait measurements for the first four tips: different tips have different number of measures and different levels of error

# Then we can use the same functions as before to fit the model to this dataset. Here we only demonstrate the inclusion of measurement error when fitting the FPK model, but the same works for the BBMV model.
ll_FPK4_with_ME=lnL_FPK(tree,TRAIT2,Npts=25,a=NULL,b=NULL,c=NULL) # the full model
fit4_with_ME=find.mle_FPK(model=ll_FPK4_with_ME)

# Now we can also compare the macroevolutionary landscapes estimated using both methods: they should be rather similar unless measurement error is large relative to the spread of the trait
par(mfrow=c(1,2))
get.landscape.FPK(fit=fit4)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))
get.landscape.FPK(fit=fit4_with_ME)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

# Estimate uncertainty in parameters when fitting the model with measurement error
fit4_with_ME$par
Uncertainty_FPK(fit=fit4_with_ME,tree,trait=TRAIT2,Npts=25,effort_uncertainty= 100,scope_a=c(0,100),scope_b=c(-15,5),scope_c=c(-5,5))

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
# Npts: the number of point on the grid
# pars_init: the initial parameters for starting the algorithm, c(log(sig2/2),a,b,c,x0). Be careful since x0 is actually the point on the grid (between 1 and Npts), not the actual root value
# prob_update: the relative frequencies of update of the different parameters of the model
# verbose: if TRUE, will print some generations of the chain to the screen
# plot: if TRUE, the chain is plotted from time to time
# save_to: file to which the chain is saved (can be useful in case the chain crashes)
# save_every: sets how often the chain is saved
# type_priors: the type of priors used, can be either normal (preferred) or uniform for log(sig2/2), a, b and c, ; and can only be discrete uniform for x0
# shape_priors: list that gives the shape for each prior. (mean,sd) for normal priors and (min,max) for continuous uniform priors. The shape is not specified for the root prior, since it is fixed to be discrete uniform on the grid.
# proposal_type: the type of proposal function, only uniform is available
# proposal_sensitivity: the width of the uniform proposal. The entire value for x0 gives how many steps at a time can be travelled on the trait grid (better to keep it to 1)

source('/Users/florianboucher/Documents/Flo_BACKUPS/Travail/BBM\ plus\ potentiel/BBMV_Github/R/MCMC_function_BBMV.r',chdir=F)
MCMC=MH_MCMC_FPK(tree,trait=TRAIT,bounds=fit4$par_fixed$bounds,Nsteps=10000,record_every=100,plot_every=100,Npts=20,pars_init=c(0,-4,-4,0,1),prob_update=c(0.2,0.25,0.25,0.25,0.05),verbose=TRUE,plot=TRUE,save_to='~/Desktop/MCMC_FPK_test.Rdata',save_every=100,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F)

# Estimate effective sample sizes in our MCMC run using the R package 'coda':
library(coda)
apply(MCMC[-c(1:50),2:9],2,effectiveSize)
par(mfrow=c(2,3))
hist(log(MCMC[-c(1:50),2]/2),breaks=20,main='log(sigsq/2)',ylab=NULL)
hist(MCMC[-c(1:50),6],breaks=20,main='root',ylab=NULL)
plot(1,1)
hist(MCMC[-c(1:50),3],breaks=20,main='a (x^4 term)',ylab=NULL)
hist(MCMC[-c(1:50),4],breaks=20,main='b (x^2 term)',ylab=NULL)
hist(MCMC[-c(1:50),5],breaks=20,main='c (x term)',ylab=NULL)

# We can plot the 95% credible interval of macroevolutionary landscapes estimated in the MCMC run
# The function plots the median value of the macroevolutionary landscape across the posterior in a solid lines and draws a polygon that streches between two quantiles, defined by probs.CI
# Here we set a burnin fraction of 50% since the MCMC run was extremely short
get.landscape.FPK.MCMC(chain=MCMC,bounds=fit4$par_fixed$bounds,Npts=100,burnin=0.5,probs.CI=c(0.025,0.975),COLOR_MEDIAN='red',COLOR_FILL='red',transparency=0.3,main='Macroevolutionary landscapes MCMC',ylab='N.exp(-V)',xlab='Trait',xlim=NULL,ylim=NULL)

# Finally, we can force the potential to be quadratic (i.e. fit an OU model). This is done by fixing the intial values of a to 0 and setting its probability of update to zero 
MCMC_OU=MH_MCMC_FPK(tree,trait=TRAIT,bounds=fit4$par_fixed$bounds,Nsteps=10000,record_every=100,plot_every=100,Npts=20,pars_init=c(0,0,-4,0,1),prob_update=c(0.25,0,0.35,0.35,0.05),verbose=TRUE,plot=TRUE,save_to='~/Desktop/MCMC_FPK_test.Rdata',save_every=100,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F)

apply(MCMC_OU[-c(1:50),2:9],2,effectiveSize)
par(mfrow=c(2,3))
hist(log(MCMC_OU[-c(1:50),2]/2),breaks=20,main='log(sigsq/2)',ylab=NULL)
hist(MCMC_OU[-c(1:50),6],breaks=20,main='root',ylab=NULL)
plot(1,1)
hist(MCMC_OU[-c(1:50),3],breaks=20,main='a (x^4 term)',ylab=NULL)
hist(MCMC_OU[-c(1:50),4],breaks=20,main='b (x^2 term)',ylab=NULL)
hist(MCMC_OU[-c(1:50),5],breaks=20,main='c (x term)',ylab=NULL)

# The MCMC function can also be run with measurement error
MCMC_ME=MH_MCMC_FPK(tree,trait=TRAIT2,bounds=fit4_with_ME$par_fixed$bounds,Nsteps=10000,record_every=100,plot_every=100,Npts=20,pars_init=c(0,-4,-4,0,1),prob_update=c(0.2,0.25,0.25,0.25,0.05),verbose=TRUE,plot=TRUE,save_to='~/Desktop/MCMC_FPK_test.Rdata',save_every=100,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F)
