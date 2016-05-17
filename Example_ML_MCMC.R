rm(list=ls())
#Below is an example of the use of BBM+V with simulated data:
library(ape)

# change the path to the following R scripts:
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/BBM+V_functions_MLoptim.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/charac_time.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/MCMC functions BBM+V.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/plot.landscape.BBMV.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/Simulate BBM+V.R', chdir = TRUE)

# Simulate data: tree + continuous trait
library(geiger) # geiger is needed for simulating the tree
tree=sim.bdtree(stop='taxa',n=20) # tree with few tips for fast tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree to a total depth of 100
TRAIT= Sim_BBMV(tree,x0=0,V=seq(from=0,to=5,length.out=50),sigma=10,bounds=c(-5, 5)) # TRAIT simulated on the tree, with a linear trend towards small values (potential increases with high values): for that you need to source the function 'Sim_BBMV.R'
hist(TRAIT,breaks=20) # the distribution of the trait at the tips of the tree: it should be rather right skewed...

###############################################
######## Maximum Likelihood estimation ########
###############################################

# Now try to fit different models:

# Fit the model to the data simulated: we use only 10 points for discretizing the trait interval to make it faster, but more points should be used on empirical datasets

# fit the model with a flat potential, i.e. BBM: multiple starting points in the optimization to ensure convergence
BBM=Optim_bBM_0_flex_pts_multiple_starts(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T) 
BBM$par # parameters estimated

# fit a model with a linear  potential: multiple starting points in the optimization to ensure convergence
BBM_x=Optim_bBM_x_flex_pts_multiple_starts(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T) # ML value of bounds estimated along other parameters
BBM_x$par

# fit a model with a quadratic potential: multiple starting points in the optimization to ensure convergence
BBM_x2x=Optim_bBM_x2x_flex_pts_multiple_starts(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T) # ML value of bounds estimated along other parameters
BBM_x2x$par

# fit the most general model with a.x^4+b.x^2+c.x potential: multiple starting points in the optimization to ensure convergence
BBM_full= Optim_bBM_x4x2x_flex_pts_multiple_starts(tree,TRAIT,Npts=20,method='Nelder-Mead',verbose=T) # ML value of bounds estimated along other parameters
BBM_full$par


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

# Now plot adaptive landscape estimated by best model
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
# You need to specify the file to which the chain is save ('save_to' parameter)
# Here we will do only a few generations so that computation time is not too long

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
