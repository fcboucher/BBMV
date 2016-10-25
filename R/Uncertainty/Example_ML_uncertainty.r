rm(list=ls())
#Below is an example of the use of BBM+V with simulated data:
library(ape)

# change the path to the following R scripts:
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/BBM+V_functions_MLoptim.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/charac_time.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/MCMC functions BBM+V.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/plot.landscape.BBMV.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/Simulate BBM+V.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/Uncertainty/Uncertainty_BBMV.R', chdir = TRUE)
# Simulate data: tree + continuous trait
library(geiger) # we will use geiger for simulating the tree
tree=sim.bdtree(stop='taxa',n=20) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree to a total depth of 100
SEQ=seq(from=-1.5,to=1.5,length.out=50) # very sensitive
a=0 ; b=-2 ; c=0.
# define V
V=a*SEQ^4+b*SEQ^2+c*SEQ
TRAIT= Sim_BBMV(tree,x0=0,V=V,sigma=1,bounds=c(-5, 5)) # TRAIT simulated on the tree, with a linear trend towards small values (potential increases with high values): for that you need to source the function 'Sim_BBMV.R'
hist(TRAIT,breaks=20) # the distribution of the trait at the tips of the tree: it should be rather left skewed...

###############################################
######## Maximum Likelihood estimation ########
###############################################

# Now try to fit different models:

# Fit the model to the data simulated using maximum-likelihood: fit_BBMV is the main function that does it.
# Multiple starting points are used in the optimization to ensure convergence
# The V_shape parameter determines the shape of potential that you want to fit
# We use only 20 points for discretizing the trait interval to make it faster, but more points should be used on empirical datasets

# fit the model with a flat potential, i.e. BBM:
BBM=fit_BBMV(tree,TRAIT,Npts=20,method='L-BFGS-B',verbose=T,V_shape='flat')
BBM$par # parameters estimated

# fit a model with a linear  potential: 
BBM_x=fit_BBMV(tree,TRAIT,Npts=20,method='L-BFGS-B',verbose=T,V_shape='linear')
BBM_x$par

# fit a model with a quadratic potential:
BBM_x2x=fit_BBMV(tree,TRAIT,Npts=20,method='L-BFGS-B',verbose=T,V_shape='quadratic')
BBM_x2x$par

# fit the most general model with a.x^4+b.x^2+c.x potential:
BBM_full=fit_BBMV(tree,TRAIT,Npts=20,method='L-BFGS-B',verbose=T,V_shape='full')
BBM_full$par

# Fit other classic models of evolution implemented in package {geiger}
BM=fitContinuous(phy=tree,dat=TRAIT,model='BM') # Brownian motion with no bounds
OU=fitContinuous(phy=tree,dat=TRAIT,model='OU') # Ornstein-Uhlenbeck process with a single optimum

# AIC comparison of all the models fitted
BBM$aicc
BBM_x$aicc
BBM_x2x$aicc
BBM_full$aicc
BM$opt$aicc
OU$opt$aicc


# Estimate uncertainty of the different models: seems to work well --> can be merged back to master
Uncertainty_BBMV(BBM,tree,trait= TRAIT,Npts=20,effort_uncertainty= 100)
Uncertainty_BBMV(BBM_x,tree,trait= TRAIT,Npts=20,effort_uncertainty= 100)
Uncertainty_BBMV(BBM_x2x,tree,trait=TRAIT,Npts=20,effort_uncertainty= 100)
Uncertainty_BBMV(BBM_full,tree,trait= TRAIT,Npts=20,effort_uncertainty= 100)

