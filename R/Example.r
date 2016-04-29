rm(list=ls()) # for tests only

#Below is an example of the use of BBM+V with simulated data:
library(ape)

source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/Cleaned functions for release/Simulate BBM+V.R', chdir = TRUE) # change the path to this script
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/Cleaned functions for release/BBM+V_functions_MLoptim.R', chdir = TRUE) # change the path to this script
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/Cleaned functions for release/plot.landscape.BBMV.R', chdir = TRUE)
source('~/Documents/Flo_BACKUPS/Travail/BBM plus potentiel/BBMV_Github/R/charac_time.R', chdir = TRUE)

# Simulate data: tree + continuous trait
library(geiger) # geiger is needed for simulating the tree
tree=sim.bdtree(stop='taxa',n=20) # tree with few tips for fast tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree to a total depth of 100
TRAIT= Sim_BBMV(tree,x0=0,V=seq(from=0,to=5,length.out=50),sigma=1,bounds=c(-5, 5)) # TRAIT simulated on the tree, with a linear trend towards small values (potential increases with high values): for that you need to source the function 'Sim_BBMV.R'
hist(TRAIT,breaks=20) # the distribution of the trait at the tips of the tree: it should be rather right skewed...


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