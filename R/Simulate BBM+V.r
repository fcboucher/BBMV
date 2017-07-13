# This function simulates a trait evolving on a phylogeny under the BBM+V model (i.e. BM between two reflective bound plus a force acting on the trait)
# It requires the functions 'DiffMat_forward', 'ConvProp_bounds' and 'VectorPos_bounds', all from the BBMV package
# Arguments are:
# tree: a phylogenetic tree in 'phylo' format (ape)
# x0: the value of the trait at the root of the tree, default to 0
# V: a vector describing the potential over the trait interval. The length of V determines the number of traits used to discretize the trait interval. All points are assumed to be equally spaced between the minimum bound (bounds[1]) and the maximum bound (bounds[2]), and V only has values for the potential in each of these points. By default we have have a flat potential, V=rep(0,100), which means BBM with 100 points used for the discretization.
# sigma: the BM evolutionary rates
# the two bounds of the trait interval: x0 must be between them of course
# The function does not simulate step by step, but uses matrix exponential to directly simulate over entire branches of the tree. It is thus very fast.
# The value returned is a vector of trait values, with names matching those of the tips of the tree

Sim_BBMV=function(tree,x0=0,V=rep(0,100),sigma,bounds){
	dCoeff=log((sigma)^2/2) # the coefficient of diffusion of the model
	dMat= DiffMat_forward(V) # the transition matrix describing the probablity of evolving between two sites in the trait grid in an infinitesimal time step.
	Npts=length(V)
	ntips=length(tree$tip.label)
	trait=rep(NA,2*ntips-1) ; names(trait)=1:(2*ntips-1)
	trait[ntips+1]=x0  # root
	pMat=prep_mat_exp(dCoeff=dCoeff,dMat,bounds) # edited
for (i in 1:length(tree$edge.length)){
	proptemp= ConvProp_bounds(X= VectorPos_bounds(trait[tree$edge[i,1]],V,bounds),t=tree$edge.length[i],prep_mat = pMat) # propagate the trait forward in time: EDITED
	trait[tree$edge[i,2]]=sample(x=seq(from=bounds[1],to=bounds[2],length.out=Npts),size=1,prob= proptemp/sum(proptemp))	# sample from this probability distribution	to get a descendent node value
}
TRAIT=trait[1:ntips] ; names(TRAIT)=tree$tip.label
	return(TRAIT)
}