# This function simulates a trait evolving on a phylogeny under the FPK model 
# It requires the functions 'DiffMat_forward', 'ConvProp_bounds' and 'VectorPos_bounds', all from the BBMV package
# Arguments are:
# tree: a phylogenetic tree in 'phylo' format (ape)
# x0: the value of the trait at the root of the tree, default to 0
# V: a vector describing the potential over the trait interval. The length of V determines the number of traits used to discretize the trait interval. All points are assumed to be equally spaced between the minimum bound (bounds[1]) and the maximum bound (bounds[2]), and V only has values for the potential in each of these points. By default we have have a flat potential, V=rep(0,100), which means BBM with 100 points used for the discretization.
# sigma: the BM evolutionary rate
# the two bounds of the trait interval: x0 must be between them of course
# The function does not simulate step by step, but uses matrix exponential to directly simulate over entire branches of the tree. It is thus very fast.
# The value returned is a vector of trait values, with names matching those of the tips of the tree

Sim_FPK=function(tree,x0=0,V=rep(0,100),sigma,bounds){
	dCoeff=log((sigma)^2/2) # the diffusion coefficient of the process
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



##################################################################
##################################################################
##################################################################
##################################################################
# Below the same simulation function, but which outputs the actual simulation traitgram too
# This is longer since the simulation is done piecewise along each branch


# two functions to transform the true interval to [-1.5,1.5] and back
trans_to_fixed=function(x,bounds){
  return(1.5*(2*(x-bounds[1])/(bounds[2]-bounds[1])-1))
}

trans_from_fixed=function(x,bounds){
  return((bounds[2]+bounds[1])/2+x*(bounds[2]-bounds[1])/3)
}

# Function to simulate trait along given tree under the FPK model
FPK_sim_traitgram=function(tree,x0,a,b,c,bounds,sigsq,time_step,res.x=200,ylim.plot=NULL,return.trait=FALSE){
  # reorder the tree
  tree=reorder.phylo(tree,order="cladewise") # this ordering is nicer for colors in the rainbow palette
  # initiate the plot window
  if (is.null(ylim.plot)){ylim.plot=bounds}
  COL=rainbow(dim(tree$edge)[1])
  plot(x=-max(branching.times(tree)),y=x0,main=NULL,xlim=c(-max(branching.times(tree)),0),ylim=ylim.plot,xlab='Time',ylab='Trait',pch=19,col=COL[1])
  # transform initial point and create a vector for collecting initial traits on branches to be inherited 
  init=rep(NA,(tree$Nnode+length(tree$tip.label))) ; names(init)=seq(from=1,to=(length(tree$tip.label)+tree$Nnode),by=1)
  init[length(tree$tip.label)+1]=trans_to_fixed(x0,bounds) # root trait transformed
  # get branching times of all nodes
  bt=c(rep(0,length(tree$tip.label)),-branching.times(tree)) ; names(bt)=seq(from=1,to=(length(tree$tip.label)+tree$Nnode),by=1)
  # get V and dMat
  x=seq(from=-1.5,to=1.5,length.out=res.x)
  V=a*x^4+b*x^2+c*x # potential
  dMat= DiffMat_forward(V)
  dCoeff=log(sigsq/2)
  pMat=prep_mat_exp(dCoeff,dMat,c(-1.5,1.5))
  # Now simulate along each edge of the tree
  for (i in 1:dim(tree$edge)[1]){
    n.slices=round(tree$edge.length[i]/time_step)
    temp_step=tree$edge.length[i]/n.slices # time step size is slightly modified so that we have an entire number of steps on each branch
    temp_x=rep(NA,(n.slices+1))
    temp_x[1]=init[which(names(init)==tree$edge[i,1])] # the initial trait value on this branch
    for (s in 2:length(temp_x)){ # propagate the FPK process one time step at a time
      ptemp= ConvProp_bounds(X= VectorPos_bounds(temp_x[s-1],V,c(-1.5,1.5)),t=temp_step,pMat)
      temp_x[s]=sample(x,size=1,prob=ptemp/sum(ptemp))
    }
    init[which(names(init)==tree$edge[i,2])]=temp_x[length(temp_x)]
    points(trans_from_fixed(temp_x,bounds)~seq(from=bt[which(names(bt)==tree$edge[i,1])],to=bt[which(names(bt)==tree$edge[i,2])],length.out=length(temp_x)),type='l',col=COL[i])
  }
  if(return.trait==T){
    trait=trans_from_fixed(init[1:length(tree$tip.label)],bounds) ; names(trait)=tree$tip.label
    return(trait)
  }
}
