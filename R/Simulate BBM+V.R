Sim_BBMV=function(tree,x0,V=rep(0,100),sigma,bounds){
	dCoeff=log((sigma)^2/2) # the coefficient of diffusion of the model
	dMat= DiffMat_forward(V)
	Npts=length(V)
	ntips=length(tree$tip.label)
	trait=rep(NA,2*ntips-1) ; names(trait)=1:(2*ntips-1)
	trait[ntips+1]=x0  # root
for (i in 1:length(tree$edge.length)){
	proptemp= ConvProp_bounds(X= VectorPos_bounds(trait[tree$edge[i,1]],V,bounds),t=tree$edge.length[i],dCoeff,dMat,bounds) # propagate the trait forward in time
	trait[tree$edge[i,2]]=sample(x=seq(from=bounds[1],to=bounds[2],length.out=Npts),size=1,prob= proptemp/sum(proptemp))	# sample from this probability distribution	to get a descendent node value
}
TRAIT=trait[1:ntips] ; names(TRAIT)=tree$tip.label
	return(TRAIT)
}