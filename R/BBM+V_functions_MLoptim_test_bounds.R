# Set of functions to fit the Bounded Brownian Motion model
# Written by F. Boucher & V. Démery, December 1st, 2015

###############################################################################
####### SET OF AUXILIARY FUNCTIONS USED BY THE MASTER FUNCTION BELOW ##########
###############################################################################
# Diagonalize the diffusion matrix that has been discretized
# returns: the diffusion matrix, for simulating traits only
DiffMat_forward=function (V){
	# V is a vector representing the potential, with 'Npts' numeric values
	Npts=length(V)
    M=matrix(0,Npts,Npts)
    for (i in 2:(Npts-1)){
    	M[i-1,i]=exp((V[i]-V[i-1])/2)
    	M[i+1,i]=exp((V[i]-V[i+1])/2)
    	M[i,i]=-(M[i-1,i]+M[i+1,i])
    	    }
    M[2,1]=exp((V[1]-V[2])/2)
    M[1,1]=-M[2,1]
    M[Npts-1,Npts]=exp((V[Npts]-V[Npts-1])/2)
    M[Npts,Npts]=-M[Npts-1,Npts]
    eig=eigen(M)
    return(list(Diff=M,diag=diag(eig$values),passage=eig$vectors))
}  

# new version with t(M) instead since we go backwards in time
DiffMat_backwards=function (V){
	# V is a vector representing the potential, with 'Npts' numeric values
	Npts=length(V)
    M=matrix(0,Npts,Npts)
    for (i in 2:(Npts-1)){
    	M[i-1,i]=exp((V[i]-V[i-1])/2)
    	M[i+1,i]=exp((V[i]-V[i+1])/2)
    	M[i,i]=-(M[i-1,i]+M[i+1,i])
    	    }
    M[2,1]=exp((V[1]-V[2])/2)
    M[1,1]=-M[2,1]
    M[Npts-1,Npts]=exp((V[Npts]-V[Npts-1])/2)
    M[Npts,Npts]=-M[Npts-1,Npts]
    M=t(M) # we go backwards in time!
    eig=eigen(M)
    return(list(Diff=M,diag=diag(eig$values),passage=eig$vectors))
}  

# write to which point of the grid a given position belongs to, 'continuous' version
VectorPos_bounds=function(x,V,bounds){
	Npts=length(V)
	X=rep(0,Npts)
	if (x==bounds[2]){X[Npts]=1}
	else {
		nx=(Npts-1)*(x-bounds[1])/(bounds[2]-bounds[1])
		ix=floor(nx)
		ux=nx-ix
		X[ix+2]=ux
		X[ix+1]=1-ux
	}	
	return(X*(Npts-1)/(bounds[2]-bounds[1]))
}

# Convolution of the matrix x with the propagator for a time t
ConvProp_bounds=function(X,t,dCoeff,dMat,bounds){
	vDiag=dMat$diag ; P=dMat$passage ; tP=solve(P) #; tP=t(P)
	Npts=dim(dMat$diag)[1]
	tau=((bounds[2]-bounds[1])/(Npts-1))^2
	expD=matrix(0,Npts,Npts)
	for (i in 1:Npts){expD[i,i]=exp(exp(dCoeff)*(t/tau)*diag(vDiag)[i])}
	a=P%*%expD%*%tP%*%X
	return(apply(a,1,function(x) max(x,0)))
}

# format tree and trait --> the tree is ordered from tips to root, with edge.length binded to the topology
# A list is also initiated, filled with the position (probabilistic) of each tip and 1 for internal nodes.
FormatTree_bounds=function(tree,trait,V,bounds){
require(ape)	
tree=reorder.phylo(tree,'postorder')
ntips=length(tree$tip.label)
tab=cbind(tree$edge,tree$edge.length) ; colnames(tab)=c('parent','children','brlen')
Pos=list() # one element per node
for (i in 1:(2*ntips-1)){
	if (i>ntips){
		Pos[[i]]=1
	}
	else {
		Pos[[i]]= VectorPos_bounds(trait[tree$tip.label[i]],V,bounds=bounds)
	}
}
return(list(tab=tab,Pos=Pos))
}

# # calculate log-likelihood over the whole tree, to be maximized
# LogLik_bounds=function(tree_formatted,dCoeff,dMat,bounds){
# #tree_formatted obtained through FormatTree ; dCoeff=log(sigsq/2)
# Npts=dim(dMat$diag)[1]
# tree_formatted2= tree_formatted
# logFactor=0
# for (i in 1:dim(tree_formatted2$tab)[1]){
	# tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=dCoeff,dMat=dMat,bounds)
	# norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	# tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
	# logFactor=logFactor+log(norm)
# }
# return(log(max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))+logFactor)
# # this is where we take the max over the root position
# }

# # wrapper for when we only optimize dCoeff
# bBM_loglik_bounds=function(tree_formatted,dMat,bounds){
	# fun=function(dCoeff){
		# return(-LogLik_bounds(tree_formatted,dCoeff,dMat,bounds))
	# }
	# return(fun)
# }

############################
# Same function but when bounds are estimated: SEEMS TO BE OK
# calculate log-likelihood over the whole tree, to be maximized
LogLik_bounds_est=function(tree,trait ,dCoeff,V,bounds){
if ((bounds[1]>min(trait))|(bounds[2]<max(trait))) {return(-Inf)} # bounds have to be outside the trait interval
	else {
Npts_tot=length(V)		
tree_formatted=FormatTree_bounds(tree,trait,V,bounds)		
dMat=DiffMat_backwards(V)
tree_formatted2= tree_formatted
logFactor=0
for (i in 1:dim(tree_formatted2$tab)[1]){
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=dCoeff,dMat=dMat,bounds)
	norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
	logFactor=logFactor+log(norm)
}
return(log(max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))+logFactor)
# this is where we take the max over the root position
	}
}

##################################################
####### LINEAR POTENTIAL + BOUNDS ESTIMATED ######
##################################################
# wrapper for when we optimize dCoeff + linear potential: SEEMS TO WORK 
bBM_loglik_linear_bounds_est=function(tree,trait,Npts){ 
	fun=function(X){ # X=c(dCoeff,a,bmin,bmax) where a is the linear term for the potential. No intercept since we are only interested in V'
		SEQ=seq(from=0,to=1,length.out=Npts)
		V=X[2]*SEQ
		return(-LogLik_bounds_est(tree,trait,X[1],V,c(X[3],X[4])))
	}
	return(fun)
}

Optim_bBM_linear_bounds_est=function(tree,trait,Npts=50,method='Nelder-Mead'){ # SEEMS TO WORK
	# Npts is fixed whatever the bounds...
	# Nelder-Mead is the default method since it seems to be more robust with this rather high number of parameters
    if (!(method%in%c('L-BFGS-B','Nelder-Mead'))){stop('Incorrect optimization method')}
    fun= bBM_loglik_linear_bounds_est(tree,trait,Npts=Npts)
	if (method=='L-BFGS-B'){
	opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='L-BFGS-B',lower=c(-10,-10,min(trait)-10*(max(trait)-min(trait)),max(trait)),upper=c(10,10,min(trait),max(trait)+10*(max(trait)-min(trait))),hessian=FALSE)
	}
    if (method=='Nelder-Mead'){
    opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='Nelder-Mead',hessian=FALSE,control=list(maxit=10000))
    }
    # dCoeff is log(sigma^2/2)
    # now retrieve the ML value of x0, using the ML of dCoeff, bounds and V
    SEQ=seq(from=0,to=1,length.out=Npts)
		V=opt$par[2]*SEQ
		dMat=DiffMat_backwards(V)
		bounds=c(opt$par[3],opt$par[4])
		tree_formatted=FormatTree_bounds(tree,trait,V,bounds)
		tree_formatted2= tree_formatted
    for (i in 1:dim(tree_formatted2$tab)[1]){
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=opt$par[1],dMat=dMat,bounds)
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    }   
    x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
    res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par[1]),a=opt$par[2],root_value=x0),lnL=-opt$value,k=5,aic=2*(5+opt$value),aicc=2*(5+opt$value)+(2*5*(5+1))/(length(trait)-5-1),method=method,convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	return(res)	
}

##################################################
## a.x^4+b.x^2+c.x POTENTIAL + BOUNDS ESTIMATED ##
##################################################
# wrapper for when we optimize dCoeff + x + x2 +x4 potential : SEEMS TO WORK
bBM_loglik_x4x2x_bounds_est=function(tree,trait,Npts){ 
	fun=function(X){ # X=c(dCoeff,a(x4),b(x2),c(x),bmin,bmax) No intercept since we are only interested in V'
		SEQ=seq(from=-1.5,to=1.5,length.out=Npts) # very sensitive
		V=X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ
		return(-LogLik_bounds_est(tree,trait,X[1],V,c(X[5],X[6])))
	}
	return(fun)
}

Optim_bBM_x4x2x_bounds=function(tree,trait,Npts=50,method='Nelder-Mead'){ # SEEMS TO WORK
	# Npts is fixed whatever the bounds...
	# Nelder-Mead is the default method since it seems to be more robust with this rather high number of parameters
    if (!(method%in%c('L-BFGS-B','Nelder-Mead'))){stop('Incorrect optimization method')}
    fun= bBM_loglik_x4x2x_bounds_est(tree,trait,Npts)
	if (method=='L-BFGS-B'){
	opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='L-BFGS-B',lower=c(-10,-10,-10,-10,min(trait)-10*(max(trait)-min(trait)),max(trait)),upper=c(10,10,10,10,min(trait),max(trait)+10*(max(trait)-min(trait))),hessian=FALSE)
	}
    if (method=='Nelder-Mead'){
    opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='Nelder-Mead',hessian=FALSE,control=list(maxit=10000))
    }
    # dCoeff is log(sigma^2/2)
    # now retrieve the ML value of x0, using the ML of dCoeff, bounds and V
    SEQ=seq(from=0,to=1,length.out=Npts)
		V=opt$par[2]*SEQ^4+opt$par[3]*SEQ^2+opt$par[4]*SEQ
		dMat=DiffMat_backwards(V)
		bounds=c(opt$par[5],opt$par[6])
		tree_formatted=FormatTree_bounds(tree,trait,V,bounds)
		tree_formatted2= tree_formatted
    for (i in 1:dim(tree_formatted2$tab)[1]){
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=opt$par[1],dMat=dMat,bounds)
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    }   
    x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
    res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par[1]),a=opt$par[2],b=opt$par[3],c=opt$par[4],root_value=x0),lnL=-opt$value,k=7,aic=2*(7+opt$value),aicc=2*(7+opt$value)+(112)/(length(trait)-8),method=method,convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	return(res)	
}

###############################################################################################
## a.x^4+b.x^2+c.x POTENTIAL + BOUNDS ESTIMATED + NPTS VARIES WITH BOUNDS, GRID STAYS FIXED ###
###############################################################################################
# wrapper for when we optimize dCoeff + x + x2 +x4 potential but grid extends: SEEMS TO WORK
bBM_loglik_x4x2x_flex_pts=function(tree,trait,Npts){ 
	fun=function(X){ # X=c(dCoeff,a(x4),b(x2),c(x),bmin,bmax) No intercept since we are only interested in V'
# Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS
	step=(max(trait)-min(trait))/(Npts-1)
	steps_below=ceiling((min(trait)-X[5])/step) 
	bmin=min(trait)-steps_below*step # we approximate X[5] by the nearest point on the grid, closer to min(trait)
	steps_above=floor((X[6]-max(trait))/step)
	bmax=max(trait)+steps_above*step # we approximate X[6] by the nearest point on the grid, closer to max(trait)
	Npts_tot=Npts+ steps_below+ steps_above
		SEQ=seq(from=-1.5,to=1.5,length.out=Npts_tot) # very sensitive
		V=X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ
		return(-LogLik_bounds_est(tree,trait,X[1],V,c(bmin,bmax)))
	}
	return(fun)
}

Optim_bBM_x4x2x_flex_pts=function(tree,trait,Npts=50,method='Nelder-Mead'){ # SEEMS TO WORK
	# Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS
	# Nelder-Mead is the default method since it seems to be more robust with this rather high number of parameters
    if (!(method%in%c('L-BFGS-B','Nelder-Mead'))){stop('Incorrect optimization method')}
    fun= bBM_loglik_x4x2x_flex_pts(tree,trait,Npts)
	if (method=='L-BFGS-B'){
	opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='L-BFGS-B',lower=c(-10,-10,-10,-10,min(trait)-10*(max(trait)-min(trait)),max(trait)),upper=c(10,10,10,10,min(trait),max(trait)+10*(max(trait)-min(trait))),hessian=FALSE)
	}
    if (method=='Nelder-Mead'){
    opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='Nelder-Mead',hessian=FALSE,control=list(maxit=10000))
    }
    # dCoeff is log(sigma^2/2)
    # now retrieve the ML value of x0, using the ML of dCoeff, bounds and V
    SEQ=seq(from=0,to=1,length.out=Npts)
		V=opt$par[2]*SEQ^4+opt$par[3]*SEQ^2+opt$par[4]*SEQ
		dMat=DiffMat_backwards(V)
		bounds=c(opt$par[5],opt$par[6])
		tree_formatted=FormatTree_bounds(tree,trait,V,bounds)
		tree_formatted2= tree_formatted
    for (i in 1:dim(tree_formatted2$tab)[1]){
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=opt$par[1],dMat=dMat,bounds)
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    }   
    x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
    res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par[1]),a=opt$par[2],b=opt$par[3],c=opt$par[4],root_value=x0),lnL=-opt$value,k=7,aic=2*(7+opt$value),aicc=2*(7+opt$value)+(112)/(length(trait)-8),method=method,convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	return(res)	
}

###############################################################################################
## a.x^4+b.x^2+c.x POTENTIAL + BOUNDS EST. + NPTS VARIES WITH BOUNDS, GRID CHANGES SLIGHTLY ###
###############################################################################################
# wrapper for when we optimize dCoeff + x + x2 +x4 potential : SEEMS TO WORK
bBM_loglik_x4x2x_flex_grid=function(tree,trait,Npts){ 
		fun=function(X){ # X=c(dCoeff,a(x4),b(x2),c(x),bmin,bmax) No intercept since we are only interested in V'
		# Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS BUT THE GRID CHANGES SLIGHTLY TO ACCOMODATE THE TRUE BOUNDS PROPOSED
		step=(max(trait)-min(trait))/(Npts-1)
		steps_below=ceiling((min(trait)-X[5])/step) 
		steps_above=floor((X[6]-max(trait))/step)
		Npts_tot=Npts+ steps_below+ steps_above
		SEQ=seq(from=-1.5,to=1.5,length.out=Npts_tot) # very sensitive
		V=X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ
		return(-LogLik_bounds_est(tree,trait,X[1],V,c(X[5],X[6])))
	}
	return(fun)
}

Optim_bBM_x4x2x_flex_grid=function(tree,trait,Npts=50,method='Nelder-Mead'){ # SEEMS TO WORK
	# Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS
	# Nelder-Mead is the default method since it seems to be more robust with this rather high number of parameters
    if (!(method%in%c('L-BFGS-B','Nelder-Mead','SANN'))){stop('Incorrect optimization method')}
    fun= bBM_loglik_x4x2x_flex_grid(tree,trait,Npts)
	if (method=='L-BFGS-B'){
	opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='L-BFGS-B',lower=c(-10,-10,-10,-10,min(trait)-10*(max(trait)-min(trait)),max(trait)),upper=c(10,10,10,10,min(trait),max(trait)+10*(max(trait)-min(trait))),hessian=FALSE)
	}
    if (method=='Nelder-Mead'){
    opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='Nelder-Mead',hessian=FALSE,control=list(maxit=10000))
    }
    if (method=='SANN'){
    opt=optim(par=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))),fn=fun,method='SANN',hessian=FALSE,control=list(maxit=10000))
    }

    # dCoeff is log(sigma^2/2)
    # now retrieve the ML value of x0, using the ML of dCoeff, bounds and V
    SEQ=seq(from=0,to=1,length.out=Npts)
		V=opt$par[2]*SEQ^4+opt$par[3]*SEQ^2+opt$par[4]*SEQ
		dMat=DiffMat_backwards(V)
		bounds=c(opt$par[5],opt$par[6])
		tree_formatted=FormatTree_bounds(tree,trait,V,bounds)
		tree_formatted2= tree_formatted
    for (i in 1:dim(tree_formatted2$tab)[1]){
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=opt$par[1],dMat=dMat,bounds)
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    }   
    x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
    res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par[1]),a=opt$par[2],b=opt$par[3],c=opt$par[4],root_value=x0),lnL=-opt$value,k=7,aic=2*(7+opt$value),aicc=2*(7+opt$value)+(112)/(length(trait)-8),method=method,convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	return(res)	
}

###############################################################################################
## a.x^4+b.x^2+c.x POTENTIAL + BOUNDS ESTIMATED + NPTS VARIES WITH BOUNDS, GRID STAYS FIXED ###
###############################################################################################
Optim_bBM_x4x2x_flex_pts_start=function(tree,trait,Npts=50,method='Nelder-Mead',start.point=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait)))){ # SEEMS TO WORK
	# Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS
	# Nelder-Mead is the default method since it seems to be more robust with this rather high number of parameters
    if (!(method%in%c('L-BFGS-B','Nelder-Mead'))){stop('Incorrect optimization method')}
    fun= bBM_loglik_x4x2x_flex_pts(tree,trait,Npts)
	if (method=='L-BFGS-B'){
	opt=optim(par=start.point,fn=fun,method='L-BFGS-B',lower=c(-10,-10,-10,-10,min(trait)-10*(max(trait)-min(trait)),max(trait)),upper=c(10,10,10,10,min(trait),max(trait)+10*(max(trait)-min(trait))),hessian=FALSE)
	}
    if (method=='Nelder-Mead'){
    opt=optim(par=start.point,fn=fun,method='Nelder-Mead',hessian=FALSE,control=list(maxit=10000))
    }
    # dCoeff is log(sigma^2/2)
    # now retrieve the ML value of x0, using the ML of dCoeff, bounds and V
    SEQ=seq(from=0,to=1,length.out=Npts)
		V=opt$par[2]*SEQ^4+opt$par[3]*SEQ^2+opt$par[4]*SEQ
		dMat=DiffMat_backwards(V)
		bounds=c(opt$par[5],opt$par[6])
		tree_formatted=FormatTree_bounds(tree,trait,V,bounds)
		tree_formatted2= tree_formatted
    for (i in 1:dim(tree_formatted2$tab)[1]){
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=opt$par[1],dMat=dMat,bounds)
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    }   
    x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
    res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par[1]),a=opt$par[2],b=opt$par[3],c=opt$par[4],root_value=x0),lnL=-opt$value,k=7,aic=2*(7+opt$value),aicc=2*(7+opt$value)+(112)/(length(trait)-8),method=method,convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	return(res)	
	
}

###############################################################################################
### a.x^4+b.x^2+c.x + BOUNDS ESTIMATED + NPTS VARIES, GRID STAYS FIXED, MULTIPLE START PTS ####
###############################################################################################
Optim_bBM_x4x2x_flex_pts_multiple_starts=function(tree,trait,Npts=50,method='Nelder-Mead',verbose=T){
	pars=matrix(c(rep(log(var(trait)/(2*max(branching.times(tree)))),3),rep(1,3),rep(c(0,0.1,0.5),2)),nrow=6,ncol=2) # different starting points for sigsq and bounds
	starts=list()
	for (i in 1:dim(pars)[1]){
		starts[[i]]=try(Optim_bBM_x4x2x_flex_pts_start(tree, trait,Npts= Npts,method=method,start.point=c(pars[i,1],0,0,0,min(trait)-pars[i,2]*(max(trait)-min(trait)),max(trait)+ pars[i,2]*(max(trait)-min(trait)))))
		if (class(starts[[i]])=='try-error') {starts[[i]]=list(lnL=-Inf)}
		if (verbose==T){cat("Finished preliminary fit, log-lik= ",starts[[i]]$lnL,sep="\n")}
	}
	# add one fit with bounds fixed?
	starts[[(length(starts)+1)]]=try(Optim_bBM_x4x2x(tree, trait,Npts= Npts,bounds=c(min(trait),max(trait)),method=method))
	if (class(starts[[length(starts)]])=='try-error') {starts[[length(starts)]]=list(lnL=-Inf)}
	if (verbose==T){cat("Finished preliminary fit with bounds fixed to min/max, log-lik= ",starts[[length(starts)]]$lnL,sep="\n")}
	lnls=c(starts[[1]]$lnL)
	for (i in 2:length(starts)){lnls=c(lnls,starts[[i]]$lnL)}
	mod=starts[[which(lnls==max(lnls))]] # the starting point which got us to the highest lnl
	# check for convergence before returning ML estimate
	if (verbose==T){cat("Starting final fit... ",sep="\n")}
	temp=Optim_bBM_x4x2x_flex_pts_start(tree, trait,Npts= Npts,method=method,start.point=c(log(mod$par$sigsq/2),mod$par$a,mod$par$b,mod$par$c,mod$par$bounds[1],mod$par$bounds[2]))
	rerun=0
	while((temp$convergence!=0)&(rerun<11)){
		if (verbose==T){cat("Trying to improve convergence... ",sep="\n")
			if (rerun==10){cat("...this is the last try!",sep="\n")}}
		temp=Optim_bBM_x4x2x_flex_pts_start(tree, trait,Npts= Npts,method= method,start.point=c(log(temp$par$sigsq/2),temp$par$a, temp$par$b, temp$par$c, temp$par$bounds[1], temp$par$bounds[2]))
		rerun=rerun+1
	}
	return(temp)
}