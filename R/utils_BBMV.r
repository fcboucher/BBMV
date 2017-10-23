###############################################################################
############### SET OF AUXILIARY FUNCTIONS USED TO FIT THE MODELS #############
###############################################################################
# Create and diagonalize the transition matrix that has been discretized
# returns: the transition matrix going forward in time, for simulating traits only
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
  passage=matrix(NA,dim(eig$vectors)[1],dim(eig$vectors)[2])
  for (col in 1:dim(eig$vectors)[2]){passage[,col]=Re(eig$vectors[,col])}
  return(list(Diff=M,diag=diag(Re(eig$values)),passage=passage))
}  

# Create and diagonalize the transition matrix that has been discretized
# returns: the transition matrix going backwards in time, used for inference
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
  passage=matrix(NA,dim(eig$vectors)[1],dim(eig$vectors)[2])
  for (col in 1:dim(eig$vectors)[2]){passage[,col]=Re(eig$vectors[,col])}
  return(list(Diff=M,diag=diag(Re(eig$values)),passage=passage))
  }  

# write to which point of the grid a given position belongs to, 'continuous' version
VectorPos_bounds=function(x,V,bounds){
  Npts=length(V)
  if (length(x)==1){ # only one value per tip
    X=rep(0,Npts)  
    if (x==bounds[2]){X[Npts]=1}
    else {
      nx=(Npts-1)*(x-bounds[1])/(bounds[2]-bounds[1])
      ix=floor(nx)
      ux=nx-ix
      X[ix+2]=ux
      X[ix+1]=1-ux
    }	
  }
  else {
    # here we treat the case in which we do not have a single value but a vector of values measured
    MAT=matrix(0,length(x),Npts)
    for (i in 1:length(x)){ # problem with the recursive form here
      if (x[i]==bounds[2]){MAT[i,Npts]=1}
      else {
        nx=(Npts-1)*(x[i]-bounds[1])/(bounds[2]-bounds[1])
        ix=floor(nx)
        ux=nx-ix
        MAT[i,ix+2]=ux
        MAT[i,ix+1]=1-ux
      }	
    }
    X=apply(MAT,2,mean) # and average them
  }
  return(X*(Npts-1)/(bounds[2]-bounds[1]))
}

# Prepare the matrix diagonal
prep_mat_exp=function(dCoeff,dMat,bounds){
  vDiag=dMat$diag ; P=dMat$passage ; tP=solve(P,tol = 1e-30) #; tP=t(P)
  Npts=dim(dMat$diag)[1]
  tau=((bounds[2]-bounds[1])/(Npts-1))^2
  diag_expD=exp(dCoeff)/tau*diag(vDiag) # faster than the for loop
  return(list(P=P,tP=tP,diag_expD=diag_expD))
}

# Convolution product over one branch
ConvProp_bounds=function(X,t,prep_mat){
  Npts=length(X)
  expD=matrix(0,Npts,Npts)
  diag(expD)=exp(t*prep_mat$diag_expD)
  a=prep_mat$P%*%expD%*%prep_mat$tP%*%X
  return(apply(a,1,function(x) max(x,0))) # prevent rounding errors for small numbers
}

# format tree and trait --> the tree is ordered from tips to root, with edge.length binded to the topology
# A list is also initiated, filled with the position (probabilistic) of each tip and 1 for internal nodes.
FormatTree_bounds=function(tree,trait,V,bounds){
  tree=reorder.phylo(tree,'postorder')
  ntips=length(tree$tip.label)
  tab=cbind(tree$edge,tree$edge.length) ; colnames(tab)=c('parent','children','brlen')
  Pos=list() # one element per node
  for (i in 1:(2*ntips-1)){
    if (i>ntips){
      Pos[[i]]=1
    }
    else {
      if (class(trait)=='numeric'){Pos[[i]]= VectorPos_bounds(trait[tree$tip.label[i]],V=V,bounds=bounds)} # only one value per tip
      else {Pos[[i]]= VectorPos_bounds(trait[[tree$tip.label[i]]],V=V,bounds=bounds)} # multiple values per tip (i.e. uncertainty)
    }
  }
  return(list(tab=tab,Pos=Pos))
}

############################
# calculate log-likelihood over the whole tree, to be maximized
LogLik_bounds=function(tree_formatted,dCoeff,dMat,bounds){
  #tree_formatted obtained through FormatTree ; dCoeff=log(sigsq/2)
  Npts=dim(dMat$diag)[1]
  tree_formatted2= tree_formatted
  pMat=prep_mat_exp(dCoeff,dMat,bounds) # edited
  logFactor=0
  for (i in 1:dim(tree_formatted2$tab)[1]){
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat)
    norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
    logFactor=logFactor+log(norm)
  }
  return(log(max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))+logFactor)
  # this is where we take the max over the root position
}