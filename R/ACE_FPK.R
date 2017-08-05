ACE_FPK <-
  function(fit,specific.point=NULL){
  # the default settings (specific.point=NULL) will produce an ACE at all internal nodes in the tree
  # Alternatively, specific.point can be used to ask for an estimation of the trait probability density at any specific point in the tree (i.e., not a node)
  # specific point must then be a vector with three elements: c(parent_node,child_node,time_from_start_of branch)
  # retrieve coefficients for the potential
  if ('a'%in%names(fit$par)){a=fit$par$a}
  else {a=fit$par_fixed$a}
  if ('b'%in%names(fit$par)){b=fit$par$b}
  else {b=fit$par_fixed$b}
  if ('c'%in%names(fit$par)){c=fit$par$c}
  else {c=fit$par_fixed$c}
  Npts=fit$Npts
  # build potential and stationary distribution of the trait
  bounds=fit$par_fixed$bounds
  sigsq=fit$par$sigsq
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
  step=(bounds[2]-bounds[1])/(Npts-1)
  tree=fit$tree ; trait=fit$trait
  dMat=DiffMat_backwards(V)
  tree_formatted= FormatTree_bounds(tree,trait,V=rep(0, Npts),bounds=bounds)
  tree_formatted2= tree_formatted
  pMat=prep_mat_exp(dCoeff=log(sigsq/2),dMat,bounds=bounds) # edited
  logFactor=0
  for (i in 1:dim(tree_formatted2$tab)[1]){
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat)
    norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
    logFactor=logFactor+log(norm)
  }
  if (is.null(specific.point)){
  ACE=lapply(tree_formatted2$Pos,function(x){cbind(seq(from=bounds[1],to=bounds[2],length.out=length(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])),x)})
  }
  else { # density at a specific point on a branch of the tree
    branch=which((tree_formatted2$tab[,1]==specific.point[1])&(tree_formatted2$tab[,2]==specific.point[2]))
    time_back=tree_formatted2$tab[branch,3]-specific.point[3]
    if (time_back<0){
      stop("Problem with your specification of specific.point: the time on the branch for which you requested an ACE is longer than the branch itself")
    }
    else {
      Conv=ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[branch,2]]],t=time_back,prep_mat = pMat)
     Conv=Conv/sum(Conv)
     ACE=cbind(seq(from=bounds[1],to=bounds[2],length.out=length(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])),Conv)
    }
  }
  return(ACE)
}
