######################################################
# FPK: the model with no bounds (i.e. bounds far away)
lnL_FPK=function(tree,trait,a=NULL,b=NULL,c=NULL,Npts){
  if (is.numeric(trait)){
    bounds=c(min(trait)-(max(trait)-min(trait))/2,max(trait)+(max(trait)-min(trait))/2)}
  if (class(trait)=='list') {
    bounds=c(min(unlist(trait))-(max(unlist(trait))-min(unlist(trait)))/2,max(unlist(trait))+(max(unlist(trait))-min(unlist(trait)))/2)
  }
  if (sum(tree$tip.label%in%names(trait))<max(length(trait),length(tree$tip.label))){stop('Tip names in tree do not match names of the trait vector')}
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  tree_formatted= FormatTree_bounds(tree,trait,V=rep(0,Npts),bounds=bounds)
  ncoeff=(is.null(a)==T)+(is.null(b)==T)+(is.null(c)==T)
  if (is.null(a)==F){
    if (is.null(b)==F){
      # all three shape parameters fixed (e.g. flat landscape if a=b=c=0)  
      if (is.null(c)==F){fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds))}}
      # only c varies (e.g. flat landscape if a=b=0)  
      else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[2]*SEQ),bounds=bounds))}}
    }
    # a is fixed (e.g. quadratic landscape if a=0)  
    else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+X[2]*SEQ^2+X[3]*SEQ),bounds=bounds))}}
  }
  # the full model: no parameter fixed
  else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ),bounds=bounds))}}
  return(list(fun=fun,ncoeff=ncoeff,par_fixed=list(a=a,b=b,c=c,bounds=bounds),tree=tree,trait=trait,Npts=Npts))
}

#################################################
# BBMV: the model with bounds defined by the user
lnL_BBMV=function(tree,trait,bounds,a=NULL,b=NULL,c=NULL,Npts){
  if (sum(tree$tip.label%in%names(trait))<max(length(trait),length(tree$tip.label))){stop('Tip names in tree do not match names of the trait vector')}
  if (is.numeric(trait)){
    if ((min(trait)<bounds[1])|(max(trait)>bounds[2])){stop('Some values in the trait vector exceed the bounds.')} 
  }
  if (class(trait)=='list') {
    if ((min(unlist(trait))<bounds[1])|(max(unlist(trait))>bounds[2])){stop('Some values in the trait data exceed the bounds.')}
  }
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  tree_formatted= FormatTree_bounds(tree,trait,V=rep(0,Npts),bounds=bounds)
  ncoeff=(is.null(a)==T)+(is.null(b)==T)+(is.null(c)==T)
  if (is.null(a)==F){
    if (is.null(b)==F){
      # all three shape parameters fixed (e.g. flat landscape if a=b=c=0)  
      if (is.null(c)==F){fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds))}}
      # only c varies (e.g. flat landscape if a=b=0)  
      else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[2]*SEQ),bounds=bounds))}}
    }
    # a is fixed (e.g. quadratic landscape if a=0)  
    else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+X[2]*SEQ^2+X[3]*SEQ),bounds=bounds))}}
  }
  # the full model: no parameter fixed
  else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ),bounds=bounds))}}
  return(list(fun=fun,ncoeff=ncoeff,par_fixed=list(a=a,b=b,c=c,bounds=bounds),tree=tree,trait=trait,Npts=Npts))
}

#############################################################
# function to fit the models prepared by lnL_FPK and lnL_BBMV
find.mle_FPK=function(model,method='Nelder-Mead',init.optim=NULL,safe=F){
  # safe=T for safer optimization starting from 3 different starting points
  if (safe==F){ # only one optimization
    if(is.null(init.optim)==T){init.optim=c(-10,rep(0,model$ncoeff))}
    else{}
    opt=optim(par=init.optim,fn=model$fun,method=method,control=list(maxit=50000))
  }
  else {
    init=c(-10,-1,0) #diffusion Coefficient
    starts=list()
    for (i in 1:3){
      starts[[i]]=optim(par=c(init[i],rep(0,model$ncoeff)),fn=model$fun,method=method,control=list(maxit=50000))
      cat("Finished preliminary fit, log-lik= ",-starts[[i]]$value,sep="\n")
    }
    lnls=c(-starts[[1]]$value)
    for (i in 2:length(starts)){lnls=c(lnls,-starts[[i]]$value)}
    opt=starts[[which(lnls==max(lnls))]] # the starting point which gave the highest lnl
  }
  par_fixed=model$par_fixed
  par=list()
  par$sigsq=2*exp(opt$par[1])
  if (model$ncoeff==3){par$a=opt$par[2] ; par$b=opt$par[3] ; par$c=opt$par[4]}
  else if (model$ncoeff==2){par$b=opt$par[2] ; par$c=opt$par[3]}
  else if (model$ncoeff==1){par$c=opt$par[2]}
  else{}
  # now retrieve x0
  # retrieve coefficients for the potential
  if ('a'%in%names(par)){a=par$a}
  else {a=par_fixed$a}
  if ('b'%in%names(par)){b=par$b}
  else {b=par_fixed$b}
  if ('c'%in%names(par)){c=par$c}
  else {c=par_fixed$c}
  bounds=par_fixed$bounds
  Npts=model$Npts ; tree=model$tree ; trait=model$trait
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  V=a*SEQ^4+b*SEQ^2+c*SEQ
  dMat=DiffMat_backwards(V)
  ll_root=as.data.frame(matrix(NA, Npts,2))
  colnames(ll_root)=c('root','density')
  ll_root$root =seq(bounds[1],to=bounds[2], length.out= Npts)
  tree_formatted= FormatTree_bounds(tree,trait,V=rep(0, Npts),bounds=bounds)
  tree_formatted2= tree_formatted
  pMat=prep_mat_exp(dCoeff=log(par$sigsq/2),dMat,bounds=bounds) # edited
  logFactor=0
  for (i in 1:dim(tree_formatted2$tab)[1]){
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat)
    norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
    logFactor=logFactor+log(norm)
  }
  ll_root$density =tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]
  return(res=list(lnL=-opt$value,aic=2*(length(init.optim)+opt$value),k=length(init.optim),par=par,par_fixed=par_fixed,root=ll_root,convergence=opt$convergence,message=opt$message,tree=tree,trait=trait,Npts=Npts)) 
}
