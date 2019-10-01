
####################################################################################################
##### Likelihood and optimization function for BBMV when  V IS SHARED between trees BUT NOT SIGMA ##
####################################################################################################
lnl_BBMV_multiclades_same_V_different_sig2=function(trees,traits,bounds,a=NULL,b=NULL,c=NULL,Npts=50){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  for (i in 1:length(trees)){
    if (sum(trees[[i]]$tip.label%in%names(traits[[i]]))<max(length(traits[[i]]),length(trees[[i]]$tip.label))){stop(paste('Tip names in tree ',i,' do not match names of corresponding trait vector'))}
    #if ((min(traits[[i]])<bounds[1])|(max(traits[[i]])>bounds[2])){stop(paste('Some values in trait vector ',i,' exceed the bounds.'))}
###### new piece of code added for Measurment error incorporation
    if (is.numeric(traits[[i]])){
      if ((min(traits[[i]])<bounds[1])|(max(traits[[i]])>bounds[2])){stop(paste('Some values in trait vector',i,'exceed the bounds.'))} 
    }
    if (class(traits[[i]])=='list') {
      if ((min(unlist(traits[[i]]))<bounds[1])|(max(unlist(traits[[i]]))>bounds[2])){stop(paste('Some values in trait vector',i,'exceed the bounds.'))}
    }
######  end new code  
  }
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  trees_formatted=list()
  for (i in 1:length(trees)){
    trees_formatted[[i]]=FormatTree_bounds(trees[[i]],traits[[i]],V=rep(0,Npts),bounds=bounds)
  }
  ncoeff=(is.null(a)==T)+(is.null(b)==T)+(is.null(c)==T) # OK
  npar=ncoeff+length(trees)
  par_names=c() 
  for (i in 1:length(trees)){par_names=c(par_names,paste('dCoeff_tree_',i,sep=''))}
  par_names=c(par_names,'a','b','c')
  if (is.null(a)==F){
    if (is.null(b)==F){
      # all three shape parameters fixed (e.g. flat landscape if a=b=c=0): seems to work 
      if (is.null(c)==F){
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }      
 
      # only c varies (e.g. flat landscape if a=b=0): seems to work   
      else {
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[',length(trees)+1,']*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }
    }
    # only a is fixed (e.g. quadratic landscape if a=0): seems to work 
    else {
      fun_text='fun=function(X){return('
      for (i in 1:length(trees)){
        fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+X[',length(trees)+1,']*SEQ^2+X[',length(trees)+2,']*SEQ),bounds=bounds)',sep='') 
      }
      fun_text=paste(fun_text,')}',sep='') # the end parenthesis
      fun=eval(parse(text=fun_text))
    } 
  }
  
  # the full model: no parameter fixed
  else {
    fun_text='fun=function(X){return('
    for (i in 1:length(trees)){
      fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(X[',length(trees)+1,']*SEQ^4+X[',length(trees)+2,']*SEQ^2+X[',length(trees)+3,']*SEQ),bounds=bounds)',sep='') 
    }
    fun_text=paste(fun_text,')}',sep='') # the end parenthesis
    fun=eval(parse(text=fun_text))
  }
  return(list(fun=fun,ncoeff=ncoeff,npar=npar,par_names=par_names,par_fixed=list(a=a,b=b,c=c,bounds=bounds),trees=trees,traits=traits,Npts=Npts))
}

####################################################################################
# Optimization function: seems to work so far
find.mle_FPK_multiple_clades_same_V_different_sig2=function(model,method='Nelder-Mead',init.optim=NULL){
  if(is.null(init.optim)==T){init.optim=c(rep(-5,length(model$trees)),rep(0,model$ncoeff))}
  else{}
  opt=optim(par=init.optim,fn=model$fun,method=method,control=list(maxit=5000))
  par_fixed=model$par_fixed
  par_names=c() 
  for (i in 1:length(model$trees)){par_names=c(par_names,paste('sigsq_tree_',i,sep=''))}
  par_names=c(par_names,'a','b','c')
  par=list()
  for (i in 1:length(model$trees)){par[[i]]=2*exp(opt$par[i]) ; names(par)[i]=par_names[i]}
  if (model$ncoeff==3){par$a=opt$par[(length(model$trees)+1)] ; par$b=opt$par[(length(model$trees)+2)] ; par$c=opt$par[(length(model$trees)+3)]}
  else if (model$ncoeff==2){par$b=opt$par[(length(model$trees)+1)] ; par$c=opt$par[(length(model$trees)+2)]}
  else if (model$ncoeff==1){par$c=opt$par[(length(model$trees)+1)]}
  else{}
  # now retrieve x0: not needed actually, since this could be done using the ACE function
  return(res=list(lnL=-opt$value,aic=2*(length(init.optim)+opt$value),k=length(init.optim),par=par,par_fixed=par_fixed,convergence=opt$convergence,message=opt$message,trees=model$trees,traits=model$traits,Npts=model$Npts)) 
}

####################################################################################################
###### Likelihood and optimization function for BBMV when ALL PARAMS ARE SHARED between trees ######
####################################################################################################
lnl_BBMV_multiclades_same_V_same_sig2=function(trees,traits,bounds,a=NULL,b=NULL,c=NULL,Npts=50){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  for (i in 1:length(trees)){
  if (sum(trees[[i]]$tip.label%in%names(traits[[i]]))<max(length(traits[[i]]),length(trees[[i]]$tip.label))){stop(paste('Tip names in tree ',i,' do not match names of corresponding trait vector'))}
  #if ((min(traits[[i]])<bounds[1])|(max(traits[[i]])>bounds[2])){stop(paste('Some values in trait vector ',i,' exceed the bounds.'))}
###### new piece of code added for Measurment error incorporation
    if (is.numeric(traits[[i]])){
      if ((min(traits[[i]])<bounds[1])|(max(traits[[i]])>bounds[2])){stop(paste('Some values in trait vector',i,'exceed the bounds.'))} 
    }
    if (class(traits[[i]])=='list') {
      if ((min(unlist(traits[[i]]))<bounds[1])|(max(unlist(traits[[i]]))>bounds[2])){stop(paste('Some values in trait vector',i,'exceed the bounds.'))}
    }
######  end new code     
  }
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  trees_formatted=list()
  for (i in 1:length(trees)){
    trees_formatted[[i]]=FormatTree_bounds(trees[[i]],traits[[i]],V=rep(0,Npts),bounds=bounds)
  }
  ncoeff=(is.null(a)==T)+(is.null(b)==T)+(is.null(c)==T) # OK, but to edit if we want some params to vary between clades
  if (is.null(a)==F){
    if (is.null(b)==F){
      # all three shape parameters fixed (e.g. flat landscape if a=b=c=0): seems to work 
      if (is.null(c)==F){
       fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds)',sep='') 
      }
      fun_text=paste(fun_text,')}',sep='') # the end parenthesis
      fun=eval(parse(text=fun_text))
        }      
      
      # only c varies (e.g. flat landscape if a=b=0): seems to work   
      else {
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[2]*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }
}
      # only a is fixed (e.g. quadratic landscape if a=0): seems to work 
    else {
      fun_text='fun=function(X){return('
      for (i in 1:length(trees)){
        fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+X[2]*SEQ^2+X[3]*SEQ),bounds=bounds)',sep='') 
      }
      fun_text=paste(fun_text,')}',sep='') # the end parenthesis
      fun=eval(parse(text=fun_text))
      } 
  }
  
  # the full model: no parameter fixed
  else {
    fun_text='fun=function(X){return('
    for (i in 1:length(trees)){
      fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[1],dMat=DiffMat_backwards(X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ),bounds=bounds)',sep='') 
    }
    fun_text=paste(fun_text,')}',sep='') # the end parenthesis
    fun=eval(parse(text=fun_text))
  }
  return(list(fun=fun,ncoeff=ncoeff,par_fixed=list(a=a,b=b,c=c,bounds=bounds),trees=trees,traits=traits,Npts=Npts))
  
}

####################################################################################
# Optimization function
find.mle_FPK_multiple_clades_same_V_same_sig2=function(model,method='Nelder-Mead',init.optim=NULL){
    if(is.null(init.optim)==T){init.optim=c(-5,rep(0,model$ncoeff))}
    else{}
    opt=optim(par=init.optim,fn=model$fun,method=method,control=list(maxit=50000))
  par_fixed=model$par_fixed
  par=list()
  par$sigsq=2*exp(opt$par[1])
  if (model$ncoeff==3){par$a=opt$par[2] ; par$b=opt$par[3] ; par$c=opt$par[4]}
  else if (model$ncoeff==2){par$b=opt$par[2] ; par$c=opt$par[3]}
  else if (model$ncoeff==1){par$c=opt$par[2]}
  else{}
  # now retrieve x0: not needed actually, since this could be done using the uncertainty function
  return(res=list(lnL=-opt$value,aic=2*(length(init.optim)+opt$value),k=length(init.optim),par=par,par_fixed=par_fixed,convergence=opt$convergence,message=opt$message,trees=model$trees,traits=model$traits,Npts=model$Npts)) 
}

#######################################################
#### Wrapper for independent models in each clade #####
#######################################################
# BBMV
fit_BBMV_multiple_clades_different_V_different_sig2=function(trees,traits,bounds,a=NULL,b=NULL,c=NULL,Npts=50,method='Nelder-Mead',init.optim=NULL){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  lnls=ks=rep(NA,length(trees))
  fits=list()
  for (i in 1:length(trees)){
    lnl_temp=lnL_BBMV(trees[[i]],traits[[i]],bounds=bounds,a=a,b=b,c=c,Npts=Npts)
    fit_temp=find.mle_FPK(model=lnl_temp,method=method,init.optim=init.optim)
    fits[[i]]=fit_temp ; names(fits)[i]=paste('fit_clade_',i,sep='')
    lnls[i]=fit_temp$lnL ; ks[i]=fit_temp$k
  }
  return(list(lnL=sum(lnls),aic=2*(sum(ks)-sum(lnls)),k=sum(ks),fits=fits))
}

# FPK
fit_FPK_multiple_clades_different_V_different_sig2=function(trees,traits,a=NULL,b=NULL,c=NULL,Npts=50,method='Nelder-Mead',init.optim=NULL){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  # define bounds far way  
  # min_trait=min(traits[[1]]) ; max_trait=max(traits[[1]])
  # for (i in 2:length(trees)){
  #   min_trait=min(min_trait,min(traits[[i]])) ; max_trait=max(max_trait,max(traits[[i]]))
  # }
  # bounds=c(min_trait-(max_trait-min_trait)/2,max_trait+(max_trait-min_trait)/2)
  # new formulation for when there is measurement error
  span=c(min(unlist(traits)),max(unlist(traits)))
  bounds=c(span[1]-(span[2]-span[1])/2,span[2]+(span[2]-span[1])/2)  
  return(fit_BBMV_multiple_clades_different_V_different_sig2(trees=trees,traits=traits,bounds=bounds,a=a,b=b,c=c,Npts=Npts,method=method,init.optim=init.optim))
}

####################################################################################################
########### Likelihood function for FPK when V IS SHARED between trees BUT NOT SIGMA ###############
####################################################################################################
lnl_FPK_multiclades_same_V_different_sig2=function(trees,traits,a=NULL,b=NULL,c=NULL,Npts=50){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_FPK instead')}
  for (i in 1:length(trees)){
    if (sum(trees[[i]]$tip.label%in%names(traits[[i]]))<max(length(traits[[i]]),length(trees[[i]]$tip.label))){stop(paste('Tip names in tree ',i,' do not match names of corresponding trait vector'))}
  }
# define bounds far way  
#   min_trait=min(traits[[1]]) ; max_trait=max(traits[[1]])
# for (i in 2:length(trees)){
#   min_trait=min(min_trait,min(traits[[i]])) ; max_trait=max(max_trait,max(traits[[i]]))
# }
#   bounds=c(min_trait-(max_trait-min_trait)/2,max_trait+(max_trait-min_trait)/2)
# new formulation for when there is measurement error
  span=c(min(unlist(traits)),max(unlist(traits)))
  bounds=c(span[1]-(span[2]-span[1])/2,span[2]+(span[2]-span[1])/2)
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  trees_formatted=list()
  for (i in 1:length(trees)){
    trees_formatted[[i]]=FormatTree_bounds(trees[[i]],traits[[i]],V=rep(0,Npts),bounds=bounds)
  }
  ncoeff=(is.null(a)==T)+(is.null(b)==T)+(is.null(c)==T) # OK
  npar=ncoeff+length(trees)
  par_names=c() 
  for (i in 1:length(trees)){par_names=c(par_names,paste('dCoeff_tree_',i,sep=''))}
  par_names=c(par_names,'a','b','c')
  if (is.null(a)==F){
    if (is.null(b)==F){
      # all three shape parameters fixed (e.g. flat landscape if a=b=c=0): seems to work 
      if (is.null(c)==F){
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }      
      
      # only c varies (e.g. flat landscape if a=b=0): seems to work   
      else {
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[',length(trees)+1,']*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }
    }
    # only a is fixed (e.g. quadratic landscape if a=0): seems to work 
    else {
      fun_text='fun=function(X){return('
      for (i in 1:length(trees)){
        fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+X[',length(trees)+1,']*SEQ^2+X[',length(trees)+2,']*SEQ),bounds=bounds)',sep='') 
      }
      fun_text=paste(fun_text,')}',sep='') # the end parenthesis
      fun=eval(parse(text=fun_text))
    } 
  }
  
  # the full model: no parameter fixed
  else {
    fun_text='fun=function(X){return('
    for (i in 1:length(trees)){
      fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(X[',length(trees)+1,']*SEQ^4+X[',length(trees)+2,']*SEQ^2+X[',length(trees)+3,']*SEQ),bounds=bounds)',sep='') 
    }
    fun_text=paste(fun_text,')}',sep='') # the end parenthesis
    fun=eval(parse(text=fun_text))
  }
  return(list(fun=fun,ncoeff=ncoeff,npar=npar,par_names=par_names,par_fixed=list(a=a,b=b,c=c,bounds=bounds),trees=trees,traits=traits,Npts=Npts))
}

####################################################################################################
############# Likelihood function for FPK when ALL PARAMS ARE SHARED between trees #################
####################################################################################################
lnl_FPK_multiclades_same_V_same_sig2=function(trees,traits,a=NULL,b=NULL,c=NULL,Npts=50){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  for (i in 1:length(trees)){
    if (sum(trees[[i]]$tip.label%in%names(traits[[i]]))<max(length(traits[[i]]),length(trees[[i]]$tip.label))){stop(paste('Tip names in tree ',i,' do not match names of corresponding trait vector'))}
  }
  # define bounds far way  
  # min_trait=min(traits[[1]]) ; max_trait=max(traits[[1]])
  # for (i in 2:length(trees)){
  #   min_trait=min(min_trait,min(traits[[i]])) ; max_trait=max(max_trait,max(traits[[i]]))
  # }
  # bounds=c(min_trait-(max_trait-min_trait)/2,max_trait+(max_trait-min_trait)/2)
  # new formulation for when there is measurement error
  span=c(min(unlist(traits)),max(unlist(traits)))
  bounds=c(span[1]-(span[2]-span[1])/2,span[2]+(span[2]-span[1])/2)  
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  trees_formatted=list()
  for (i in 1:length(trees)){
    trees_formatted[[i]]=FormatTree_bounds(trees[[i]],traits[[i]],V=rep(0,Npts),bounds=bounds)
  }
  ncoeff=(is.null(a)==T)+(is.null(b)==T)+(is.null(c)==T) # OK, but to edit if we want some params to vary between clades
  if (is.null(a)==F){
    if (is.null(b)==F){
      # all three shape parameters fixed (e.g. flat landscape if a=b=c=0): seems to work 
      if (is.null(c)==F){
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }      
      
      # only c varies (e.g. flat landscape if a=b=0): seems to work   
      else {
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[2]*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }
    }
    # only a is fixed (e.g. quadratic landscape if a=0): seems to work 
    else {
      fun_text='fun=function(X){return('
      for (i in 1:length(trees)){
        fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+X[2]*SEQ^2+X[3]*SEQ),bounds=bounds)',sep='') 
      }
      fun_text=paste(fun_text,')}',sep='') # the end parenthesis
      fun=eval(parse(text=fun_text))
    } 
  }
  
  # the full model: no parameter fixed
  else {
    fun_text='fun=function(X){return('
    for (i in 1:length(trees)){
      fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[1],dMat=DiffMat_backwards(X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ),bounds=bounds)',sep='') 
    }
    fun_text=paste(fun_text,')}',sep='') # the end parenthesis
    fun=eval(parse(text=fun_text))
  }
  return(list(fun=fun,ncoeff=ncoeff,par_fixed=list(a=a,b=b,c=c,bounds=bounds),trees=trees,traits=traits,Npts=Npts))
  
}

####################################################################################################
#### Reformat results of multiclade fit so that ACE, charac_time and get_landscape can be used #####
####################################################################################################
reformat_multiclade_results=function(fit){
  n_clades=length(fit$trees)
  fits=list()
  if ('sigsq_tree_1' %in% names(fit$par)){ # then we have different sigmas in each clade
    for (i in 1:n_clades){
      eval(parse(text=paste('fits$fit_clade_',i,'=list()',sep='')))
      fits[[i]]$par=list()
      fits[[i]]$par$sigsq=fit$par[[i]]
      if ('a'%in%names(fit$par)){fits[[i]]$par$a=fit$par$a}
      if ('b'%in%names(fit$par)){fits[[i]]$par$b=fit$par$b}
      if ('c'%in%names(fit$par)){fits[[i]]$par$c=fit$par$c}
      fits[[i]]$par_fixed=fit$par_fixed
      fits[[i]]$tree=fit$trees[[i]]
      fits[[i]]$trait=fit$trait[[i]]
      fits[[i]]$Npts=fit$Npts
    }
  } 
  else {
    for (i in 1:n_clades){
      eval(parse(text=paste('fits$fit_clade_',i,'=list()',sep='')))
      fits[[i]]$par=list()
      fits[[i]]$par$sigsq=fit$par$sigsq
      if ('a'%in%names(fit$par)){fits[[i]]$par$a=fit$par$a}
      if ('b'%in%names(fit$par)){fits[[i]]$par$b=fit$par$b}
      if ('c'%in%names(fit$par)){fits[[i]]$par$c=fit$par$c}
      fits[[i]]$par_fixed=fit$par_fixed
      fits[[i]]$tree=fit$trees[[i]]
      fits[[i]]$trait=fit$trait[[i]]
      fits[[i]]$Npts=fit$Npts
    }  
  }
  return(fits)
}

# Set of functions to fit the Bounded Brownian Motion model + a potential with shape a.x^4+b.x^2+c.x in a MCMC framework, using a simple Metropolis-Hastings algorithm
# Written by F. Boucher October 2015 - May 2016


########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################





################################################
##########  MCMC FIT FOR THE MULTICLADE ########
##########   WITH DIFFERENT SIGMA2 BUT  ########
########## SAME POTENTIAL ACROSS CLADES ########
################################################

###############################################
######## PRIORS AND PROPOSALS FUNCTIONS #######
###############################################
# log prior on all params (n_clades x sigsqs, a, b, c) but not root position: seems to work
log_prior_nclades_plus_3_pars=function(type=NULL,shape=NULL,pars,n_clades){
  # params: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c)
  # type: either uniform of normal prior for each param, only uniform (discrete) for root position
  # the prior is on log(sigsq/2)=dCoeff, not sigsq
  # shape: params of the prior, min/max for uniform , mean/sd for normal, no shape for root position yet (uniform discrete between 1 and Npts)
  npars=n_clades+3
  if (is.null(type)){
    type=list()
    shape=list()
    for (i in 1:npars){type[[i]]='Normal'; shape[[i]]=c(0,10)}
  }
  p=rep(NA,npars)
  for (i in 1:npars){
    if (type[i]=='Normal'){p[i]=dnorm(x=pars[i],mean=shape[[i]][1],sd=shape[[i]][2])}
    if (type[i]=='Uniform'){p[i]=dunif(x=pars[i],min=shape[[i]][1],max=shape[[i]][2])}
  }
  return(sum(log(p)))
}

# proposal function to move from one value of the parameters to the next: seems to work too
proposal_nclades_plus_3_pars=function(type='Uniform',sensitivity,pars,n_clades){ 
  # same params: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c)
  npars=n_clades+3
  if (type=='Uniform'){
    par_temp=rep(NA,length(pars))
    for (i in 1:npars){
      par_temp[i]=runif(n=1,min=pars[i]-sensitivity[i],max=pars[i]+sensitivity[i])
    }
  }
  if (type=='Normal'){
    par_temp=rep(NA,length(pars))
    for (i in 1:npars){
      par_temp[i]=rnorm(n=1,mean=pars[i],sd= sensitivity[i])
    }
  }
  return(par_temp)
}

#####################################
############ MCMC SAMPLER ###########
#####################################
# MCMC sampler for a.x^4+b.x^2+c.x
MH_MCMC_FPK_multiclades=function(trees,traits,bounds,Nsteps=500000,record_every=100,plot_every=500,Npts=50,pars_init=NULL,prob_update=NULL,verbose=TRUE,plot=TRUE,save_to='MCMC_FPK_test.Rdata',save_every=10000,type_priors=NULL,shape_priors=NULL,proposal_type='Normal',proposal_sensitivity=NULL,prior.only=F,burnin.plot=0.1){
  # the oder of parameters is the same for pars_init, prob_update,type_priors,shape_priors and proposal_sensitivity. It is: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c) , with dCoeff=log(sigsq/2)
  # prior.only to sample from prior only (check that MCMC algorithm mixes well). Default to F for actual posterior exploration	
  # burnin.plot gives the proportion burnin for plots only (the whole chain is actually saved)  
  # we update parameters separately: prob_update gives the probability that each param is updated
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  n_clades=length(trees) ; npars=n_clades+3
  for (i in 1:n_clades){
    if (sum(trees[[i]]$tip.label%in%names(traits[[i]]))<max(length(traits[[i]]),length(trees[[i]]$tip.label))){stop(paste('Tip names in tree ',i,' do not match names of corresponding trait vector'))}
    ###### new piece of code added for Measurment error incorporation
    if (is.numeric(traits[[i]])){
      if ((min(traits[[i]])<bounds[1])|(max(traits[[i]])>bounds[2])){stop(paste('Some values in trait ',i,' vector exceed the bounds.'))} 
    }
    if (class(traits[[i]])=='list') {
      if ((min(unlist(traits[[i]]))<bounds[1])|(max(unlist(traits[[i]]))>bounds[2])){stop(paste('Some values in trait ',i,' vector exceed the bounds.'))}
    }
  }  
  if (is.null(prob_update)){prob_update=rep(1/npars,npars)}
  if (is.null(pars_init)){pars_init=c(rnorm(n=n_clades,mean=-8,sd=3),rnorm(n=3,mean=0,sd=2))}
  if (is.null(proposal_sensitivity)){proposal_sensitivity=rep(0.1,npars)}
  ######  end new code  
  SEQ=seq(from=-1.5,to=1.5,length.out= Npts) # the potential V is modelled as a quadratic function over [-1.5,1.5], but in real data space, this corresponds to [bounds[1],bounds[2]]
  V_init= pars_init[(n_clades+1)]*SEQ^4+pars_init[(n_clades+2)]*SEQ^2+pars_init[(n_clades+3)]*SEQ
  temp= pars_init
  chain=matrix(NA,Nsteps/record_every,(n_clades+9))
  colnames(chain)[c(1,(n_clades+2):(n_clades+9))]=c('step','a','b','c','lnprior','lnlik','quasi-lnpost','Accept','Par_updated') 
  for (clade in 1:n_clades){eval(parse(text=paste("colnames(chain)[",clade,"+1]='sigsq_clade_",clade,"'",sep='')))}
  if (prior.only==T){lnlik=1}
  else {
    NEG_LNL_func=lnl_BBMV_multiclades_same_V_different_sig2(trees=trees,traits=traits,bounds=bounds,a=NULL,b=NULL,c=NULL,Npts=50)$fun
    LNL_func=function(X){return(-NEG_LNL_func(X))}
    lnlik= LNL_func(X=temp)
  }
  lnprior= log_prior_nclades_plus_3_pars(type=type_priors,shape=shape_priors,pars=temp,n_clades=n_clades)
  lnpost=lnlik+ lnprior
  if ((is.na(lnpost))|(lnpost==(-Inf))){stop('Likelihood cannot be estimated at initial parameters. Please change them')}
  for (i in 1:Nsteps){
    par_to_update=sample(1:length(pars_init),size=1,prob=prob_update) # sample which parameter will be updated in this step
    sensitivity_temp=rep(0,length(pars_init)) # set all sensitivities to 0 so that parameters are not updated...
    sensitivity_temp[par_to_update]= proposal_sensitivity[par_to_update] #... except the one chosen
    prop= proposal_nclades_plus_3_pars(type=proposal_type,sensitivity=sensitivity_temp,pars=temp,n_clades=n_clades)
    lnprior_proposed=log_prior_nclades_plus_3_pars(type=type_priors,shape=shape_priors,pars=prop,n_clades=n_clades) 
    if (lnprior_proposed ==(-Inf)){lnpost_proposed=-Inf} # no lnl calculation when prior is null
    else {
      if (prior.only==T){lnlik_proposed=1}
      else {
        #	V_proposed= prop[(n_clades+1)]*SEQ^4+ prop[(n_clades+2)]*SEQ^2+ prop[(n_clades+3)]*SEQ # proposed potential
        lnlik_proposed= try(LNL_func(X=prop))
      }
      if (class(lnlik_proposed)=='try-error'){lnlik_proposed=NaN} # redo steps where likelihood calculation failed
      lnpost_proposed= lnlik_proposed + lnprior_proposed # un-normalized log-posterior
    }
    ALPHA=exp(lnpost_proposed-lnpost) # acceptance ratio (ratio of un-norm. posteriors)
    if (is.nan(ALPHA)){i=i-1} # re-do this step if it produces an NaN
    else {
      U=runif(1)
      if (U<ALPHA){
        temp= prop
        lnlik= lnlik_proposed
        lnpost= lnpost_proposed
        lnprior= lnprior_proposed
        accept=1
      }
      else {accept=0}
      if (i%%record_every==0){
        chain[(i/record_every),c(1,(n_clades+2):(n_clades+9))]=c(i,temp[(n_clades+1):(n_clades+3)],lnprior,lnlik, lnpost,accept, par_to_update)
        for (clade in 1:n_clades){chain[(i/record_every),(clade+1)]=2*exp(temp[clade])}
      }
      if (i%%plot_every==0){
        if (verbose==T){
          print(chain[(i/record_every),])
        }
        if (plot==T){
          par(mfrow=c(ceiling(sqrt(n_clades+6)),ceiling(sqrt(n_clades+6))))
          for (clade in 1:n_clades){plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(clade+1)],type='l',main=paste('sigsq_clade_',clade,sep=''),log='y',ylab='',xlab='')}
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+2)],type='l',main='a (x^4 term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+3)],type='l',main='b (x^2 term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+4)],type='l',main='c (x term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+5)],type='l',main='lnprior',ylab='',xlab='')
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+6)],type='l',main='lnlik',ylab='',xlab='')
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+7)],type='l',main='quasi-lnpost',ylab='',xlab='')
        }
      }
      if (i%%save_every==0){
        save(chain,file=save_to)
      }	
    }
  }
  return(chain)
}

















################################################
##########  MCMC FIT FOR THE MULTICLADE ########
##########       HIERARCHICAL MODEL     ########
########## SAME POTENTIAL ACROSS CLADES ########
################################################

###############################################
######## PRIORS AND PROPOSALS FUNCTIONS #######
###############################################
# log prior on all params (n_clades x sigsqs, a, b, c) but not root position: seems to work
log_prior_nclades_plus_3_pars_hierarchical=function(type=NULL,shape=NULL,pars,n_clades){
  # params: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c,HYPERMU,HYPERSIG)
  # HYPERMU is the mean of the prior distribution for dCoeffs
  # HYPERSIG is the sd of the prior distribution for dCoeffs, which is gaussian
  # type: either uniform of normal prior for each param, only uniform (discrete) for root position
  # the prior is on log(sigsq/2)=dCoeff, not sigsq
  # shape: params of the prior, min/max for uniform , mean/sd for normal, no shape for root position yet (uniform discrete between 1 and Npts)
  npars=n_clades+5
  if (is.null(type)){
    type=list()
    shape=list()
    for (i in 1:npars){type[[i]]='Normal'; shape[[i]]=c(0,10)}
  }
  p=rep(NA,npars)
  for (i in 1:npars){
    if (type[i]=='Normal'){p[i]=dnorm(x=pars[i],mean=shape[[i]][1],sd=shape[[i]][2])}
    if (type[i]=='Uniform'){p[i]=dunif(x=pars[i],min=shape[[i]][1],max=shape[[i]][2])}
  }
  return(sum(log(p)))
}

# proposal function to move from one value of the parameters to the next: seems to work too
proposal_nclades_plus_3_pars_hierarchical=function(type='Uniform',sensitivity,pars,n_clades){ 
  # same params: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c)
  npars=n_clades+3
  if (type=='Uniform'){
    par_temp=rep(NA,length(pars))
    for (i in 1:npars){
      par_temp[i]=runif(n=1,min=pars[i]-sensitivity[i],max=pars[i]+sensitivity[i])
    }
  }
  if (type=='Normal'){
    par_temp=rep(NA,length(pars))
    for (i in 1:npars){
      par_temp[i]=rnorm(n=1,mean=pars[i],sd= sensitivity[i])
    }
  }
  return(par_temp)
}

#####################################
############ MCMC SAMPLER ###########
#####################################
# MCMC sampler for a.x^4+b.x^2+c.x
MH_MCMC_FPK_multiclades=function(trees,traits,bounds,Nsteps=500000,record_every=100,plot_every=500,Npts=50,pars_init=NULL,prob_update=NULL,verbose=TRUE,plot=TRUE,save_to='MCMC_FPK_test.Rdata',save_every=10000,type_priors=NULL,shape_priors=NULL,proposal_type='Normal',proposal_sensitivity=NULL,prior.only=F,burnin.plot=0.1){
  # the oder of parameters is the same for pars_init, prob_update,type_priors,shape_priors and proposal_sensitivity. It is: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c) , with dCoeff=log(sigsq/2)
  # prior.only to sample from prior only (check that MCMC algorithm mixes well). Default to F for actual posterior exploration	
  # burnin.plot gives the proportion burnin for plots only (the whole chain is actually saved)  
  # we update parameters separately: prob_update gives the probability that each param is updated
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  n_clades=length(trees) ; npars=n_clades+3
  for (i in 1:n_clades){
    if (sum(trees[[i]]$tip.label%in%names(traits[[i]]))<max(length(traits[[i]]),length(trees[[i]]$tip.label))){stop(paste('Tip names in tree ',i,' do not match names of corresponding trait vector'))}
    ###### new piece of code added for Measurment error incorporation
    if (is.numeric(traits[[i]])){
      if ((min(traits[[i]])<bounds[1])|(max(traits[[i]])>bounds[2])){stop(paste('Some values in trait ',i,' vector exceed the bounds.'))} 
    }
    if (class(traits[[i]])=='list') {
      if ((min(unlist(traits[[i]]))<bounds[1])|(max(unlist(traits[[i]]))>bounds[2])){stop(paste('Some values in trait ',i,' vector exceed the bounds.'))}
    }
  }  
  if (is.null(prob_update)){prob_update=rep(1/npars,npars)}
  if (is.null(pars_init)){pars_init=c(rnorm(n=n_clades,mean=-8,sd=3),rnorm(n=3,mean=0,sd=2))}
  if (is.null(proposal_sensitivity)){proposal_sensitivity=rep(0.1,npars)}
  ######  end new code  
  SEQ=seq(from=-1.5,to=1.5,length.out= Npts) # the potential V is modelled as a quadratic function over [-1.5,1.5], but in real data space, this corresponds to [bounds[1],bounds[2]]
  V_init= pars_init[(n_clades+1)]*SEQ^4+pars_init[(n_clades+2)]*SEQ^2+pars_init[(n_clades+3)]*SEQ
  temp= pars_init
  chain=matrix(NA,Nsteps/record_every,(n_clades+9))
  colnames(chain)[c(1,(n_clades+2):(n_clades+9))]=c('step','a','b','c','lnprior','lnlik','quasi-lnpost','Accept','Par_updated') 
  for (clade in 1:n_clades){eval(parse(text=paste("colnames(chain)[",clade,"+1]='sigsq_clade_",clade,"'",sep='')))}
  if (prior.only==T){lnlik=1}
  else {
    NEG_LNL_func=lnl_BBMV_multiclades_same_V_different_sig2(trees=trees,traits=traits,bounds=bounds,a=NULL,b=NULL,c=NULL,Npts=50)$fun
    LNL_func=function(X){return(-NEG_LNL_func(X))}
    lnlik= LNL_func(X=temp)
  }
  lnprior= log_prior_nclades_plus_3_pars(type=type_priors,shape=shape_priors,pars=temp,n_clades=n_clades)
  lnpost=lnlik+ lnprior
  if ((is.na(lnpost))|(lnpost==(-Inf))){stop('Likelihood cannot be estimated at initial parameters. Please change them')}
  for (i in 1:Nsteps){
    par_to_update=sample(1:length(pars_init),size=1,prob=prob_update) # sample which parameter will be updated in this step
    sensitivity_temp=rep(0,length(pars_init)) # set all sensitivities to 0 so that parameters are not updated...
    sensitivity_temp[par_to_update]= proposal_sensitivity[par_to_update] #... except the one chosen
    prop= proposal_nclades_plus_3_pars(type=proposal_type,sensitivity=sensitivity_temp,pars=temp,n_clades=n_clades)
    lnprior_proposed=log_prior_nclades_plus_3_pars(type=type_priors,shape=shape_priors,pars=prop,n_clades=n_clades) 
    if (lnprior_proposed ==(-Inf)){lnpost_proposed=-Inf} # no lnl calculation when prior is null
    else {
      if (prior.only==T){lnlik_proposed=1}
      else {
        #	V_proposed= prop[(n_clades+1)]*SEQ^4+ prop[(n_clades+2)]*SEQ^2+ prop[(n_clades+3)]*SEQ # proposed potential
        lnlik_proposed= try(LNL_func(X=prop))
      }
      if (class(lnlik_proposed)=='try-error'){lnlik_proposed=NaN} # redo steps where likelihood calculation failed
      lnpost_proposed= lnlik_proposed + lnprior_proposed # un-normalized log-posterior
    }
    ALPHA=exp(lnpost_proposed-lnpost) # acceptance ratio (ratio of un-norm. posteriors)
    if (is.nan(ALPHA)){i=i-1} # re-do this step if it produces an NaN
    else {
      U=runif(1)
      if (U<ALPHA){
        temp= prop
        lnlik= lnlik_proposed
        lnpost= lnpost_proposed
        lnprior= lnprior_proposed
        accept=1
      }
      else {accept=0}
      if (i%%record_every==0){
        chain[(i/record_every),c(1,(n_clades+2):(n_clades+9))]=c(i,temp[(n_clades+1):(n_clades+3)],lnprior,lnlik, lnpost,accept, par_to_update)
        for (clade in 1:n_clades){chain[(i/record_every),(clade+1)]=2*exp(temp[clade])}
      }
      if (i%%plot_every==0){
        if (verbose==T){
          print(chain[(i/record_every),])
        }
        if (plot==T){
          par(mfrow=c(ceiling(sqrt(n_clades+6)),ceiling(sqrt(n_clades+6))))
          for (clade in 1:n_clades){plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(clade+1)],type='l',main=paste('sigsq_clade_',clade,sep=''),log='y',ylab='',xlab='')}
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+2)],type='l',main='a (x^4 term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+3)],type='l',main='b (x^2 term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+4)],type='l',main='c (x term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+5)],type='l',main='lnprior',ylab='',xlab='')
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+6)],type='l',main='lnlik',ylab='',xlab='')
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+7)],type='l',main='quasi-lnpost',ylab='',xlab='')
        }
      }
      if (i%%save_every==0){
        save(chain,file=save_to)
      }	
    }
  }
  return(chain)
}
