# Set of functions to fit the Bounded Brownian Motion model + a potential with shape a.x^4+b.x^2+c.x in a MCMC framework, using a simple Metropolis-Hastings algorithm
# Written by F. Boucher October 2015 - May 2016

###############################################################################
########### LOG-LIKELIHOOD FUNCTION INCL. ROOT & BOUNDS AS PARAMS #############
###############################################################################
# calculate log-likelihood over the whole tree with root value as a parameter, to be maximized
LogLik_bounds_est_root=function(tree,trait ,dCoeff,V, x0_pos,bounds){
Npts_tot=length(V)		
tree_formatted=FormatTree_bounds(tree,trait,V,bounds)		
dMat=DiffMat_backwards(V)
pMat=prep_mat_exp(dCoeff,dMat,bounds) # edited
tree_formatted2= tree_formatted
logFactor=0
for (i in 1:dim(tree_formatted2$tab)[1]){
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat=pMat) # edited with prep_mat for faster calculations
	norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
	logFactor=logFactor+log(norm)
}
return(log(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]][x0_pos])+logFactor) # extraction of the likelihood for a given root value (in the ML optimization version, we take max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))
}

###############################################
######## PRIORS AND PROPOSALS FUNCTIONS #######
###############################################
# log prior on all 5 params including root position: seems to work
log_prior_5pars_root_bounds=function(type=c(rep('Normal',4),'Uniform'),shape=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),pars){
	# pars: the actual parameters for which to calculate the prior (dCoeff,a,b,c,x0).
	# type: either uniform of normal prior for each param, only uniform (discrete) for root position
	# the prior is on log(sigsq/2)=dCoeff, not sigsq
	# shape: params of the prior, min/max for uniform , mean/sd for normal, no shape for root position yet (uniform discrete between 1 and Npts)
	p=list() 
	for (i in 1:4){
		if (type[i]=='Normal'){p[[i]]=dnorm(x=pars[i],mean=shape[[i]][1],sd=shape[[i]][2])}
		if (type[i]=='Uniform'){p[[i]]=dunif(x=pars[i],min=shape[[i]][1],max=shape[[i]][2])}
	}
	# root should also be treated as a continuous value???? change the call in log-lik function then to retrieve the exact point. Or we shift it when we update bmin... better probably
	if (type[5]=='Uniform'){
		if (pars[5]%in%c(1:Npts)){
			p[[5]]=1/Npts
		}
		else {p[[5]]=0}
	}
	
 	else { stop('Only uniform prior on root position is supported.')}
	return(log(p[[1]]*p[[2]]*p[[3]]*p[[4]]*p[[5]]))
}

# proposal function: to move from one value of the parameters to the next SHOULD BE OK
proposal_5pars_root_bounds=function(type='Uniform',sensitivity,pars){ 
	# same 5 params: (dCoeff,a,b,c,x0)
	# for the root we draw a point on the grid nearby from root-sensitivity, to root+sensitivity, excluding root (i.e. we HAVE to move). Probably better to move one step at the time, i.e. sensitivity[5]=1. Same for the bounds
	if (type=='Uniform'){
		par_temp=c(runif(n=1,min=pars[1]-sensitivity[1],max=pars[1]+sensitivity[1]),runif(n=1,min=pars[2]-sensitivity[2],max=pars[2]+sensitivity[2]),runif(n=1,min=pars[3]-sensitivity[3],max=pars[3]+sensitivity[3]),runif(n=1,min=pars[4]-sensitivity[4],max=pars[4]+sensitivity[4]))
	}
	if (type=='Normal'){
	par_temp=c(rnorm(n=1,mean=pars[1],sd= sensitivity[1]),rnorm(n=1,mean=pars[2],sd= sensitivity[2]),rnorm(n=1,mean=pars[3],sd= sensitivity[3]),rnorm(n=1,mean=pars[4],sd= sensitivity[4]))
	}
	# update the root
	if (sensitivity[5]==0){root=pars[5]}
	else{
		root=pars[5]+sample(size=1,x=c(seq(from=-sensitivity[5],to=-1,by=1),seq(from=1,to=sensitivity[5],by=1))) }
	return(c(par_temp,root))
}

#####################################
############ MCMC SAMPLER ###########
#####################################
# MCMC sampler for a.x^4+b.x^2+c.x
# to do: save chain by chunks (every 1000 steps?) and write a new matrix ? as in .store.bayou??? For the moment we save the whole table with most rows empty, waiting to be filled
MH_MCMC_FPK=function(tree,trait,bounds,Nsteps=500000,record_every=100,plot_every=500,Npts=50,pars_init=c(0,0,0,0,25),prob_update=c(0.2,0.2,0.2,0.2,0.2),verbose=TRUE,plot=TRUE,save_to='MCMC_FPK_test.Rdata',save_every=10000,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F,burnin.plot=0.1){
# prior.only to sample from prior only (check that MCMC algorithm mixes well). Default to F for actual posterior exploration	
# burnin.plot gives the proportion burnin for plots only (the whole chain is actually saved)  
# we update parameters separately: prob_update gives the probability that each param is updated
  if (is.numeric(trait)){
    if ((min(trait)<bounds[1])|(max(trait)>bounds[2])){stop('Some values in the trait vector exceed the bounds.')} 
  }
  if (class(trait)=='list') {
    if ((min(unlist(trait))<bounds[1])|(max(unlist(trait))>bounds[2])){stop('Some values in the trait data exceed the bounds.')}
  }
SEQ=seq(from=-1.5,to=1.5,length.out= Npts) # the potential V is modelled as a quadratic function over [-1.5,1.5], but in real data space, this corresponds to [bounds[1],bounds[2]]
V_init= pars_init[2]*SEQ^4+pars_init[3]*SEQ^2+pars_init[4]*SEQ
temp= pars_init
chain=matrix(NA,Nsteps/record_every,11)
colnames(chain)=c('step','sigsq','a','b','c','root','lnprior','lnlik','quasi-lnpost','Accept','Par_updated')  
if (prior.only==T){lnlik=1}
else {
lnlik= LogLik_bounds_est_root(tree, trait,dCoeff=temp[1],x0_pos=temp[5],V= V_init,bounds=bounds)
}
lnprior= log_prior_5pars_root_bounds(type=type_priors,shape=shape_priors,pars=temp,Npts=Npts,bounds=bounds,trait=trait)
lnpost=lnlik+ lnprior
if ((is.na(lnpost))|(lnpost==(-Inf))){stop('Likelihood cannot be estimated at initial parameters. Please change them')}
for (i in 1:Nsteps){
	par_to_update=sample(1:length(pars_init),size=1,prob=prob_update) # sample which parameter will be updated in this step
	sensitivity_temp=rep(0,length(pars_init)) # set all sensitivities to 0 so that parameters are not updated...
	sensitivity_temp[par_to_update]= proposal_sensitivity[par_to_update] #... except the one chosen
	prop= proposal_5pars_root_bounds(type='Uniform',sensitivity= sensitivity_temp,pars=temp,trait)
	lnprior_proposed= log_prior_5pars_root_bounds(type=type_priors,shape=shape_priors,pars=prop,Npts=Npts,bounds=bounds,trait)
	if (lnprior_proposed ==(-Inf)){lnpost_proposed=-Inf} # no lnl calculation when prior is null
	else {
	if (prior.only==T){lnlik_proposed=1}
	else {
	#SEQ=seq(from=-1.5,to=1.5,length.out= Npts)		
	V_proposed= prop[2]*SEQ^4+ prop[3]*SEQ^2+ prop[4]*SEQ # proposed potential
	lnlik_proposed= try(LogLik_bounds_est_root(tree,trait, dCoeff=prop[1],x0_pos=prop[5], V=V_proposed,bounds=bounds)) # test with try
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
	chain[(i/record_every),]=c(i,2*exp(temp[1]),temp[2:4],bounds[1]+(temp[5]-1)*(bounds[2]-bounds[1])/(Npts-1),lnprior,lnlik, lnpost,accept, par_to_update)
	}
	if (i%%plot_every==0){
	if (verbose==T){
		print(chain[(i/record_every),])
	}
	if (plot==T){
par(mfrow=c(3,3))
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),2],type='l',main='sigsq',log='y',ylab='',xlab='')
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),6],type='l',main='root',ylab='',xlab='')
if (is.numeric(trait)){
  abline(h=min(trait)+(max(trait)-min(trait))/2,col=2)
}
if (class(trait)=='list') {
  abline(h=min(unlist(trait))+(max(unlist(trait))-min(unlist(trait)))/2,col=2)
}
plot(1,1,ylab='',xlab='')
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),3],type='l',main='a (x^4 term)',ylab='',xlab='')
abline(h=0,col=2)
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),4],type='l',main='b (x^2 term)',ylab='',xlab='')
abline(h=0,col=2)
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),5],type='l',main='c (x term)',ylab='',xlab='')
abline(h=0,col=2)
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),7],type='l',main='lnprior',ylab='',xlab='')
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),8],type='l',main='lnlik',ylab='',xlab='')
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),9],type='l',main='quasi-lnpost',ylab='',xlab='')
	}
	}
	if (i%%save_every==0){
		save(chain,file=save_to)
		}	
	}
}
return(chain)
}

#####################################
########## LANDSCAPE MCMC ###########
#####################################
get.landscape.FPK.MCMC=function(chain,bounds,Npts=100,burnin=0.1,probs.CI=c(0.05,0.95),COLOR_MEDIAN='red',COLOR_FILL='red',transparency=0.3,main='Macroevolutionary landscapes MCMC',ylab='N.exp(-V)',xlab='Trait',xlim=NULL,ylim=NULL){
  chain2=chain[-c(1:floor(dim(chain)[1]*burnin)),]  # remove burnin
  all_V=as.data.frame(matrix(NA,dim(chain2)[1],Npts))
  step=(bounds[2]-bounds[1])/(Npts-1)
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  for (i in 1:dim(chain2)[1]){
    temp=chain2[i,'a']*SEQ^4+chain2[i,'b']*SEQ^2+chain2[i,'c']*SEQ #potential
    all_V[i,]=exp(-temp)/sum(exp(-temp)*step)
  }
  if (is.null(xlim)){xlim=bounds}
  if (is.null(ylim)){ylim=c(0,max(all_V))}
  par(mfrow=c(1,1))
  plot(apply(all_V,2,function(x){quantile(x,probs=c(0.5))})~seq(from=bounds[1],to=bounds[2],length.out=Npts),col=COLOR_MEDIAN,lwd=3,type='l',main=main,ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim)
  polygon(x=c(seq(from=bounds[1],to=bounds[2],length.out=100),seq(from=bounds[2],to=bounds[1],length.out=100)),y=c(apply(all_V,2,function(x){quantile(x,probs=probs.CI[1])}),rev(apply(all_V,2,function(x){quantile(x,probs=probs.CI[2])}))),col=adjustcolor(col=COLOR_FILL,alpha.f=transparency))
}