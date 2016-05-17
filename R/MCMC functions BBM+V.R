# Set of functions to fit the Bounded Brownian Motion model + a potential with shape a.x^4+b.x^2+c.x in a MCMC framework, using a simple Metropolis-Hastings algorithm
# Written by F. Boucher October 2015 - May 2016

###############################################################################
########### LOG-LIKELIHOOD FUNCTION INCL. ROOT & BOUNDS AS PARAMS #############
###############################################################################
# calculate log-likelihood over the whole tree with root value as a parameter, to be maximized
LogLik_bounds_est_root=function(tree,trait ,dCoeff,V, x0_pos,bounds){
if ((bounds[1]>min(trait))|(bounds[2]<max(trait))) {return(-Inf)} # bounds have to be outside the trait interval
	else {
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
}

###############################################
######## PRIORS AND PROPOSALS FUNCTIONS #######
###############################################
# log prior on all 7 params including root position and bounds SHOULD BE OK
log_prior_7pars_root_bounds=function(type=c(rep('Normal',4),rep('Uniform',3)),shape=list(c(0,10),c(0,10),c(0,10),c(0,10),NA,10,10),pars,Npts_int,trait){
	# pars: the actual parameters for which to calculate the prior (dCoeff,a,b,c,x0,bmin,bmax). Although bmin and bmax only move between points of the grid, they are treated as continuous for practical purposes (all likelihood functions treat them as such...) 
	# type: either uniform of normal prior for each param, only uniform (discrete) for root position, and only uniform continuous for bounds
	# Npts_int is Npts in the interval, but if bounds lie further than min/max of the trait, then we have a real Npts that takes these extra points into account
	# the prior is on log(sigsq/2)=dCoeff, not sigsq
	# shape: params of the prior, min/max for uniform , mean/sd for normal, no shape for root position yet (uniform discrete between 1 and Npts)
	# shape[6,7]: the maximum number of steps explored away from the bounds
	step=(max(trait)-min(trait))/(Npts_int-1)
	Npts=floor((pars[7]-pars[6])/step+1) # to prevent rounding errors
	Npts_max= Npts_int+shape[[6]][1]+shape[[7]][1]
	p=list() 
	for (i in 1:4){
		if (type[i]=='Normal'){p[[i]]=dnorm(x=pars[i],mean=shape[[i]][1],sd=shape[[i]][2])}
		if (type[i]=='Uniform'){p[[i]]=dunif(x=pars[i],min=shape[[i]][1],max=shape[[i]][2])}
	}
	# root should also be treated as a continuous value???? change the call in log-lik function then to retrieve the exact point. Or we shift it when we update bmin... better probably
	if (type[5]=='Uniform'){
		if (pars[5]%in%c(1:Npts)){
			p[[5]]=1/Npts_max
		}
		else {p[[5]]=0}
	}
	# change that and make bounds only move on the grid I think
	if (type[6]=='Uniform'){p[[6]]=dunif(x=pars[6],min=min(trait)-shape[[6]][1]*step,max=min(trait))
		}
	if (type[7]=='Uniform'){p[[7]]=dunif(x=pars[7],min=max(trait),max=max(trait)+shape[[7]][1]*step)
		}	
	
# # 	else { stop('Only uniform prior on root position supported so far. More fancy things coming soon.')}
	return(log(p[[1]]*p[[2]]*p[[3]]*p[[4]]*p[[5]]*p[[6]]*p[[7]]))
}

# proposal function: to move from one value of the parameters to the next SHOULD BE OK
proposal_7pars_root_bounds=function(type='Uniform',sensitivity,pars,trait,Npts_int){ 
	# same 5 params: (dCoeff,a,b,c,x0,bmin,bmax)
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
	# update the bounds
	step=(max(trait)-min(trait))/(Npts_int-1)
	if (sensitivity[6]==0){bmin=pars[6]}
	else{
		move_bmin=sample(size=1,x=c(seq(from=-sensitivity[6],to=-1,by=1),seq(from=1,to=sensitivity[6],by=1)))
		bmin=pars[6]+ move_bmin*step 
		root=root-move_bmin # change the index of the root if we moved the bounds (the index on the grid is calculated from the minimum bound, hence we only need to change the root when we update bmin)
		}
	if (sensitivity[7]==0){bmax=pars[7]}
	else{
		move_bmax=sample(size=1,x=c(seq(from=-sensitivity[7],to=-1,by=1),seq(from=1,to=sensitivity[7],by=1)))
		bmax=pars[7]+ move_bmax*step }

	return(c(par_temp,root,bmin,bmax))
}

#####################################
############ MCMC SAMPLER ###########
#####################################
# MCMC sampler for a.x^4+b.x^2+c.x
# to do: save chain by chunks (every 1000 steps?) and write a new matrix ? as in .store.bayou??? For the moment we save the whole table with most rows empty, waiting to be filled
MH_MCMC_V_ax4bx2cx_root_bounds=function(tree,trait,Nsteps=500000,record_every=100,plot_every=500,Npts_int=50,pars_init=c(0,0,0,0,25,-10,10),prob_update=c(0.25,0.2,0.2,0.2,0.05,0.05,0.05),verbose=TRUE,plot=TRUE,save_to='~/Desktop/MCMC5_test_ax4bx2cx_root_bounds_est.R',save_every=10000,type_priors=c(rep('Normal',4),rep('Uniform',3)),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA,10,10),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1,1,1),prior.only=F){
# prior.only to sample from prior only (check MCMC algo mixes well). Default to F for actual posterior exploration	
# we update parameters separately: prob_update gives the probability that each param is updated
step=(max(trait)-min(trait))/(Npts_int-1)
bounds_init= pars_init[c(6,7)] # change that: this is a parameter now
SEQ=seq(from=-1.5,to=1.5,length.out= Npts_int) # the potential V is modelled as a quadratic function over [-1.5,1.5], but in real data space, this corresponds to [min(trait),max(trait)]
V_init= pars_init[2]*SEQ^4+pars_init[3]*SEQ^2+pars_init[4]*SEQ
#dMat_init=DiffMat_backwards(V_init)
#tree_formatted= FormatTree_bounds(tree, trait, V_init, bounds_init) # needs not be updated (only depends on length(V))
temp= pars_init
chain=matrix(NA,Nsteps/record_every,13)
colnames(chain)=c('step','sigsq','a','b','c','root','bmin','bmax','lnprior','lnlik','quasi-lnpost','Accept','Par_updated')  
if (prior.only==T){lnlik=1}
else {
lnlik= LogLik_bounds_est_root(tree, trait,dCoeff=temp[1],x0_pos=temp[5],V= V_init,bounds=temp[c(6,7)])
}
lnprior= log_prior_7pars_root_bounds(type=type_priors,shape=shape_priors,pars=temp,Npts_int,trait=trait)
lnpost=lnlik+ lnprior
if ((is.na(lnpost))|(lnpost==(-Inf))){stop('Likelihood cannot be estimated at initial parameters. Please change them')}
for (i in 1:Nsteps){
	par_to_update=sample(1:length(pars_init),size=1,prob=prob_update) # sample which parameter will be updated in this step
	sensitivity_temp=rep(0,length(pars_init)) # set all sensitivities to 0 so that parameters are not updated...
	sensitivity_temp[par_to_update]= proposal_sensitivity[par_to_update] #... except the one chosen
	prop= proposal_7pars_root_bounds(type='Uniform',sensitivity= sensitivity_temp,pars=temp,trait,Npts_int)
	lnprior_proposed= log_prior_7pars_root_bounds(type=type_priors,shape=shape_priors,pars=prop, Npts_int,trait)
	if (lnprior_proposed ==(-Inf)){lnpost_proposed=-Inf} # no lnl calculation when prior is null
	else {
	if (prior.only==T){lnlik_proposed=1}
	else {
	Npts=floor((prop[7]-prop[6])/step+1)
	SEQ=seq(from=-1.5,to=1.5,length.out= Npts)		
	V_proposed= prop[2]*SEQ^4+ prop[3]*SEQ^2+ prop[4]*SEQ # proposed potential
	lnlik_proposed= try(LogLik_bounds_est_root(tree,trait, dCoeff=prop[1],x0_pos=prop[5], V=V_proposed,bounds=prop[c(6,7)])) # test with try
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
	chain[(i/record_every),]=c(i,2*exp(temp[1]),temp[2:4],temp[6]+(temp[5]-1)*step,temp[6],temp[7],lnprior,lnlik, lnpost,accept, par_to_update)
	}
	if (i%%plot_every==0){
	if (verbose==T){
		print(chain[(i/record_every),])
	}
	if (plot==T){
par(mfrow=c(2,5))
plot(chain[1:(i/record_every),2],type='l',main='sigsq',log='y',ylab=NULL)
plot(chain[1:(i/record_every),3],type='l',main='a (x^4 term)',ylab=NULL)
plot(chain[1:(i/record_every),4],type='l',main='b (x^2 term)',ylab=NULL)
plot(chain[1:(i/record_every),5],type='l',main='c (x term)',ylab=NULL)
plot(chain[1:(i/record_every),6],type='l',main='root',ylab=NULL)
abline(h=min(trait)+(max(trait)-min(trait))/2,col=2)
plot(chain[1:(i/record_every),7],type='l',main='bmin',ylab=NULL)
abline(h=min(trait),col=2)
plot(chain[1:(i/record_every),8],type='l',main='bmax',ylab=NULL)
abline(h=max(trait),col=2)
plot(chain[1:(i/record_every),9],type='l',main='lnprior',ylab=NULL)
plot(chain[1:(i/record_every),10],type='l',main='lnlik',ylab=NULL)
plot(chain[1:(i/record_every),11],type='l',main='quasi-lnpost',ylab=NULL)
	}
	}
	if (i%%save_every==0){
		save(chain,file=save_to)
		}	
	}
}
return(chain)
}
