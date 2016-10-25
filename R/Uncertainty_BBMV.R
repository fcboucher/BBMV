Uncertainty_BBMV=function(model,tree,trait,Npts=50,effort_uncertainty= 100){
	npar=length(model[[1]]) # 3: flat ; 4: linear ; 5: quadratic ; 6: x4
	# retrieve coefficients for the potential
	if (npar==3){model$par$a=model$par$b=model$par$c=0}
	if (npar==4){model$par$a=model$par$b=0 ; pars=model[[1]][-which(names(model[[1]])%in%c('bounds','sigsq','root_value'))] ; model$par$c=pars[[1]]}
	if (npar==5){model$par$a=0 ; pars=model[[1]][-which(names(model[[1]])%in%c('bounds','sigsq','root_value'))] ; model$par$b=pars[[1]] ; model$par$c=pars[[2]]}
	if (npar==6){model$par$a=model[[1]]$a ; model$par$b=model[[1]]$b ; model$par$c=model[[1]]$c}
	
	par(mfrow=c(2,4))
	
	# for each of the parameters (from 4 to 7) fix all other params to their MLE and then do n= effort_uncertainty estimations of the log-lik, varying the value of this param
	Npts_tot=round(Npts*(model$par$bounds[2]-model$par$bounds[1])/(max(trait)-min(trait)))
	tree_formatted= FormatTree_bounds(tree,trait,V=rep(0, Npts_tot),bounds=model$par$bounds)
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts_tot) # very sensitive
	
	# sig2
	# define V
	V=model$par$a*SEQ^4+model$par$b*SEQ^2+model$par$c*SEQ
	dMat=DiffMat_backwards(V)
	ll_sig2=as.data.frame(matrix(NA, effort_uncertainty,2))
    colnames(ll_sig2)=c('sig2','loglik')
    ll_sig2$sig2 =seq(from=(model$par$sigsq/10),to=(model$par$sigsq*10), length.out=effort_uncertainty)
    ll_sig2$loglik=sapply(ll_sig2$sig2,FUN=function(x){LogLik_bounds(tree_formatted,log(x/2),dMat,bounds=model$par$bounds)})
    if (sum(is.nan(ll_sig2$loglik))>0) {ll_sig2=ll_sig2[-which(is.nan(ll_sig2$loglik)),]}
    ll_cond3=exp(ll_sig2[,2])/sum(exp(ll_sig2[,2]))
    CDF_right=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    CDF_left=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    CI95_sig2=c(ll_sig2[min(which(CDF_right>0.025)),1],ll_sig2[max(which(CDF_left>0.025)),1])
	plot(ll_sig2, main='log-likelihood vs sigsq',type='l') ; abline(v=model$par$sigsq,col=2)
	
	# a
	if ((model$par$a)!=0){
		ll_a=as.data.frame(matrix(NA, effort_uncertainty,2))
    	colnames(ll_a)=c('a','loglik')
    	ll_a$a =sort(seq(from=-10,to=10, length.out=effort_uncertainty))
		loglik_a_temp=function(x){
			V=x*SEQ^4+model$par$b*SEQ^2+model$par$c*SEQ
			dMat=DiffMat_backwards(V)
			return(LogLik_bounds(tree_formatted,log(model$par$sigsq/2),dMat,bounds=model$par$bounds))
		}
		ll_a$loglik=sapply(ll_a$a,FUN= loglik_a_temp)
		if (sum(is.nan(ll_a$loglik))>0) {ll_a=ll_a[-which(is.nan(ll_a$loglik)),]}
		ll_cond3=exp(ll_a[,2])/sum(exp(ll_a[,2]))
    	CDF_right=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    	CDF_left=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    	CI95_a=c(ll_a[min(which(CDF_right>0.025)),1], ll_a[max(which(CDF_left>0.025)),1])
    	plot(ll_a, main='log-likelihood vs a (x^4 term)',type='l') ; abline(v=model$par$a,col=2)
	}
	else {CI95_a=NA}
	
	# b
	if ((model$par$b)!=0){
		ll_b=as.data.frame(matrix(NA, effort_uncertainty,2))
    	colnames(ll_b)=c('b','loglik')
    	ll_b$b =sort(seq(from=-10,to=10, length.out=effort_uncertainty))
		loglik_b_temp=function(x){
			V=model$par$a*SEQ^4+x*SEQ^2+model$par$c*SEQ
			dMat=DiffMat_backwards(V)
			return(LogLik_bounds(tree_formatted,log(model$par$sigsq/2),dMat,bounds=model$par$bounds))
		}
		ll_b$loglik=sapply(ll_b$b,FUN= loglik_b_temp)
		if (sum(is.nan(ll_b$loglik))>0) {ll_b=ll_b[-which(is.nan(ll_b$loglik)),]}
		ll_cond3=exp(ll_b[,2])/sum(exp(ll_b[,2]))
    	CDF_right=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    	CDF_left=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    	CI95_b=c(ll_b[min(which(CDF_right>0.025)),1], ll_b[max(which(CDF_left>0.025)),1])
    	plot(ll_b, main='log-likelihood vs b (x^2 term)',type='l') ; abline(v=model$par$b,col=2)
	}
	else {CI95_b=NA}

	# c
	if ((model$par$c)!=0){
		ll_c=as.data.frame(matrix(NA, effort_uncertainty,2))
    	colnames(ll_c)=c('c','loglik')
    	ll_c$c =sort(seq(from=-10,to=10, length.out=effort_uncertainty))
		loglik_c_temp=function(x){
			V=model$par$a*SEQ^4+model$par$b*SEQ^2+x*SEQ
			dMat=DiffMat_backwards(V)
			return(LogLik_bounds(tree_formatted,log(model$par$sigsq/2),dMat,bounds=model$par$bounds))
		}
		ll_c$loglik=sapply(ll_c$c,FUN= loglik_c_temp)
		if (sum(is.nan(ll_c$loglik))>0) {ll_c=ll_c[-which(is.nan(ll_c$loglik)),]}
		ll_cond3=exp(ll_c[,2])/sum(exp(ll_c[,2]))
    	CDF_right=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    	CDF_left=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    	CI95_c=c(ll_c[min(which(CDF_right>0.025)),1], ll_c[max(which(CDF_left>0.025)),1])
    	plot(ll_c, main='log-likelihood vs c (x term)',type='l') ; abline(v=model$par$c,col=2)
	}
	else {CI95_c=NA}
	
	# x0
	V=model$par$a*SEQ^4+model$par$b*SEQ^2+model$par$c*SEQ
	dMat=DiffMat_backwards(V)
	ll_root=as.data.frame(matrix(NA, Npts_tot,2))
    colnames(ll_root)=c('root','density')
    ll_root$root =seq(from=model$par$bounds[1],to=model$par$bounds[2], length.out= Npts_tot)
    tree_formatted2= tree_formatted
	pMat=prep_mat_exp(dCoeff=log(model$par$sigsq/2),dMat,bounds=model$par$bounds) # edited
	logFactor=0
	for (i in 1:dim(tree_formatted2$tab)[1]){
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat)
	norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
	logFactor=logFactor+log(norm)
}
    ll_root$density =tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]
    CDF_root_right=sapply(1:length(ll_root$density),FUN=function(x){sum(ll_root$density[c(1:x)])})
    CDF_root_left=sapply(1:length(ll_root$density),FUN=function(x){sum(ll_root$density[c(x:length(ll_root$density))])})
    cells_root=c(min(which(CDF_root_right>0.025)),max(which(CDF_root_left>0.025)))
    if (model$par$root_value==model$par$bounds[1]){cells_root=c(1,max(which(CDF_root_left>0.05)))}
    if (model$par$root_value== model$par$bounds[2]){cells_root=c(min(which(CDF_root_right>0.05)),length(ll_root$density))}
	CI95_root=model$par$bounds[1]+(model$par$bounds[2]-model$par$bounds[1])/(Npts_tot-1)*(cells_root-1)
	plot(ll_root, main='density vs root trait value',type='l',ylim=c(0,1)) ; abline(v=model$par$root_value,col=2)
	
		
	# for the bounds we have to re-format the tree at each estimation
	# Bmin
	ll_bmin=as.data.frame(matrix(NA, Npts,2))
    colnames(ll_bmin)=c('bmin','loglik')
    ll_bmin$bmin =seq(from=min(trait)-(max(trait)-min(trait)),to=min(trait), length.out= Npts)
	loglik_bmin_temp=function(x){
			SEQ=seq(from=-1.5,to=1.5,length.out=(2*Npts+1-x)) # very sensitive
			V=model$par$a*SEQ^4+model$par$b*SEQ^2+model$par$c*SEQ
			tree_formatted_temp= FormatTree_bounds(tree,trait,V,bounds=c(ll_bmin$bmin[x],model$par$bounds[2]))
			dMat=DiffMat_backwards(V)
			return(LogLik_bounds(tree_formatted_temp,dCoeff=log(model$par$sigsq/2),dMat,bounds=c(ll_bmin$bmin[x],model$par$bounds[2])))
		}
		ll_bmin$loglik=sapply(c(1:Npts),FUN= loglik_bmin_temp)
		if (sum(is.nan(ll_bmin$loglik))>0) {ll_bmin=ll_bmin[-which(is.nan(ll_bmin$loglik)),]}
		ll_cond3=exp(ll_bmin[,2])/sum(exp(ll_bmin[,2]))
    	CDF=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    	CI95_bmin=c(ll_bmin[max(which(CDF>0.95)),1], ll_bmin[dim(ll_bmin)[1],1])
    	plot(ll_bmin, main='log-likelihood vs Bmin',type='l') ; abline(v=model$par$bounds[1],col=2)
	
	# Bmax
	ll_bmax=as.data.frame(matrix(NA, Npts,2))
    colnames(ll_bmax)=c('bmax','loglik')
    ll_bmax$bmax =seq(from=max(trait),to=max(trait)+(max(trait)-min(trait)), length.out= Npts)
	loglik_bmax_temp=function(x){
			SEQ=seq(from=-1.5,to=1.5,length.out=(Npts-1+x)) # very sensitive
			V=model$par$a*SEQ^4+model$par$b*SEQ^2+model$par$c*SEQ
			tree_formatted_temp= FormatTree_bounds(tree,trait,V,bounds=c(model$par$bounds[1],ll_bmax$bmax[x]))
			dMat=DiffMat_backwards(V)
			return(LogLik_bounds(tree_formatted_temp,dCoeff=log(model$par$sigsq/2),dMat,bounds=c(model$par$bounds[1],ll_bmax$bmax[x])))
		}
		ll_bmax$loglik=sapply(c(1:Npts),FUN= loglik_bmax_temp)
		if (sum(is.nan(ll_bmax$loglik))>0) {ll_bmax=ll_bmax[-which(is.nan(ll_bmax$loglik)),]}
		ll_cond3=exp(ll_bmax[,2])/sum(exp(ll_bmax[,2]))
    	CDF=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    	CI95_bmax=c(ll_bmax[1,1], ll_bmax[min(which(CDF>0.95)),1])
    	plot(ll_bmax, main='log-likelihood vs Bmax',type='l') ; abline(v=model$par$bounds[2],col=2)
	
	return(list(CI95_sig2= CI95_sig2, CI95_a= CI95_a, CI95_b= CI95_b, CI95_c= CI95_c, CI95_root= CI95_root, CI95_bmin = CI95_bmin, CI95_bmax = CI95_bmax ))
}