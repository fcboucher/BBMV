Uncertainty_BBMV=function(fit,tree,trait,Npts=50,effort_uncertainty= 100,scope_a=c(-10,10),scope_b=c(-10,10),scope_c=c(-10,10)){
  # retrieve coefficients for the potential
  if ('a'%in%names(fit$par)){a=fit$par$a}
  else {a=fit$par_fixed$a}
  if ('b'%in%names(fit$par)){b=fit$par$b}
  else {b=fit$par_fixed$b}
  if ('c'%in%names(fit$par)){c=fit$par$c}
  else {c=fit$par_fixed$c}
  bounds=fit$par_fixed$bounds
  
	par(mfrow=c(2,3))
	# for each of the parameters fix all other params to their MLE and then do n= effort_uncertainty estimations of the log-lik, varying the value of this param
	tree_formatted= FormatTree_bounds(tree,trait,V=rep(0, Npts),bounds=bounds)
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts) # very sensitive
	
	# sig2
	# define V
	V=a*SEQ^4+b*SEQ^2+c*SEQ
	dMat=DiffMat_backwards(V)
	ll_sig2=as.data.frame(matrix(NA, effort_uncertainty,2))
    colnames(ll_sig2)=c('sig2','loglik')
    ll_sig2$sig2 =seq(from=(fit$par$sigsq/10),to=(fit$par$sigsq*10), length.out=effort_uncertainty)
    ll_sig2$loglik=sapply(ll_sig2$sig2,FUN=function(x){LogLik_bounds(tree_formatted,log(x/2),dMat,bounds=bounds)})
    if (sum(is.nan(ll_sig2$loglik))>0) {ll_sig2=ll_sig2[-which(is.nan(ll_sig2$loglik)),]}
    ll_cond3=exp(ll_sig2[,2])/sum(exp(ll_sig2[,2]))
    CDF_right=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    CDF_left=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    CI95_sig2=c(ll_sig2[min(which(CDF_right>0.025)),1],ll_sig2[max(which(CDF_left>0.025)),1])
	plot(round(ll_sig2,digits=4), main='log-likelihood vs sigsq',type='l') ; abline(v=fit$par$sigsq,col=2)
	
	# a
	if ('a'%in%names(fit$par)){
		ll_a=as.data.frame(matrix(NA, effort_uncertainty,2))
    	colnames(ll_a)=c('a','loglik')
    	ll_a$a =sort(seq(from=scope_a[1],to=scope_a[2], length.out=effort_uncertainty)) # problem with [-10,10] interval
		loglik_a_temp=function(x){
			V=x*SEQ^4+b*SEQ^2+c*SEQ
			dMat=DiffMat_backwards(V)
			return(LogLik_bounds(tree_formatted,log(fit$par$sigsq/2),dMat,bounds=bounds))
		}
		ll_a$loglik=sapply(ll_a$a,FUN= loglik_a_temp)
		if (sum(is.nan(ll_a$loglik))>0) {ll_a=ll_a[-which(is.nan(ll_a$loglik)),]}
		ll_cond3=exp(ll_a[,2])/sum(exp(ll_a[,2]))
    	CDF_right=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    	CDF_left=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    	CI95_a=c(ll_a[min(which(CDF_right>0.025)),1], ll_a[max(which(CDF_left>0.025)),1])
    	plot(round(ll_a,digits=4), main='log-likelihood vs a (x^4 term)',type='l') ; abline(v=fit$par$a,col=2)
	}
	else {CI95_a=NA}
	
	# b
	if ('b'%in%names(fit$par)){
		ll_b=as.data.frame(matrix(NA, effort_uncertainty,2))
    	colnames(ll_b)=c('b','loglik')
    	ll_b$b =sort(seq(from=scope_b[1],to=scope_b[2], length.out=effort_uncertainty))
		loglik_b_temp=function(x){
			V=a*SEQ^4+x*SEQ^2+c*SEQ
			dMat=DiffMat_backwards(V)
			return(LogLik_bounds(tree_formatted,log(fit$par$sigsq/2),dMat,bounds=bounds))
		}
		ll_b$loglik=sapply(ll_b$b,FUN= loglik_b_temp)
		if (sum(is.nan(ll_b$loglik))>0) {ll_b=ll_b[-which(is.nan(ll_b$loglik)),]}
		ll_cond3=exp(ll_b[,2])/sum(exp(ll_b[,2]))
    	CDF_right=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    	CDF_left=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    	CI95_b=c(ll_b[min(which(CDF_right>0.025)),1], ll_b[max(which(CDF_left>0.025)),1])
    	plot(round(ll_b,digits=4), main='log-likelihood vs b (x^2 term)',type='l') ; abline(v=fit$par$b,col=2)
	}
	else {CI95_b=NA}

	# c
	if ('c'%in%names(fit$par)){
		ll_c=as.data.frame(matrix(NA, effort_uncertainty,2))
    	colnames(ll_c)=c('c','loglik')
    	ll_c$c =sort(seq(from=scope_c[1],to=scope_c[2], length.out=effort_uncertainty))
		loglik_c_temp=function(x){
			V=a*SEQ^4+b*SEQ^2+x*SEQ
			dMat=DiffMat_backwards(V)
			return(LogLik_bounds(tree_formatted,log(fit$par$sigsq/2),dMat,bounds=bounds))
		}
		ll_c$loglik=sapply(ll_c$c,FUN= loglik_c_temp)
		if (sum(is.nan(ll_c$loglik))>0) {ll_c=ll_c[-which(is.nan(ll_c$loglik)),]}
		ll_cond3=exp(ll_c[,2])/sum(exp(ll_c[,2]))
    	CDF_right=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(1:x)])})
    	CDF_left=sapply(1:length(ll_cond3),FUN=function(x){sum(ll_cond3[c(x:length(ll_cond3))])})
    	CI95_c=c(ll_c[min(which(CDF_right>0.025)),1], ll_c[max(which(CDF_left>0.025)),1])
    	plot(round(ll_c,digits=4), main='log-likelihood vs c (x term)',type='l') ; abline(v=fit$par$c,col=2)
	}
	else {CI95_c=NA}
	
	# x0
	V=a*SEQ^4+b*SEQ^2+c*SEQ
	dMat=DiffMat_backwards(V)
	ll_root=as.data.frame(matrix(NA, Npts,2))
    colnames(ll_root)=c('root','density')
    ll_root$root =seq(bounds[1],to=bounds[2], length.out= Npts)
    tree_formatted2= tree_formatted
	pMat=prep_mat_exp(dCoeff=log(fit$par$sigsq/2),dMat,bounds=bounds) # edited
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
    # TO EDIT WHEN ROOT VALUE IS ADDED IN THERE
    #if (fit$par$root_value==bounds[1]){cells_root=c(1,max(which(CDF_root_left>0.05)))}
    #if (fit$par$root_value== bounds[2]){cells_root=c(min(which(CDF_root_right>0.05)),length(ll_root$density))}
	CI95_root=bounds[1]+(bounds[2]-bounds[1])/(Npts-1)*(cells_root-1)
	plot(round(ll_root,digits=4), main='density vs root trait value',type='l') ; abline(v=fit$par$root_value,col=2)
	
		
	return(list(CI95_sig2= CI95_sig2, CI95_a= CI95_a, CI95_b= CI95_b, CI95_c= CI95_c, CI95_root= CI95_root))
}