# This function plots the adaptive landscape estimated by the BBM+V model.
# It takes a model fitted by BBMV as its main argument
# Npts determines the number of points used to plot the adaptive landscape or potential
plot.landscape.BBMV=function(model,Npts=50,main='Adaptive landscape',ylab='-V',xlab='Trait',xlim=NULL,ylim=NULL){
	# determine number of parameters of the model
	npar=length(model[[1]]) # 3: flat ; 4: linear ; 5: quadratic ; 6: x4
	# retrieve coefficients for the potential
	if (npar==3){a=b=c=0}
	if (npar==4){a=b=0 ; pars=model[[1]][-which(names(model[[1]])%in%c('bounds','sigsq','root_value'))] ; c=pars[[1]]}
	if (npar==5){a=0 ; pars=model[[1]][-which(names(model[[1]])%in%c('bounds','sigsq','root_value'))] ; b=pars[[1]] ; c=pars[[2]]}
	if (npar==6){a=model[[1]]$a ; b=model[[1]]$b ; c=model[[1]]$c}
	# build potential and stationary distribution of the trait
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
	bounds=model[[1]]$bounds
	if (is.null(xlim)){xlim=c(bounds[1],bounds[2])}
	if (is.null(ylim)){ylim=c(min(-V),max(-V))}
plot(-V~seq(from=bounds[1],to=bounds[2],length.out=Npts),type='l',col=2,lwd=3,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
}

# Same function than above but allows for plotting all adaptive landscapes from a list of models fitted. This is rather experimental... and probably not very useful anyway.
plot.multiple.landscapes.BBMV=function(models,Npts=50,ylim=c(0,0.2),main='Adaptive landscapes',ylab='-V',xlab='Trait'){
	trm=c()
	for (i in 1:length(models)){if (class(models[[i]])=='try-error'){trm=c(trm,i)}}
	if (length(trm)>0){models=models[-trm]} 
	model=models[[1]]	
	# determine number of parameters of the model and retrieve coefficients
	npar=length(model[[1]]) # 3: flat ; 4: linear ; 5: quadratic ; 6: x4
	if (npar==3){a=b=c=0}
	if (npar==4){a=b=0 ; pars=model[[1]][-which(names(model[[1]])%in%c('bounds','sigsq','root_value'))] ; c=pars[[1]]}
	if (npar==5){a=0 ; pars=model[[1]][-which(names(model[[1]])%in%c('bounds','sigsq','root_value'))] ; b=pars[[1]] ; c=pars[[2]]}
	if (npar==6){a=model[[1]]$a ; b=model[[1]]$b ; c=model[[1]]$c}
	# build potential and stationary distribution of the trait
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
	bounds=model[[1]]$bounds
plot(-V~seq(from=bounds[1],to=bounds[2],length.out=Npts),type='l',col=2,lwd=1,main=main,xlab=xlab,ylab=ylab,ylim=ylim)
	for (i in 2:length(models)){
	model=models[[i]]	
	npar=length(model[[1]]) # 3: flat ; 4: linear ; 5: quadratic ; 6: x4
	# retrieve coefficients for the potential
	if (npar==3){a=b=c=0}
	if (npar==4){a=b=0 ; pars=model[[1]][-which(names(model[[1]])%in%c('bounds','sigsq','root_value'))] ; c=pars[[1]]}
	if (npar==5){a=0 ; pars=model[[1]][-which(names(model[[1]])%in%c('bounds','sigsq','root_value'))] ; b=pars[[1]] ; c=pars[[2]]}
	if (npar==6){a=model[[1]]$a ; b=model[[1]]$b ; c=model[[1]]$c}
	# build potential and stationary distribution of the trait
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
	bounds=model[[1]]$bounds
lines(-V~seq(from=bounds[1],to=bounds[2],length.out=Npts),col=2)
	}
}