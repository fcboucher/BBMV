# This function plots the adaptive landscape estimated by the BBM+V model.
# It takes a model fitted by BBMV as its main argument
# if landscape=TRUE (default), then the adaptive landscape is plotted
# if landscape=FALSE, then the potential is plotted instead
# Npts determines the number of points used to plot the adaptive landscape or potential
plot.landscape.BBMV=function(model,landscape=T,Npts=50,main='Adaptive landscape'){
	# determine number of parameters of the model
	npar=length(model[[1]]) # 3: flat ; 4: linear ; 5: quadratic ; 6: x4
	# retrieve coefficients for the potential
	if (npar==3){a=b=c=0}
	if (npar==4){a=b=0 ; c=model[[1]]$c}
	if (npar==5){a=0 ; b=model[[1]]$b ; c=model[[1]]$c}
	if (npar==6){a=model[[1]]$a ; b=model[[1]]$b ; c=model[[1]]$c}
	# build potential and stationary distribution of the trait
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
	L=exp(-V)/sum(exp(-V)) # landscape = stationary distribution = exp(-V) normalized
	bounds=model[[1]]$bounds
	if (landscape==T){
		plot(L~seq(from=bounds[1],to=bounds[2],length.out=Npts),type='l',col=2,lwd=3,main=main,xlab='Trait')
	}
	else {
plot(V~seq(from=bounds[1],to=bounds[2],length.out=Npts),type='l',col=2,lwd=3,main='Potential',xlab='Trait')
	}
}

# Same function than above but allows for plotting all adaptive landscapes from a list of models fitted. This is rather experimental... and probably not very useful anyway.
plot.multiple.landscapes.BBMV=function(models,landscape=T,Npts=50,ylim=c(0,0.2),main='Adaptive landscapes'){
	trm=c()
	for (i in 1:length(models)){if (class(models[[i]])=='try-error'){trm=c(trm,i)}}
	if (length(trm)>0){models=models[-trm]} 
	model=models[[1]]	
	# determine number of parameters of the model and retrieve coefficients
	npar=length(model[[1]]) # 3: flat ; 4: linear ; 5: quadratic ; 6: x4
	if (npar==3){a=b=c=0}
	if (npar==4){a=b=0 ; c=model[[1]]$c}
	if (npar==5){a=0 ; b=model[[1]]$b ; c=model[[1]]$c}
	if (npar==6){a=model[[1]]$a ; b=model[[1]]$b ; c=model[[1]]$c}
	# build potential and stationary distribution of the trait
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
	L=exp(-V)/sum(exp(-V)) # landscape = stationary distribution = exp(-V) normalized
	bounds=model[[1]]$bounds
	if (landscape==T){
		plot(L~seq(from=bounds[1],to=bounds[2],length.out=Npts),type='l',col=2,lwd=1,main=main,xlab='Trait',ylim=ylim)
	}
	else {
plot(V~seq(from=bounds[1],to=bounds[2],length.out=Npts),type='l',col=2,lwd=1,main='Potential',xlab='Trait',ylim=ylim)
	}
	for (i in 2:length(models)){
	model=models[[i]]	
	npar=length(model[[1]]) # 3: flat ; 4: linear ; 5: quadratic ; 6: x4
	# retrieve coefficients for the potential
	if (npar==3){a=b=c=0}
	if (npar==4){a=b=0 ; c=model[[1]]$c}
	if (npar==5){a=0 ; b=model[[1]]$b ; c=model[[1]]$c}
	if (npar==6){a=model[[1]]$a ; b=model[[1]]$b ; c=model[[1]]$c}
	# build potential and stationary distribution of the trait
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
	L=exp(-V)/sum(exp(-V)) # landscape = stationary distribution = exp(-V) normalized
	bounds=model[[1]]$bounds
	if (landscape==T){
		lines(L~seq(from=bounds[1],to=bounds[2],length.out=Npts),col=2)
	}
	else {
lines(V~seq(from=bounds[1],to=bounds[2],length.out=Npts),col=2)
	}
	}
}