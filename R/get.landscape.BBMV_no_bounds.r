# This function plots the adaptive landscape estimated by the FPK and BBM+V models.
# It takes a model fitted by the function 'find.mle_FPK' as its main argument
# Npts determines the number of points used to plot the adaptive landscape or potential
get.landscape.FPK=function(fit,Npts=100,main='Macroevolutionary landscape',ylab='N.exp(-V)',xlab='Trait',xlim=NULL,ylim=NULL){
	# retrieve coefficients for the potential
	if ('a'%in%names(fit$par)){a=fit$par$a}
  else {a=fit$par_fixed$a}
  if ('b'%in%names(fit$par)){b=fit$par$b}
  else {b=fit$par_fixed$b}
  if ('c'%in%names(fit$par)){c=fit$par$c}
  else {c=fit$par_fixed$c}
	# build potential and stationary distribution of the trait
  bounds=fit$par_fixed$bounds
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
	step=(bounds[2]-bounds[1])/(Npts-1)
	if (is.null(xlim)){xlim=c(bounds[1],bounds[2])}
	if (is.null(ylim)){ylim=c(0,max(exp(-V)/sum(exp(-V)*step)))}
plot((exp(-V)/sum(exp(-V)*step))~seq(from=bounds[1],to=bounds[2],length.out=Npts),type='l',col=2,lwd=3,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
}