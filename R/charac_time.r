# Calculate the characteristic time it takes for the process to reach its stationary distribution
charac_time=function(Npts=100,fit){
	# retrieve coefficients for the potential
	if ('a'%in%names(fit$par)){a=fit$par$a}
	else {a=fit$par_fixed$a}
	if ('b'%in%names(fit$par)){b=fit$par$b}
	else {b=fit$par_fixed$b}
	if ('c'%in%names(fit$par)){c=fit$par$c}
	else {c=fit$par_fixed$c}
  bounds=fit$par_fixed$bounds
	# build potential and get second largest eigenvalue
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	Vq=a*SEQ^4+b*SEQ^2+c*SEQ
	Mat=DiffMat_forward(Vq)
	vp=diag(Mat$diag)
	Tc=(2*(bounds[2]-bounds[1])^2/fit$par$sigsq)/(Npts-1)^2/abs(sort(Re(vp),decreasing=T)[2]) # the first one is 0 (or something very close due to numerical approximations)
	return(Tc)
}
