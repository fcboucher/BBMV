# Calculate the characteristic time it takes for the process to reach its stationary distribution
charac_time=function(Npts,model){ # X=c(Npts,bound,a,b,c,sigma)
	npar=length(model$par) # 3: flat ; 4: linear ; 5: quadratic ; 6: x4
	# retrieve coefficients for the potential
	if (npar==3){a=b=c=0}
	if (npar==4){a=b=0 ; c=model$par$c}
	if (npar==5){a=0 ; b=model$par$b ; c=model$par$c}
	if (npar==6){a=model$par$a ; b=model$par$b ; c=model$par$c}
	# build potential and stationary distribution of the trait
	SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
	V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
	Tc=exp(max(V)-min(V))*(model$par$bounds[2]-model$par$bounds[1])^2/model$par$sigsq
	return(Tc)
}
