# still problem with the transformation in sigsq
posterior_vs_prior=function(chain,param='a',burnin=0.2,type_prior='Normal',shape_prior=c(0,10)){
  if (type_prior=='Normal'){
  empirical_density=rnorm(n=10000,mean=shape_prior[1],sd=shape_prior[2])
  }
  if (type_prior=='Uniform'){
    empirical_density=runif(n=10000,min=shape_prior[1],max=shape_prior[2])
  }
  prior=density(empirical_density) ; posterior=density(chain[-c(1:floor(dim(chain)[1]*burnin)),param])
  if (param=='sigsq'){posterior=density(log(chain[-c(1:floor(dim(chain)[1]*burnin)),param]/2))}
  plot(prior,col='grey55',main='Prior vs. posterior comparison',lwd=0.4,xlab='',ylim=c(0,max(c(prior$y,posterior$y))),xlim=c(min(c(prior$x,posterior$x)),max(c(prior$x,posterior$x))))
  polygon(prior,col=adjustcolor(col='grey',alpha.f=0.5),border='grey55')
  polygon(posterior,col=adjustcolor(col='red',alpha.f=0.5),border='red')
  legend(legend=c('prior','posterior'),fill=c(adjustcolor(col='grey',alpha.f=0.5),adjustcolor(col='red',alpha.f=0.5)),x='topright')
}