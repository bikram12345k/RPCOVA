
u = runif(10000,0,2)	# unobserved confounder
		x = rmvnorm(10000, rep(3,d))

		#z = rbinom(n,1, 1/(1+exp(-(.4*u+.5*rowMeans(x)-1.5))))
		z = rbinom(10000,1, 1/(1+exp(-(u/2+.5*rowMeans(x)-1.5))))
		y0 = 1+u^2+d*rowMeans(x)+rnorm(1000,0,1)*rowMeans(x-3)/6
		if(d==1)
			taufn <- function(x) .5*x[1]
		if(d>1)
			taufn <- function(x) .5*x[1]+mean(x[-1])/3
		tau = apply(x, 1, taufn)
		#tau = rowSums(taufn(x))
		
		y1 = y0 + (tau+rnorm(10000,0,.25))
		y0 = y0 + rnorm(10000,0,.25)
		
## Correlation calculations
library(ppcor)	
pcor(data.frame(z, u, x), 'kendall')

library(ppcor)	
pcor(data.frame(y0, u, x), 'kendall')

mean(y1-y0)
mean((y1-y0)[z==1])