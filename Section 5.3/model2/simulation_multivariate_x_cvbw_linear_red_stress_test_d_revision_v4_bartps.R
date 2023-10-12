#############
######## Simulation for ATT using the proposed method
#############

args = commandArgs(trailingOnly=TRUE)
id=as.numeric(args)[1]


seedseq <- 2*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721) #, 4241, 1921, 5693, 1218, 8216) 
seedseq <- c(seedseq, 3*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))
seedseq <- c(seedseq, 4*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))
seedseq <- c(seedseq, 5*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))

dseq = 2:5
nseq = c(500, 1000, 1500, 2000)

foo <- function(x, y) c(x, y)
setup = outer(seedseq , dseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup , nseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)


thisseed = setup[1,id]
set.seed(thisseed)
d = setup[2,id]
n = setup[3,id]

print(c(n,d,thisseed))

Itr = 10 #25 #100 #60

library(Matrix)
library(mvtnorm)
library(KernSmooth)
library(neuralnet)
library(slam)
library(gurobi)	# See directions for installation in the document
library(np, lib.loc="/blue/bkarmakar/bkarmakar/struncatedP/oldversion/")

library(grf)
#library(bcf)
library(causalToolbox)
library(BART, quietly = TRUE)

source('supporting_fn.R')

deprintize<-function(f){
 return(function(...) {capture.output(w<-f(...));return(w);});
}	
rm_outliers <- function(xtemp){
	uq = quantile(xtemp, .999, na.rm = TRUE)
	lq = quantile(xtemp, .000, na.rm = TRUE)
	xtemp_out1 = which(xtemp<lq | xtemp >uq)
	xtemp_out2 = which(xtemp %in% boxplot.stats(xtemp)$out)
	
	xtemp_out = intersect(xtemp_out1, xtemp_out2)
	if(length(xtemp_out)>0){
		return(xtemp[-xtemp_out])
	} else return(xtemp)
}

#####################################
## CALCULATE TRUE ATT
		u = runif(1000000,0,2)	# unobserved confounder
		x = rmvnorm(1000000, rep(3,d))

		#z = rbinom(n,1, 1/(1+exp(-(.4*u+.5*rowMeans(x)-1.5))))
		z = rbinom(1000000,1, 1/(1+exp(-(u/2+.5*rowMeans(x)-1.5))))
		y0 = 1+u^2+d*rowMeans(x)+rnorm(1000000,0,1)*rowMeans(x-3)/6
		if(d==1)
			taufn <- function(x) .5*x[1]
		if(d>1)
			taufn <- function(x) .5*x[1]+mean(x[-1])/3
		tau = apply(x, 1, taufn)
		#tau = rowSums(taufn(x))
		
		y1 = y0 + (tau+rnorm(1000000,0,.25))
		y0 = y0 + rnorm(1000000,0,.25)
		
		mean( (y1-y0))
		mean( (y1-y0)[z==1])
		mean( (y1-y0)[z==0])
		
		att = mean( (y1-y0)[z==1])

rm(x, u, z, y1, y0, tau)
####################################


res = NULL
res_mse = NULL

for(itr in 1:Itr){


	# Data generating model
	ok=0
	while(!ok){
		u = runif(n,0,2)	# unobserved confounder
		x = rmvnorm(n, rep(3,d))

		#z = rbinom(n,1, 1/(1+exp(-(.4*u+.5*rowMeans(x)-1.5))))
		z = rbinom(n,1, 1/(1+exp(-(u/2+.5*rowMeans(x)-1.5))))
		y0 = 1+u^2+d*rowMeans(x)+rnorm(n,0,1)*rowMeans(x-3)/6
		if(d==1)
			taufn <- function(x) .5*x[1]
		if(d>1)
			taufn <- function(x) .5*x[1]+mean(x[-1])/3
		tau = apply(x, 1, taufn)
		#tau = rowSums(taufn(x))
		y = rep(0, n)
		
		y[z==1] = y0[z==1] + (tau[z==1]+rnorm(sum(z),0,.25))
		y[z==0] = y0[z==0] + rnorm(sum(1-z),0,.25)
		
		#nnfit = deprintize(neuralnet::neuralnet)(z~., data=data.frame(z,x), hidden=8)
		psfit = deprintize(pbart)(data.frame(x), z)
		ok=1
		#if( !("weights" %in% names(nnfit)) )
		#	ok = 0
	}

	pixnn = colMeans(deprintize(predict)(psfit, data.frame(x))$prob.test)
	# Truncate away from 0 and 1
	pixnn = pmax(pixnn, .015)
	pixnn = pmin(pixnn, .985)
	#hist(pixnn)

	## Resistant POpulation VAriance (RPOVA)
	sigma02 = (rowMeans(x-3))^2/36+16*(4/45) + .25^2 #1+16*(4/45)

	## Estimate pi(x)
	# there is some issue with convergence sometimes(!?)
	newx = x # rmvnorm(2*n, rep(3,d)) #seq(1,5,.025)

	# estimate E(y | x)
	yvec = as.numeric(y)
	xdat = data.frame(x)
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat,regtype='ll')
	mx <- deprintize(npreg)(h1)

	# estimate Var(y | x)
	rx = yvec - mx$mean
	
	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat, bws=h1$bw,regtype='ll',bandwidth.compute = FALSE)
	sigma2x = deprintize(npreg)(h2, exdat=data.frame(newx))

	
	Delta2x2 = (sigma2x$mean - sigma02)/( pixnn*(1-pixnn) )

	hc = deprintize(npregbw)(ydat=yvec[z==0],xdat=xdat[z==0,], bws=h2$bw, bandwidth.compute = TRUE,regtype='ll')#$bw
	ht = deprintize(npregbw)(ydat=yvec[z==1],xdat=xdat[z==1,], bws=h2$bw, bandwidth.compute = TRUE,regtype='ll')#$bw

	mx0 = deprintize(npreg)(hc, exdat=data.frame(newx))
	mx1 = deprintize(npreg)(ht, exdat=data.frame(newx))


	# run constraint optimization
	Delta2x1 = sapply(1:nrow(newx), function(id) 
			deprintize(const_loclin2_d)(newx[id,], h=ht$bw, hs=hc$bw, x, d, const=Delta2x2[id], y))

	if(is.list(Delta2x1))
		Delta2x1 = sapply(Delta2x1, function(x) ifelse(length(x)==0, NA, x))
	Delta2x = pmax(Delta2x1^2 - Delta2x2, 0)


	# Calculate MSE
	res_mse_temp = c( (mean(rm_outliers(Delta2x1 - sqrt(Delta2x)), na.rm = TRUE)-att)^2 ,
		(mean(rm_outliers(Delta2x1 + sqrt(Delta2x)), na.rm = TRUE)-att)^2,
		(mean(rm_outliers(Delta2x1), na.rm = TRUE)-att)^2 )
	
	res_mse_temp
	
	restemp = rbind(cbind(newx,NA,NA),
		c(tau, NA, NA),
		c(Delta2x1 - sqrt(Delta2x), NA, res_mse_temp[1]),
		c(Delta2x1 + sqrt(Delta2x), NA, res_mse_temp[2]),
		c(Delta2x1, NA, res_mse_temp[3]) )
		
	res = rbind(res, c(newx,NA,
		Delta2x1 - sqrt(Delta2x), NA, 
		Delta2x1 + sqrt(Delta2x), NA, 
		Delta2x1, NA, NA, res_mse_temp ))
	
	
	res_mse = rbind(res_mse, res_mse_temp)
	
	print(c(itr, res_mse_temp ))
	
	## Write results.
	write.csv(cbind(thisseed,d,n,restemp), file=paste0('linear_stresstest_d_',n,'_',d,'_',thisseed,'_',itr,'_revision_v4_bartps.csv'), row.names = FALSE, col.names = FALSE)

}

print('final')
print(cbind(thisseed,d,n,res_mse))
write.csv(cbind(thisseed,d,n,res), file=paste0('linear_stresstest_d_',n,'_',d,'_',thisseed,'_revision_v4_bartps.csv'), append=TRUE, row.names = FALSE, col.names = FALSE)
