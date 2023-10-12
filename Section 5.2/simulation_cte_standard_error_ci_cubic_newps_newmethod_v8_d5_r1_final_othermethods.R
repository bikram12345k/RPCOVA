## Simulation for Section 5.2 the competing methods for dimension d = 5
args = commandArgs(trailingOnly=TRUE)
id=as.numeric(args)[1]

dseq = 5 #c(5,5,5,5,5,5,5,5,5)
seedseq = as.numeric(sapply((1:5), function(x) x*c(432, 140, 645, 236, 513, 142, 2635, 3616, 6231, 6216, 321, 408, 451)))
nseq = c(2000, 4000, 6000, 8000, 10000, 12000)


foo <- function(x, y) c(x, y)
setup = outer(dseq, seedseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup, nseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
#setup = outer(setup, cscaleseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)


library(Matrix, lib.loc='/blue/bkarmakar/bkarmakar/struncatedP/')

#library(Matrix)

#install.packages('ks')
library(ks)

######################
n = setup[3,id]; d = setup[1,id]; thisseed = setup[2,id]; #cscale = setup[4,id]

set.seed(thisseed)
print(c(n, d, thisseed)) #, cscale

######################

library(mvtnorm)
library(KernSmooth)
library(neuralnet)
library(slam)
library(gurobi)	# See directions for installation in the document
library(np)
library(mclust)

library(grf)
#library(bcf)
library(causalToolbox)
library(BART, quietly = TRUE)

source('supporting_fn.R')

deprintize<-function(f){
 return(function(...) {capture.output(w<-f(...));return(w);});
}	
rm_outliers <- function(xtemp){
	uq = quantile(xtemp, .995, na.rm = TRUE)
	lq = quantile(xtemp, .005, na.rm = TRUE)
	xtemp_out1 = which(xtemp<=lq | xtemp >=uq)
	xtemp_out2 = which(xtemp %in% boxplot.stats(xtemp)$out)
	
	xtemp_out = intersect(xtemp_out1, xtemp_out2)
	xtemp[-xtemp_out]
}

	
	##########################
	
	#newx = x	# evaluation points
	
	xiseq = seq(1,5,.1*d)
	if(d>3) xiseq = seq(1.2,3.8,.5) 
	if(d>1){
		newx = outer(xiseq, xiseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
		if(d>2)
		for(i in 1:(d-2))
			newx = outer(newx, xiseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
	}
	if(d==1){ 
		newx = matrix(xiseq,ncol=1)
	} else newx = t(sapply(newx, cbind))
	
	set.seed(0841)	# Same evaluation points across iterations
	newx = newx[seq(1,nrow(newx),4),]
	newx = newx + matrix(rnorm(nrow(newx)*d, 0,.2), ncol=d)

sigma02tr = var(runif(n*2,0,2))/4+.25*(newx[,1]-3)^2+1

	newx = newx[sigma02tr>1.3,]
	set.seed(thisseed)
	
#for(itr in 1:40){
for(itr in 1:5){
	print(itr)

	ok=0
	while(!ok){
		# unobserved confounder
		u = runif(n,0,2)
		x = rmvnorm(n, rep(3,d)) #3 + matrix(runif(n*d, -1, 1), nrow=n, ncol=d)

		#ps = 1/(1+exp(-(.5*rowMeans(x)-u+.5))) 
		ps=1/(1+exp(-(.5*rowMeans(x)+2*u-3)))#*( 1/(1+exp(-2*u+2)) )
		summary(ps)
		z = rbinom(n,1, ps)

		y0 = 1+u+rowSums(x)+rnorm(n,0,.5)*(x[,1]-3)+rnorm(n,0,1)

		if(d==1)
			taufn <- function(x) .5*(x[1]-3)^2 #.5*x[1]
		if(d>1)
			taufn <- function(x) .5*(x[1]-3)+mean(x[-1])/3 #.5*x[1]+mean(x[-1])/3
		tau = apply(x, 1, taufn)

		y = y0 + tau*z
		#nnfit = deprintize(neuralnet::neuralnet)(z~., data=data.frame(z,x), hidden=4, linear.output=TRUE)
		#nnfit = bartMachine(data.frame(x), as.factor(z), verbose=FALSE)
		psfit = deprintize(pbart)(data.frame(x), z)


		ok=1
		#if( !("weights" %in% names(nnfit)) )
		#	ok = 0
		
	}
	
	yvec = as.numeric(y)
	xdat = data.frame(x)
	

	##################################################################
	## Causal forest  
	c.forest <- causal_forest(xdat, yvec, z, num.trees = 4000)
	c.pred <- predict(c.forest, newdata = data.frame(newx), estimate.variance = TRUE)

	res_cf = cbind(c.pred$predictions, sqrt(c.pred$variance.estimates), 
			(c.pred$predictions - 1.96*sqrt(c.pred$variance.estimates) < tau) &
				(tau < c.pred$predictions + 1.96*sqrt(c.pred$variance.estimates)))
	print(c('track =', 6))
	
	##################################################################
	## Bayesian causal tree: not relavent for CI
	#pixnn = predict(nnfit, newdata=xdat)[,1]
	#bcf_fit = bcf(yvec, z, as.matrix(xdat), as.matrix(xdat), 
	#		pixnn, nburn=1000, nsim=1000)
	# Get posterior of treatment effects
	#tau_post = bcf_fit$tau
	#tauhat = colMeans(tau_post)
	#mapnewx = sapply(newx, findInterval, vec=sort(x))
	#tauhatnewx = (tauhat[order(as.numeric(x))])[mapnewx]

	##################################################################
	## X-Learner
	if(d==1){
		xl_rf <- X_RF(feat = data.frame(x=x, x2=x^2), tr = z, yobs = yvec)
		cate_esti_xrf <- EstimateCate(xl_rf, data.frame(x=newx, x2=newx^2))
		#xl_ci_rf <- CateCI(xl_rf, data.frame(x=newx, x2=newx^2), B = 500)
	} else {
		xl_rf <- X_RF(feat = data.frame(x=x), tr = z, yobs = yvec)
		cate_esti_xrf <- EstimateCate(xl_rf, data.frame(x=newx))
		#xl_ci_rf <- CateCI(xl_rf, data.frame(x=newx), B = 500)
	}
	print(c('track =', 7))
	## bootstrap ci
	
	#boot_ci = xl_ci_rf
	
	rm(xl_rf)

	boot_cate = NULL
	for(b in 1:200){
		sampb = sample(1:n, n)
		if(d==1){
			xl_rfb <- X_RF(feat = data.frame(x=x[sampb], x2=x[sampb]^2), tr = z[sampb], yobs = yvec[sampb])
			
			#xl_ci_rf <- CateCI(xl_rf, feature_test, B = 500)
			
			cate_esti_xrfb <- EstimateCate(xl_rfb, data.frame(x=newx, x2=newx^2))
		} else {
			xl_rfb <- X_RF(feat = data.frame(x=x[sampb,]), tr = z[sampb], yobs = yvec[sampb])
			cate_esti_xrfb <- EstimateCate(xl_rfb, data.frame(x=newx))
		}
		boot_cate = cbind(boot_cate, cate_esti_xrfb)
		rm(xl_rfb)
		gc()
	}
	print(c('track =', 8))
	boot_ci = apply(boot_cate, 1, quantile, c(.025, .975))
	
	res_xrf = cbind(cate_esti_xrf, boot_ci[1,], boot_ci[2,], (boot_ci[1,] < tau) & (tau < boot_ci[2,]))


	## Results file
	thisrun = 0
	fname = paste0('/blue/bkarmakar/bkarmakar/RPCOVA_ci/new_ll_fixed/simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw_othermethods.csv')
	while(file.exists(fname)){
		thisrun = thisrun + 1
		fname = paste0('/blue/bkarmakar/bkarmakar/RPCOVA_ci/new_ll_fixed/simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw_othermethods.csv')
	}
	## Write results
	write.csv( cbind(n, d, res_cf, res_xrf),	#res_proposal), #
		file = fname, row.names = FALSE)

	rm(boot_ci, boot_cate, cate_esti_xrfb, xl_rfb, cate_esti_xrf, xl_rf,
			c.forest, c.pred, psfit, fx, term1, term2, term3, mx1, mx0,
			tauplusest,	tauminusest, setauplus,	setauminus, Delta2x, Delta2x1,
			cutoffn,seplus,seminus,Delta2xeps,Delta2x2,
			rx0, rx1, v02, v12, seDelta2x, pixnn,
			sigma4lambda2, sigma2nu0eta0, sigma2nu1eta1,sebetaU, mx
			
			)
	gc()
}




#install.packages('rgl')
#library(rgl)
#plot3d(x=x[,1], y=x[,2], z=tauminus) 
