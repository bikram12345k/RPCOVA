args = commandArgs(trailingOnly=TRUE)
id=as.numeric(args)[1]


dseq = c(3,3,3,4,4,4,5,5,5)
nseq = c(2000, 3000, 4000, 5000)
#seedseq = c(1234, 8041, 1546, 1362, 2315, 2141, 5362, 6163, 2326, 6126)
#seedseq = c(seedseq, c(123, 804, 154))
seedseq = c(4321, 1408, 6451, 2362, 5132, 1421, 2635, 3616, 6231, 6216, 321, 408, 451)

foo <- function(x, y) c(x, y)
setup = outer(dseq, nseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup, seedseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)


library(Matrix, lib.loc='/blue/bkarmakar/bkarmakar/struncatedP/newversion')

#library(Matrix)

#install.packages('ks')
library(ks)

######################
n = setup[2,id]; d = setup[1,id]; thisseed = setup[3,id]

set.seed(thisseed)
print(c(n, d, thisseed))

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


for(itr in 1:40){
	print(itr)

	ok=0
	while(!ok){
		u = runif(n,0,2)	# unobserved confounder
		x = rmvnorm(n, rep(3,d))

		ps = 1/(1+exp(-(.75*rowMeans(x))))#*( 1/(1+exp(-2*u+2)) )
		summary(ps)
		z = rbinom(n,1, ps)

		y0 = 1+u^2+rowSums(x)+rnorm(n,0,1)

		if(d==1)
			taufn <- function(x) .5*(x[1]-3)^2 #.5*x[1]
		if(d>1)
			taufn <- function(x) .5*(x[1]-3)^2 +mean(x[-1])/3 #.5*x[1]+mean(x[-1])/3
		tau = apply(x, 1, taufn)

		y = y0 + tau*z
		nnfit = deprintize(neuralnet::neuralnet)(z~., data=data.frame(z,x), hidden=4, linear.output=TRUE)
		ok=1
		if( !("weights" %in% names(nnfit)) )
			ok = 0
		
	}
	#Sys.time()
	#pixnn = predict(nnfit, newdata=data.frame(x))[,1]
	
	print(c('track =', 1))
	
	## Resistant Population Calibration Of VAriance (RPCOVA)
	sigma02 = 1+var(runif(3*n,0,2)^2)

	#newx = x	# evaluation points
	
	xiseq = seq(1,5,.1*d)
	if(d>3) xiseq = seq(1,5,.2*d)
	if(d>1){
		newx = outer(xiseq, xiseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
		if(d>2)
		for(i in 1:(d-2))
			newx = outer(newx, xiseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
	}
	if(d==1){ 
		newx = matrix(xiseq,ncol=1)
	} else newx = t(sapply(newx, cbind))

	#dim(newx)
	
	## Estimate pi(x)
	pixnn = predict(nnfit, newdata=data.frame(newx))[,1]
	# Truncate away from 0 and 1
	pixnn = pmax(pixnn, .015)
	pixnn = pmin(pixnn, .985)
	#hist(pixnn)
	
	yvec = as.numeric(y)
	xdat = data.frame(x)
	
	#### Estimation	####
	# Step 1
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat,regtype='lc')
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat,bws=(h1$bw)^(1.1), bandwidth.compute = FALSE,regtype='ll')

	mx <- deprintize(npreg)(h1)

	# Step 2
	rx = yvec - mx$mean

	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat,regtype='lc')
	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat, bws=(h2$bw)^(1.1), bandwidth.compute = FALSE,regtype='ll') #undersmooth
	sigma2x = deprintize(npreg)(h2, exdat=data.frame(newx))

	# Step 3
	Delta2x2 = (sigma2x$mean - sigma02)/( pixnn*(1-pixnn) )

	hc = deprintize(npregbw)(ydat=yvec[z==0],xdat=xdat[z==0,], bws=h2$bw, bandwidth.compute = TRUE,regtype='lc')
	hc = deprintize(npregbw)(ydat=yvec[z==0],xdat=xdat[z==0,], bws=hc$bw^(1.1), bandwidth.compute = FALSE,regtype='ll')	#undersmooth

	ht = deprintize(npregbw)(ydat=yvec[z==1],xdat=xdat[z==1,], bws=h2$bw, bandwidth.compute = TRUE,regtype='lc')
	ht = deprintize(npregbw)(ydat=yvec[z==1],xdat=xdat[z==1,], bws=ht$bw^(1.1), bandwidth.compute = FALSE,regtype='ll')	#undersmooth
	
	# run constraint optimization
	Delta2x1 = sapply(1:nrow(newx), function(id) {			
			if(id %% 500 == 0) {
				print(id)
				deprintize(gc)()}
			deprintize(const_loclin2_d)(newx[id,], h=ht$bw, hs=hc$bw, x, d, const=  Delta2x2[id], y)
			})
	deprintize(gc)()
	print(c('track =', 2))	

	# Step 4
	if(is.list(Delta2x1))
		Delta2x1 = sapply(Delta2x1, function(x) ifelse(length(x)==0, NA, x))
	Delta2x = pmax(Delta2x1^2 - Delta2x2, 0)
	print(c('track =', 3))
	
	
	#### Calculating terms for the standard errors	####
	
	mx0 = deprintize(npreg)(hc, exdat=data.frame(x[z==0,]))
	mx1 = deprintize(npreg)(ht, exdat=data.frame(x[z==1,]))


	rx0 = yvec[z==0] - mx0$mean
	rx1 = yvec[z==1] - mx1$mean

	v02 <- deprintize(npregbw)(ydat=rx0^2,xdat=xdat[z==0,], bws=h2$bw, bandwidth.compute = TRUE)
	v02 = deprintize(npreg)(v02, exdat=data.frame(newx))

	v12 <- deprintize(npregbw)(ydat=rx1^2,xdat=xdat[z==1,], bws=h2$bw, bandwidth.compute = TRUE)
	v12 = deprintize(npreg)(v12, exdat=data.frame(newx))



	sigma4lambda2 = deprintize(npregbw)(ydat=rx^4,xdat=xdat, regtype = 'lc')
	sigma4lambda2 = deprintize(npreg)(sigma4lambda2, exdat=data.frame(newx))$mean
	sigma4lambda2 = pmax(sigma4lambda2 - sigma2x$mean^2, 0)

	sigma2nu0eta0 = deprintize(npregbw)(ydat=rx[z==0]^2*rx0,xdat=xdat[z==0,], regtype = 'lc')
	sigma2nu0eta0 = deprintize(npreg)(sigma2nu0eta0, exdat=data.frame(newx))$mean

	sigma2nu1eta1 = deprintize(npregbw)(ydat=rx[z==1]^2*rx1,xdat=xdat[z==1,], regtype = 'lc')
	sigma2nu1eta1 = deprintize(npreg)(sigma2nu1eta1, exdat=data.frame(newx))$mean

	print(c('track =', 4))
	
	## Estimate density
	# By ks
	#if(d <= 4){
	#	fx = kde(x = x, h=hpi(x), eval.points = newx, density=TRUE)
	#} else fx = kde(x = x, h=hpi(x), eval.points = newx, binned=TRUE, density=TRUE)
	#fx = fx$estimate 
	
	# By gaussian mixture model
	fx = deprintize(densityMclust)(x, plot=FALSE)
	fx = predict(fx, newdata=newx)
	

	
	#### Variance calculations  ####
	
	# se of tauplus and tauminus
	thetaK = (1/(2*sqrt(pi)))^d
	#hprime = (h2$bw+ht$bw+hc$bw)/3
	hprime = (prod(h2$bw)+prod(ht$bw)+prod(hc$bw))/3
	#alphad = (c(h1$bw,h2$bw,ht$bw,hc$bw)/hprime)^d
	alphad = c(prod(h1$bw),prod(h2$bw),prod(ht$bw),prod(hc$bw))/hprime

	Delta2xeps = Delta2x + .00001

	## betaU
	mx0 = deprintize(npreg)(hc, exdat=data.frame(newx))
	mx1 = deprintize(npreg)(ht, exdat=data.frame(newx))
	
	Delta2x10 = Delta2x1
	
	Delta2x1 = mx1$mean - mx0$mean
	
	tauminus = Delta2x1 - sqrt(Delta2x)

	term1 = tauminus^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))/Delta2xeps
	term2 = sigma4lambda2/(4*alphad[2]*Delta2xeps*pixnn^2*(1-pixnn)^2)
	term3 = tauminus*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(Delta2xeps*pixnn*(1-pixnn))
	
	negterm3 = (term1+term2+term3 < 0)
	term3[negterm3] = 2*sqrt(term1*term2)[negterm3]
	
	seminus = sqrt( (term1+term2+term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans

	tauplus = Delta2x1 + sqrt(Delta2x)

	term1 = tauplus^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))/Delta2xeps
	term2 = sigma4lambda2/(4*alphad[2]*Delta2xeps*pixnn^2*(1-pixnn)^2)
	term3 = tauplus*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(Delta2xeps*pixnn*(1-pixnn))
	
	negterm3 = (term1+term2+term3 < 0)
	term3[negterm3] = 2*sqrt(term1*term2)[negterm3]

	seplus = sqrt( (term1+term2+term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans

	
	
	# se of Delta2x
	term1 = (mx1$mean - mx0$mean)^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))
	term2 = sigma4lambda2/(4*alphad[2]*pixnn^2*(1-pixnn)^2)
	term3 = (mx1$mean - mx0$mean)*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(pixnn*(1-pixnn))
	
	negterm3 = (term1+term2+term3 < 0)
	term3[negterm3] = 2*sqrt(term1*term2)[negterm3]

	seDelta2x = sqrt( 4*(term1 + term2 + term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans
	

	# se of betaU
	term1 = v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4])
	sebetaU = sqrt( term1*thetaK/fx )/sqrt(n*hprime) 
	
	
	## Create confidence intervals
	# For plus
	tauplusest = NULL
	tauminusest = NULL
	setauplus = NULL
	setauminus = NULL
	cutoffn = pmin(seDelta2x/2, 1)*max(3, (sqrt(n*hprime))^(max(1/3,d/6)))
	#cutoffn = (2/sqrt(n*hprime))*max(3, (sqrt(n*hprime))^(1/3))
	deltaused = 0
	for(i in 1:nrow(newx)){
		if(is.nan(seDelta2x[i]) | is.na(Delta2x[i])){
			setauplus = c(setauplus, sebetaU[i])
			tauplusest = c(tauplusest, mx1$mean[i] - mx0$mean[i])
			setauminus = c(setauminus, sebetaU[i])
			tauminusest = c(tauminusest, mx1$mean[i] - mx0$mean[i])
		} else if(Delta2x[i] < cutoffn[i]){
			setauplus = c(setauplus, sebetaU[i])
			tauplusest = c(tauplusest, mx1$mean[i] - mx0$mean[i])
			setauminus = c(setauminus, sebetaU[i])
			tauminusest = c(tauminusest, mx1$mean[i] - mx0$mean[i])
		} else {
			deltaused = deltaused + 1
			setauplus = c(setauplus, seplus[i])
			tauplusest = c(tauplusest, tauplus[i])
			setauminus = c(setauminus, seminus[i])
			tauminusest = c(tauminusest, tauminus[i])			
		}
	}
	
	tau = apply(newx, 1, taufn)
	
	# coverage
	mean((mx1$mean - mx0$mean - 2*sebetaU) < tau & tau < (mx1$mean - mx0$mean + 2*sebetaU))
	
	mean(((mx1$mean - mx0$mean - 2*sebetaU) < tau & tau < (mx1$mean - mx0$mean + 2*sebetaU))[Delta2x<cutoffn], na.rm=TRUE)
	
	mean( (((tauminusest - 2*setauminus) < tau & tau < (tauminusest + 2*setauminus)))[Delta2x>cutoffn], na.rm=TRUE)
	
	mean((((tauplusest - 2*setauplus) < tau & tau < (tauplusest + 2*setauplus)))[Delta2x>cutoffn], na.rm=TRUE)
	
	mean(((tauplusest - 2*setauplus) < tau & tau < (tauplusest + 2*setauplus)), na.rm=TRUE)
	
	summary(setauplus[Delta2x>cutoffn])
	
	mean(Delta2x>cutoffn)
	
	mean(((tauminusest - 2*setauminus) < tau & tau < (tauminusest + 2*setauminus))
	| ((tauplusest - 2*setauplus) < tau & tau < (tauplusest + 2*setauplus)), na.rm=TRUE)
	
	
	
	### Proposed method
	res_proposal = cbind(pixnn, tau, betaU = (mx1$mean - mx0$mean), Delta2x10, absDeltax = sqrt(Delta2x), cutoffn, nhd = n*hprime, seDelta2x, seplus, seminus, sebetaU, setauminus, setauplus, (tauminusest - 2*setauminus) < tau & tau < (tauminusest + 2*setauminus),
			(tauplusest - 2*setauplus) < tau & tau < (tauplusest + 2*setauplus))
				
	print(c('track =', 5))
	
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


	## Write results
	write.csv( cbind(n, d, res_proposal, res_cf, res_xrf),
		file = paste0('/blue/bkarmakar/bkarmakar/RPCOVA_ci/new_ll_fixed/simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v6_null.csv'), row.names = FALSE)

}




#install.packages('rgl')
#library(rgl)
#plot3d(x=x[,1], y=x[,2], z=tauminus) 

