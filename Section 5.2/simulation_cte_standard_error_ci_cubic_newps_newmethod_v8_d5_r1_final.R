## Simulation for Section 5.2 the proposed method for dimension d = 5
args = commandArgs(trailingOnly=TRUE)
id=as.numeric(args)[1]

dseq = 5 
seedseq = c(8041, 8042, 8043, 8044, 8045)
#seedseq = as.numeric(sapply((1:9), function(x) x*c(432, 140, 645, 236, 513, 142, 2635, 3616, 6231, 6216, 321, 408, 451))) #
nseq = c(2000, 4000, 6000, 8000, 10000, 12000)
cscaleseq1 = c(1,1.5,2,2.5,3, 3.5, 4)
cscaleseq2 = c(.1, .25, .5, 4.5, 5)
cscaleseq3 = c(.02, .05, .1, .25, .5, .75, 10, 25)/10	
cscaleseq = sort(c(cscaleseq1, cscaleseq2, cscaleseq3))	# Varianble bandwidth, see line 112

foo <- function(x, y) c(x, y)
setup = outer(dseq, seedseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup, nseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup, cscaleseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)


library(Matrix, lib.loc='/blue/bkarmakar/bkarmakar/struncatedP/')

#library(Matrix)

#install.packages('ks')
library(ks)

######################
n = setup[3,id]; d = setup[1,id]; thisseed = setup[2,id]; cscale = setup[4,id]

set.seed(thisseed)
print(c(n, d, thisseed, cscale))

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


## Estimation of sigma02
	
		u0 = runif(floor(n^(1.1)),0,2)	# unobserved confounder
		x0 = rmvnorm(floor(n^(1.1)), rep(3,d)) #3 + matrix(runif(floor(n^(1.1))*d, -1, 1), nrow=floor(n^(1.1)), ncol=d)
		#rmvnorm(floor(n^(1.1)), rep(3,d))
		y00 = 1+u0+rowSums(x0)+rnorm(floor(n^(1.1)),0,.5)*(x0[,1]-3)+rnorm(floor(n^(1.1)),0,1)
		mx0 = deprintize(gbart)(data.frame(x0), y00)
		mx02 = deprintize(gbart)(data.frame(x0), y00^2)
		
		sigma0fit = deprintize(gbart)(data.frame(x0[,1]), colMeans(mx02$yhat.train) - colMeans(mx0$yhat.train)^2) 

sigma0fit1 = deprintize(gbart)(data.frame(x0[,1]), (y00 - colMeans(mx0$yhat.train))^2) 

	
	##########################
	# Calculate the evaluation points
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
	
	set.seed(0841)	# Set the points same across iterations
	newx = newx[seq(1,nrow(newx),4),]
	newx = newx + matrix(rnorm(nrow(newx)*d, 0,.2), ncol=d)

	sigma02tr = var(runif(n*2,0,2))/4+.25*(newx[,1]-3)^2+1

	newx = newx[sigma02tr>1.3,]
	set.seed(thisseed)
	#dim(newx)
	sigma02 = (colMeans(deprintize(predict)(sigma0fit, data.frame(newx[,1]))) + 
		colMeans(deprintize(predict)(sigma0fit1, data.frame(newx[,1]))))/2

## Local linear regression model setup
npregbw_my <- function(ydat,xdat,regtype='ll', bws,
			ftol=1.5e-01, itmax=2000, tol=1e-1, ckertype='gaussian',bandwidth.compute = TRUE,...){
			
	bws = rep(cscale*n^(-2/(2*d+1)), d)
	res = deprintize(npregbw)(ydat=ydat,xdat=xdat,regtype=regtype, bws=bws,
				 itmax=itmax, ckertype=ckertype,bandwidth.compute = FALSE,...)
	res
}
	


#for(itr in 1:40){
for(itr in 1:ifelse(thisseed==8041,100,20)){
	print(itr)

	ok=0
	while(!ok){
		# unobserved confounder
		u = runif(n,0,2)
		x = rmvnorm(n, rep(3,d)) #cpvariates

		#Propensity score model
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
	#Sys.time()
	#pixnn = predict(nnfit, newdata=data.frame(newx))[,1]
	
	print(c('track =', 1))
	
	## Resistant Population Calibration Of VAriance (RPCOVA)
	# sigma02 = 1+var(runif(3*n,0,2)^2)/25+(x[1]-3)^2/4

	
	## Estimate pi(x)
	pixnn = colMeans(deprintize(predict)(psfit, data.frame(newx))$prob.test)#[,1]
	#pstr = apply(newx, 1, function(x) mean(1/(1+exp(-(.5*mean(x)+2*runif(1000,0,2)-3)))))
	#plot(pstr, pixnn)
	#abline(a=0,b=1)
	
		
	# Truncate away from 0 and 1
	pixnn = pmax(pixnn, .015)
	pixnn = pmin(pixnn, .985)
	#hist(pixnn)
	
	yvec = as.numeric(y)
	xdat = data.frame(x)
	
	## Bandwidth under smoothing
	us <- function(bw){ pw = 1.1*(bw<1)+.9*(bw>1); bw^pw }
	
		
	kern = 'gaussian'
	
	#### Estimation	####
	# Step 1: estimate m(x) = E(Y|x) and E(Y^2|x)
	h1 <- deprintize(npregbw_my)(ydat=yvec,xdat=xdat,regtype='ll', bws=rep(.1,d),
			 itmax=2000, ckertype=kern)
	cat('h1:', h1$bw)
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat,bws=us(h1$bw), bandwidth.compute = FALSE,regtype='ll', ckertype=kern)

	mx <- deprintize(npreg)(h1)

	h12 <- deprintize(npregbw)(ydat=yvec^2,xdat=xdat,bws=us(h1$bw), bandwidth.compute = FALSE,regtype='ll', ckertype=kern)
	mx2 <- deprintize(npreg)(h12)
	
	# Step 2: Estimation of Var(Y|x) using both E( (y-m(x))^2 | x) and E(y^2|x) - m(x)^2
	rx = yvec - mx$mean

	h2 <- deprintize(npregbw_my)(ydat=rx^2,xdat=xdat,regtype='ll', bws=rep(.1,d), 
				itmax=2000, ckertype=kern, remin=FALSE)
	cat('\nh2:', h2$bw)

	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat, bws=us(h2$bw), bandwidth.compute = FALSE,regtype='ll', ckertype=kern) #undersmooth
	sigma2x = deprintize(npreg)(h2, exdat=data.frame(newx))
	
	mx <- deprintize(npreg)(h1, exdat=data.frame(newx))
	mx2 <- deprintize(npreg)(h12, exdat=data.frame(newx))
	sigmatemp = pmax(mx2$mean - mx$mean^2, sigma2x$mean)
	
	sigma2x$mean[sigmatemp>0] = sigmatemp[sigmatemp>0]
	
	# Step 3: Estimation of Delta^2(x) by Mixed integer programming
	Delta2x2 = (sigma2x$mean - sigma02)/( pixnn*(1-pixnn) )

	# ll model setup for E(Y | x, Z=0)
	hc = deprintize(npregbw_my)(ydat=yvec[z==0],xdat=xdat[z==0,], bws=h1$bw, bandwidth.compute = TRUE,regtype='ll', itmax=2000, ckertype=kern)
	cat('\nhc:', hc$bw)
	
	hc = deprintize(npregbw)(ydat=yvec[z==0],xdat=xdat[z==0,], bws=us(hc$bw), bandwidth.compute = FALSE,regtype='ll', ckertype=kern)	#undersmooth

	# ll model setup for E(Y | x, Z=1)
	ht = deprintize(npregbw_my)(ydat=yvec[z==1],xdat=xdat[z==1,], bws=h1$bw, bandwidth.compute = TRUE,regtype='ll', itmax=2000, remin=FALSE, ckertype=kern)
	cat('\nht:', ht$bw)
	
	ht = deprintize(npregbw)(ydat=yvec[z==1],xdat=xdat[z==1,], bws=us(ht$bw), bandwidth.compute = FALSE,regtype='ll', ckertype=kern)	#undersmooth
	
	# run constraint optimization
	Delta2x1 = sapply(1:nrow(newx), function(id) {			
			if(id %% 500 == 0) {
				print(id)
				deprintize(gc)()}
			deprintize(const_loclin2_d)(newx[id,], h=ht$bw, hs=hc$bw, x, d, const=  Delta2x2[id], y)
			})
	deprintize(gc)()
	print(c('track =', 2))	

	# Step 4: Final estimates
	if(is.list(Delta2x1))
		Delta2x1 = sapply(Delta2x1, function(x) ifelse(length(x)==0, NA, x))
	Delta2x = pmax(Delta2x1^2 - Delta2x2, 0)
	print(c('track =', 3))
	
	mx0 = deprintize(npreg)(hc, exdat=data.frame(newx))
	mx1 = deprintize(npreg)(ht, exdat=data.frame(newx))
	Delta2x = pmin((mx1$mean-mx0$mean)^2 - Delta2x2, Delta2x)
	Delta2x = pmax(Delta2x, 0)

	tauminus = Delta2x1 - sqrt(Delta2x)
	tau = apply(newx, 1, taufn)
	plot(tau, tauminus); abline(a=0,b=1)



	#### Calculating terms for the standard errors	####
	
	mx0 = deprintize(npreg)(hc, exdat=data.frame(x[z==0,]))
	mx1 = deprintize(npreg)(ht, exdat=data.frame(x[z==1,]))

	rx0 = yvec[z==0] - mx0$mean
	rx1 = yvec[z==1] - mx1$mean

	# Using bart, estimate the relavant quantities. See remark 6 
	v02fit = deprintize(gbart)(data.frame(xdat[z==0,]), rx0^2) 
		v02 = colMeans(deprintize(predict)(v02fit, data.frame(newx)))
	v12fit = deprintize(gbart)(data.frame(xdat[z==1,]), rx1^2) 
		v12 = colMeans(deprintize(predict)(v12fit, data.frame(newx)))	
	sigma4lambda2fit = deprintize(gbart)(data.frame(xdat), rx^4) 
		sigma4lambda2 = colMeans(deprintize(predict)(sigma4lambda2fit, data.frame(newx)))	

	sigma2nu0eta0fit = deprintize(gbart)(data.frame(xdat[z==0,]), rx[z==0]^2*rx0) 
		sigma2nu0eta0 = colMeans(deprintize(predict)(sigma2nu0eta0fit, data.frame(newx)))	
	sigma2nu1eta1fit = deprintize(gbart)(data.frame(xdat[z==1,]), rx[z==1]^2*rx1) 
		sigma2nu1eta1 = colMeans(deprintize(predict)(sigma2nu1eta1fit, data.frame(newx)))	
			
	v02t = v02; v12t = v12
	v02  =list();v02$mean = v02t
	v12  =list();v12$mean = v12t

	v02$mean = pmax(0, v02$mean)
	v02$mean[v02$mean==0] = median(v02$mean[v02$mean!=0])
	v12$mean = pmax(0, v12$mean)
	v12$mean[v12$mean==0] = median(v12$mean[v12$mean!=0])
	sigma4lambda2 = pmax(0, sigma4lambda2)
	sigma4lambda2[sigma4lambda2==0] = median(sigma4lambda2[sigma4lambda2!=0])

	print(c('track =', 4))
	
	## Estimate density	
	# By gaussian mixture model
	fx = deprintize(densityMclust)(x, plot=FALSE)
	fx = predict(fx, newdata=newx)
	#fx = pmax(fx, .0001)

	
	#### Variance calculations  #### variance formulas (6), (7)
	
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
	Delta2x1 = Delta2x10	
	
	hprime = min(c(.8^d, hprime))

	tauminus = Delta2x1 - sqrt(Delta2x)

	term1 = tauminus^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))/Delta2xeps	
	term2 = sigma4lambda2/(4*alphad[2]*Delta2xeps*pixnn^2*(1-pixnn)^2)
	term3 = tauminus*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(Delta2xeps*pixnn*(1-pixnn))
	
	term3[is.na(term3)] = 2*sqrt(pmax(0,term1*term2))[is.na(term3)]

	negterm3 = (term1+term2+term3 < 0)
	term3[negterm3] = 2*sqrt(term1*term2)[negterm3]
	
	seminus = sqrt( (term1+term2+term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans

	tauplus = Delta2x1 + sqrt(Delta2x)

	term1 = tauplus^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))/Delta2xeps
	term2 = sigma4lambda2/(4*alphad[2]*Delta2xeps*pixnn^2*(1-pixnn)^2)
	term3 = tauplus*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(Delta2xeps*pixnn*(1-pixnn))

	term3[is.na(term3)] = 2*sqrt(pmax(0,term1*term2))[is.na(term3)]

	negterm3 = (term1+term2+term3 < 0)
	term3[negterm3] = 2*sqrt(term1*term2)[negterm3]

	seplus = sqrt( (term1+term2+term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans

	
	
	# se of Delta2x
	term1 = (mx1$mean - mx0$mean)^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))
	term2 = sigma4lambda2/(4*alphad[2]*pixnn^2*(1-pixnn)^2)
	term3 = (mx1$mean - mx0$mean)*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(pixnn*(1-pixnn))
	
	term3[is.na(term3)] = 2*sqrt(pmax(0,term1*term2))[is.na(term3)]
	
	negterm3 = (term1+term2+term3 < 0)
	term3[negterm3] = 2*sqrt(term1*term2)[negterm3]

	seDelta2x = sqrt( 4*(term1 + term2 + term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans
	

	# se of betaU
	term1 = v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4])
	sebetaU = sqrt( term1*thetaK/fx )/sqrt(n*hprime) 
	
	plot(tau[Delta2x>.1], (tauminus+2*seminus)[Delta2x>.1]); abline(a=0,b=1)
	#mean(tau < tauminus + 2*seminus, na.rm=TRUE)


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
	print(mean(tau<tauminus+2*setauminus, na.rm=TRUE))
	print(mean(tau<tauminus+2*seminus, na.rm=TRUE))
	# coverage
	print(mean((mx1$mean - mx0$mean - 2*sebetaU) < tau & tau < (mx1$mean - mx0$mean + 2*sebetaU)))
	
	print(mean(((mx1$mean - mx0$mean - 2*sebetaU) < tau & tau < (mx1$mean - mx0$mean + 2*sebetaU))[Delta2x<cutoffn], na.rm=TRUE))
	
	print(mean( (((tauminusest - 2*setauminus) < tau & tau < (tauminusest + 2*setauminus)))[Delta2x>cutoffn], na.rm=TRUE))
	
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
	
	
	# ## Results file
	thisrun = 0
	if(cscale %in% cscaleseq1) cscalef = cscale*10
	if(cscale %in% cscaleseq2) cscalef = cscale*100
	if(cscale %in% cscaleseq3) cscalef = cscale*1000+99
	if(thisseed == 8041) cscalef = cscalef
	if(thisseed != 8041) cscalef = cscale*1000+99
	
	fname = paste0('/blue/bkarmakar/bkarmakar/RPCOVA_ci/new_ll_fixed/simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw',cscalef,'.csv')
	while(file.exists(fname)){
		thisrun = thisrun + 1
		fname = paste0('/blue/bkarmakar/bkarmakar/RPCOVA_ci/new_ll_fixed/simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw',cscalef,'.csv')
	}
	## Write results
	write.csv( cbind(n, d, res_proposal), #res_cf, res_xrf),
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

