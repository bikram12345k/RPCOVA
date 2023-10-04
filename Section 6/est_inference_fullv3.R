# args = commandArgs(trailingOnly=TRUE)
# id=as.numeric(args)[1]

# cigarTRUE = id==1

d1 = read.csv('compiled_data_2010.csv')
d1 = d1[,-1]



library(Matrix, lib.loc='/blue/bkarmakar/bkarmakar/struncatedP/newversion')

#library(Matrix)

#install.packages('ks')
library(ks)

######################
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



d2 = d1[which( (d1$age <= 43) & (d1$nprenvisits<=26 | d1$nprenvisits==99) & (d1$moncare<= 10 | is.na(d1$moncare)) ),]
d2$nprenvisits[d2$nprenvisits==99] = NA

x = d2[d2$treated | d2$control, c('age', 'nprenvisits', 'moncare')]
x[is.na(x$moncare), 'moncare'] = mean(x$moncare, na.rm=TRUE)
x[is.na(x$nprenvisits), 'nprenvisits'] = mean(x$nprenvisits, na.rm=TRUE)

z = d2[d2$treated | d2$control, 'treated']

y = d2[d2$treated | d2$control, 'birthweight']

y0 = d2[d2$resistant, 'birthweight']
x0 = d2[d2$resistant, c('age', 'nprenvisits', 'moncare')]
x0[is.na(x0$moncare), 'moncare'] = mean(x0$moncare, na.rm=TRUE)
x0[is.na(x0$nprenvisits), 'nprenvisits'] = mean(x0$nprenvisits, na.rm=TRUE)

n = nrow(x)#;
samp.1 = sample(n, floor(n*.4))#; if(id==3) samp.2 = sample(nrow(x0), 8000)
x = x[samp.1,] 
#cigarid=x$cigar==cigarTRUE; if(id==3) cigarid=rep(TRUE, nrow(x))
#x=x[cigarid,-5]
y = y[samp.1]#; y=y[cigarid]
z = z[samp.1]#; z=z[cigarid]

samp.2 = sample(nrow(x0), 20000)#; if(id==3) samp.2 = sample(nrow(x0), 12000)
x0 = x0[samp.2,]
#cigarid0 = x0$cigar==cigarTRUE; if(id==3) cigarid=rep(TRUE, nrow(x))
#x0=x0[cigarid0,-5]
y0 = y0[samp.2]#; y0=y0[cigarid0]

n = nrow(x); dim.d = ncol(x)
d = dim.d

#######################################
foo <- function(x, y) c(x, y)
deprintize<-function(f){
 return(function(...) {capture.output(w<-f(...));return(w);});
}

#######################################
xiseq = sort(unique(x[,1]))
for(i in 2:dim.d){
	if(i==2) newx = outer(xiseq,  sort(unique(x[,i])), FUN = Vectorize(foo, SIMPLIFY = FALSE))
	if(i > 2)
		newx = outer(newx, sort(unique(x[,i])), FUN = Vectorize(foo, SIMPLIFY = FALSE))
}
newx = t(sapply(newx, cbind))
dim(newx)

newx = newx[sort(sample(1:nrow(newx),5000)),]

#######################################
## Resistant Population Calibration Of VAriance (RPCOVA)

	yvec = as.numeric(y0)
	xdat = data.frame(x0)
	
	#### Estimation	####
	# Step 1
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat)
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat,bws=(h1$bw)^(1.1+(d/18)), bandwidth.compute = FALSE)

	mx <- deprintize(npreg)(h1)

	# Step 2
	rx = yvec - mx$mean

	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat)
	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat, bws=(h2$bw)^(1.1+(d/18)), bandwidth.compute = FALSE) #undersmooth
	sigma02x = deprintize(npreg)(h2, exdat=data.frame(newx))


	#sigma02 = 1+var(runif(3*n,0,2)^2)

#########################
	ok=0
	while(!ok){
	nnfit = deprintize(neuralnet::neuralnet)(z~., data=data.frame(z,x), hidden=8, linear.output=TRUE)
	ok=1
	if( !("weights" %in% names(nnfit)) )
		ok = 0
		
	}
## Estimate pi(x)
	pixnn = predict(nnfit, newdata=data.frame(newx))[,1]
	# Truncate away from 0 and 1
	pixnn = pmax(pixnn, .015)
	pixnn = pmin(pixnn, .985)
	#hist(pixnn)
#######################
	yvec = as.numeric(y)
	xdat = data.frame(x)
	
	#### Estimation	####
	# Step 1
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat)
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat,bws=(h1$bw)^(1.1+(d/18)), bandwidth.compute = FALSE)

	mx <- deprintize(npreg)(h1)

	# Step 2
	rx = yvec - mx$mean

	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat)
	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat, bws=(h2$bw)^(1.1+(d/18)), bandwidth.compute = FALSE) #undersmooth
	sigma2x = deprintize(npreg)(h2, exdat=data.frame(newx))

	# Step 3
	Delta2x2 = (sigma2x$mean - sigma02x$mean)/( pixnn*(1-pixnn) )

	hc = deprintize(npregbw)(ydat=yvec[z==0],xdat=xdat[z==0,], bws=h2$bw, bandwidth.compute = TRUE)
	hc = deprintize(npregbw)(ydat=yvec[z==0],xdat=xdat[z==0,], bws=hc$bw^(1.1+(d/18)), bandwidth.compute = FALSE)	#undersmooth

	ht = deprintize(npregbw)(ydat=yvec[z==1],xdat=xdat[z==1,], bws=h2$bw, bandwidth.compute = TRUE)
	ht = deprintize(npregbw)(ydat=yvec[z==1],xdat=xdat[z==1,], bws=ht$bw^(1.1+(d/18)), bandwidth.compute = FALSE)	#undersmooth
	
	# run constraint optimization
	Delta2x1 = sapply(1:nrow(newx), function(id) {			
			if(id %% 500 == 0) {
				print(id)
				deprintize(gc)()}
			deprintize(const_loclin2_d)(newx[id,], h=ht$bw, hs=hc$bw, as.matrix(x), d, const=  Delta2x2[id], y)
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

negterm3 = which(term1+term2+term3 < 0)
term3[negterm3] = 2*sqrt(term1*term2)[negterm3]

seminus = sqrt( (term1+term2+term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans

tauplus = Delta2x1 + sqrt(Delta2x)

term1 = tauplus^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))/Delta2xeps
term2 = sigma4lambda2/(4*alphad[2]*Delta2xeps*pixnn^2*(1-pixnn)^2)
term3 = tauplus*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(Delta2xeps*pixnn*(1-pixnn))

negterm3 = which(term1+term2+term3 < 0)
term3[negterm3] = 2*sqrt(term1*term2)[negterm3]

seplus = sqrt( (term1+term2+term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans



# se of Delta2x
term1 = (mx1$mean - mx0$mean)^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))
term2 = sigma4lambda2/(4*alphad[2]*pixnn^2*(1-pixnn)^2)
term3 = (mx1$mean - mx0$mean)*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(pixnn*(1-pixnn))

negterm3 = which(term1+term2+term3 < 0)
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

save(list = ls(all.names = TRUE), file =  paste0("Est_INFERENCE_032323_fullv3.RData"),  envir = .GlobalEnv)
