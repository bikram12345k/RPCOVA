rm(list=ls(all=TRUE))


library(mvtnorm)
library(KernSmooth)
library(neuralnet)
library(gurobi)	# See directions for installation in the document
library(np)

library(grf)
library(bcf)
library(causalToolbox)

source('supporting_fn.R')

library(foreach)
library(doSNOW)
library(parallel)


cl = makeCluster(5)
registerDoSNOW(cl)


Itr = 100
res = foreach(itr = 1:Itr, .combine='rbind', 
	.packages=c('gurobi', 'np', 'KernSmooth', 'mvtnorm', 'neuralnet',
			'grf','bcf','causalToolbox')) %dopar% {

set.seed(10+itr*5)

source('supporting_fn.R')

########################################################
## Generate data
ok=0
while(!ok){
n=2000
d = 1 # dimension of x
u = runif(n,0,2)	# unobserved confounder
x = rmvnorm(n, rep(0,d))

z = rbinom(n,1, 1/(1+exp(-(.5*u+u^2+.5*x))))
y0 = 1+2*u*sqrt(abs(x+4))+x+rnorm(n,0,.5)
taufn <- function(x) .5*abs(x)^3
tau = taufn(x)

y = y0 + tau*z

nnfit = try(neuralnet::neuralnet(z~x, data=data.frame(z,x), hidden=8))
ok=1
if( !("weights" %in% names(nnfit)) )
	ok = 0
}

########################################################

## Resistant POpulation VAriance (RPOVA)
sigma02 = function(x) .25+4*4*(1/12)*abs(x+4)

## Estimate pi(x)
# there is some issue with convergence sometimes(!?)
newx = seq(-2,2,.025)

yvec = as.numeric(y)
xdat = data.frame(x)
h1 <- npregbw(ydat=yvec,xdat=xdat,regtype='ll')
mx <- npreg(h1)



rx = yvec - mx$mean

#plot(x, rx^2)




h2 <- npregbw(ydat=rx^2,xdat=xdat,regtype='ll')
sigma2x = npreg(h2, exdat=data.frame(newx))

pixnn = predict(nnfit, newdata=data.frame(sigma2x$eval))[,1]

Delta2x2 = (sigma2x$mean - sigma02(sigma2x$eval))/( pixnn*(1-pixnn) )
Delta2x2 = unlist(Delta2x2)

hc = npregbw(ydat=yvec[z==0],xdat=xdat[z==0,])$bw
ht = npregbw(ydat=yvec[z==1],xdat=xdat[z==1,])$bw


# run constraint optimization
Delta2x1 = sapply(1:length(newx), function(id) const_loclin2(newx[id], h=ht, hs=hc, x, d, const=Delta2x2[id]))
if(is.list(Delta2x1))
	Delta2x1 = sapply(Delta2x1, function(x) ifelse(length(x)==0, NA, x))
Delta2x = pmax(Delta2x1^2 - Delta2x2, 0)

###


c.forest <- causal_forest(xdat, yvec, z, num.trees = 4000)
c.pred <- predict(c.forest, data.frame(newx), estimate.variance = FALSE)


pixnn = predict(nnfit, newdata=xdat)[,1]
bcf_fit = bcf(yvec, z, as.matrix(xdat), as.matrix(xdat), 
		pixnn, nburn=1000, nsim=1000)

# Get posterior of treatment effects
tau_post = bcf_fit$tau
tauhat = colMeans(tau_post)

mapnewx = sapply(newx, findInterval, vec=sort(x))
tauhatnewx = (tauhat[order(as.numeric(x))])[mapnewx]

xl_rf <- X_RF(feat = data.frame(x=x, x2=x^2), tr = z, yobs = yvec)
cate_esti_xrf <- EstimateCate(xl_rf, data.frame(x=newx, x2=newx^2))



c(newx,NA,
Delta2x1 - sqrt(Delta2x), NA, 
Delta2x1 + sqrt(Delta2x), NA, 
Delta2x1, NA,
c.pred$predictions,NA,
tauhatnewx, NA,
cate_esti_xrf)

}

save.image(file='expcubic_varysigma0_12_19_22.RData')
