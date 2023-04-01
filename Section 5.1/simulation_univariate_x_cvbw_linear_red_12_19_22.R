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
n=3000
d = 1 # dimension of x
u = runif(n,0,2)	# unobserved confounder
x = rmvnorm(n, rep(3,d))

z = rbinom(n,1, 1/(1+exp(-(.5*x*u+.5*x-1.5))))
y0 = 1-u*3+x+rnorm(n,0,.5)
taufn <- function(x) .5*x
tau = taufn(x)

y = y0 + tau*z
nnfit = neuralnet::neuralnet(z~x, data=data.frame(z,x), hidden=4)
ok=1
if( !("weights" %in% names(nnfit)) )
	ok = 0
}

########################################################

## Resistant POpulation VAriance (RPOVA)
sigma02 = .25+4*9*(1/12)

## Estimate pi(x)
# there is some issue with convergence sometimes(!?)
newx = seq(1,5,.025)

yvec = as.numeric(y)
xdat = data.frame(x)
h1 <- npregbw(ydat=yvec,xdat=xdat,regtype='ll')
mx <- npreg(h1)



rx = yvec - mx$mean

#plot(x, rx^2)




h2 <- npregbw(ydat=rx^2,xdat=xdat,regtype='ll')
sigma2x = npreg(h2, exdat=data.frame(newx))

pixnn = predict(nnfit, newdata=data.frame(sigma2x$eval))[,1]

Delta2x2 = (sigma2x$mean - sigma02)/( pixnn*(1-pixnn) )

hc = npregbw(ydat=yvec[z==0],xdat=xdat[z==0,],regtype='ll')$bw
ht = npregbw(ydat=yvec[z==1],xdat=xdat[z==1,],regtype='ll')$bw


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

save.image(file='linear_red_12_19_22.RData')


dev.set(2)
plot(x=res[1,1:161], y = res[1,1*162+(1:161)], type='l', xlim=c(1,5),
			ylim=c(-5,5), col='gray', lwd=1.5, ylab='')

for(itr in 2:100){
par(new=TRUE)
plot(x=res[itr,1:161], y = res[itr,1*162+(1:161)], type='l', xlim=c(1,5),
			ylim=c(-5,5), col='gray', lwd=1.5, ylab='')
}
par(new=TRUE)
for(itr in 1:100){
par(new=TRUE)
plot(x=res[itr,1:161], y = res[itr,2*162+(1:161)], type='l', xlim=c(1,5),
			ylim=c(-5,5), col='gray50', lwd=1.5, ylab='')
}
par(new=TRUE)
plot(x=res[1,1:161], y = taufn(res[1,1:161]), type='l', xlim=c(1,5),
			ylim=c(-5,5),  lwd=1.5, col='black', ylab='')


########################################################
## Plot
dev.set(2)
plot(x=newx, y = Delta2x1 - sqrt(Delta2x), type='l', xlim=c(1,5),
			ylim=c(-5,5), col='blue', lwd=1.5, ylab='')
par(new=TRUE)
plot(x=newx, y = Delta2x1 + sqrt(Delta2x), type='l', xlim=c(1,5),
			ylim=c(-5,5), col='red', lwd=1.5, ylab='')
par(new=TRUE)
plot(x=newx, y = Delta2x1, type='l', xlim=c(1,5),
			ylim=c(-5,5),  lwd=1.5, ylab='')
par(new=TRUE)
plot(x=newx, y = taufn(newx), type='l', xlim=c(1,5),
			ylim=c(-5,5),  lwd=1.5, col='gray', ylab='')
legend(x=1,y=-3.1,legend=c('True', 'Naive', 'Est1', 'Est2'), col=c('gray', 'black', 'red', 'blue'), lty=1, lwd=2)

par(new=TRUE)
# causal forest
plot(x=newx, y = c.pred$predictions, type='l', xlim=c(1,5),
			ylim=c(-5,5), col='green', lwd=1.5, ylab='')

par(new=TRUE)
# bcf
plot(x=newx, y = tauhatnewx , type='l', xlim=c(1,5),
			ylim=c(-5,5), col='purple', lwd=1.5, ylab='')

par(new=TRUE)
# x-learner rf
plot(x=newx, y = cate_esti_xrf, type='l', xlim=c(1,5),
			ylim=c(-5,5), col='orange', lwd=1.5, ylab='')

