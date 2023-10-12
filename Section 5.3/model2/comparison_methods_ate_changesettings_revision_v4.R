#############
######## Simulation for Average treatment effect using linear regression and double robustness
#############

args = commandArgs(trailingOnly=TRUE)
id=as.numeric(args)[1]

#idseq = (1:240)[-c(36:40, 61:80, 91:120, 126:240)]
#idseq = (1:240)#[-c(36:40, 61:80, 91:120, 126:240)]
#idseq = c(131:135, 140:240)
#id = idseq[id]

seedseq <- 2*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721) #, 4241, 1921, 5693, 1218, 8216) 
seedseq <- c(seedseq, 3*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))
seedseq <- c(seedseq, 4*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))
seedseq <- c(seedseq, 5*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))

dseq = 2:5
nseq = c(500, 1000, 1500, 2000)

foo <- function(x, y) c(x, y)
setup = outer(seedseq , dseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup , nseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)	#800

thisseed = setup[1,id]
set.seed(thisseed)
d = setup[2,id]
n = setup[3,id]




library(mvtnorm)
#library(DoubleML)
library(drtmle)
library(SuperLearner)
library(xgboost)
#library(bartMachine)

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
## Calculate true ATE and ATT
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
		
		ate = mean(y1 -y0)
		att = mean( (y1-y0)[z==1])

rm(x, u, z, y1, y0, tau)
####################################



res = NULL

#for(d in 3:5){
for(itr in 1:10){
	
	# Data generating model

		u = runif(n,0,2)	# unobserved confounder
		x = rmvnorm(n, rep(3,d))

		z = rbinom(n,1, 1/(1+exp(-(u/2+.5*rowMeans(x)-1.5))))
		y0 = 1+u^2+d*rowMeans(x)+rnorm(n,0,1)*rowMeans(x-3)/6
		if(d==1)
			taufn <- function(x) .5*x[1]
		if(d>1)
			taufn <- function(x) .5*x[1]+mean(x[-1])/3
		tau = apply(x, 1, taufn)

		y = rep(0, n)
		
		y[z==1] = y0[z==1] + (tau[z==1]+rnorm(sum(z),0,.25))
		y[z==0] = y0[z==0] + rnorm(sum(1-z),0,.25)
		




dat <- data.frame(y=y,z=z,x=x)

## Linear regression
res_ols_att = (att - lm(y~.,data=dat)$coef[2])^2
res_ols_ate = (ate - lm(y~.,data=dat)$coef[2])^2


## Doubly robust method
fit1 <- drtmle(
	W = data.frame(x), A = z, Y = y, a_0 = c(1, 0),
family = gaussian(),
stratify = FALSE,
SL_Q = c("SL.glm", "SL.xgboost"),
SL_g = c("SL.glm", "SL.xgboost"),
SL_Qr =	 c("SL.glm", "SL.xgboost"),
SL_gr = c("SL.glm", "SL.xgboost"), maxIter = 1
)

res_dr_att = (att - diff(rev(fit1$drtmle$est)))^2
res_dr_ate = (ate - diff(rev(fit1$drtmle$est)))^2

res = rbind(res, c(n,d,res_ols_att, res_ols_ate, res_dr_att, res_dr_ate))
print(c(itr, n, d, res_ols_att, res_ols_ate, res_dr_att, res_dr_ate))
}
#}
