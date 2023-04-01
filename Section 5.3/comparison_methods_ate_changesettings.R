args = commandArgs(trailingOnly=TRUE)
id=as.numeric(args)[1]

#idseq = (1:240)[-c(36:40, 61:80, 91:120, 126:240)]
#idseq = (1:240)[-c(36:40, 61:80, 91:120, 126:240)]
#id = idseq[id]

seedseq <- c(1111, 6640, 7042, 8674, 4168, 9448, 2666,  515,  839, 8155, 3096, 7441, 2976, 6240, 1624,
		1237, 4384, 4191, 2688, 6721, 4241, 1921, 5693, 1218, 8216)
#c(1237, 4384, 4191, 2688, 6721, 4241, 1921, 5693, 1218, 8216) #c(1237, 4384, 4191, 2688, 6721) #, 4241, 1921, 5693, 1218, 8216) 
dseq = 5:8
nseq = c(500, 1000, 1500, 2000)#c(500, 1000, 1500, 2000, 3000, 4000)

foo <- function(x, y) c(x, y)
setup = outer(seedseq , dseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup , nseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)	#800

thisseed = setup[1,id]
set.seed(thisseed)
d = setup[2,id]
n = setup[3,id]




#install.packages("DoubleML")
#install.packages('drtmle')
#install.packages('bartMachine')
#install.packages('xgboost')

library(mvtnorm)
#library(DoubleML)
library(drtmle)
library(SuperLearner)
library(xgboost)
#library(bartMachine)


res = NULL

#for(d in 3:5){
for(itr in 1:25){
#if(itr<=500) n = 1500
#if(itr>500 & itr <=1000) n = 2000
#if(itr>1000 & itr <=1500) n = 3000



		u = runif(n,0,2)	# unobserved confounder
		x = rmvnorm(n, rep(3,d))

		#z = rbinom(n,1, 1/(1+exp(-(.4*u+.5*rowMeans(x)-1.5))))
		z = rbinom(n,1, 1/(1+exp(-(rowMeans(x)*u+.5*rowMeans(x)-1.5))))
		y0 = 1+u^2+rowMeans(x)+rnorm(n,0,1)
		if(d==1)
			taufn <- function(x) .5*x[1]
		if(d>1)
			taufn <- function(x) .5*x[1]+mean(x[-1])/3
		tau = apply(x, 1, taufn)
		#tau = rowSums(taufn(x))

		y = y0 + tau*z



dat <- data.frame(y=y,z=z,x=x)
#colnames(d)

#lm(y~.,data=d)

res_ols = (mean(tau) - lm(y~.,data=dat)$coef[2])^2


fit1 <- drtmle(
	W = data.frame(x), A = z, Y = y, a_0 = c(1, 0),
family = gaussian(),
stratify = FALSE,
SL_Q = c("SL.glm", "SL.xgboost"),
SL_g = c("SL.glm", "SL.xgboost"),
SL_Qr =	 c("SL.glm", "SL.xgboost"),
SL_gr = c("SL.glm", "SL.xgboost"), maxIter = 1
)

res_dr = (mean(tau) - diff(rev(fit1$drtmle$est)))^2

res = rbind(res, c(n,d,res_ols, res_dr))
print(c(itr, n, d, res_ols, res_dr))
}
#}

#save(list = ls(all.names = TRUE), file = "comparison_results_ate_926.RData", envir = .GlobalEnv)