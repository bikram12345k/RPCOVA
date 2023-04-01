args = commandArgs(trailingOnly=TRUE)
id=as.numeric(args)[1]

#idseq = (1:240)[-c(36:40, 61:80, 91:120, 126:240)]
idseq = (1:240)#[-c(36:40, 61:80, 91:120, 126:240)]
id = idseq[id]

seedseq <- c(4241, 1921, 5693, 1218, 8216) #c(1237, 4384, 4191, 2688, 6721) #, 4241, 1921, 5693, 1218, 8216) 
dseq = 1:8
nseq = c(500, 600, 1000, 1500, 2000, 3000)

foo <- function(x, y) c(x, y)
setup = outer(seedseq , dseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup , nseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)

thisseed = setup[1,id]
set.seed(thisseed)
d = setup[2,id]
n = setup[3,id]

print(c(n,d,thisseed))

Itr = 100 #60

library(Matrix)
library(mvtnorm)
library(KernSmooth)
library(neuralnet)
library(slam)
library(gurobi)	# See directions for installation in the document
library(np)

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


res = NULL
res_mse = NULL

for(itr in 1:Itr){


	ok=0
	while(!ok){
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
		nnfit = deprintize(neuralnet::neuralnet)(z~., data=data.frame(z,x), hidden=8)
		ok=1
		if( !("weights" %in% names(nnfit)) )
			ok = 0
	}

	pixnn = predict(nnfit, newdata=x)[,1]
	#hist(pixnn)

	## Resistant POpulation VAriance (RPOVA)
	sigma02 = 1+16*(4/45)

	## Estimate pi(x)
	# there is some issue with convergence sometimes(!?)
	newx = x # rmvnorm(2*n, rep(3,d)) #seq(1,5,.025)

	yvec = as.numeric(y)
	xdat = data.frame(x)
	h1 <- deprintize(npregbw)(ydat=yvec,xdat=xdat,regtype='ll')
	mx <- deprintize(npreg)(h1)


	rx = yvec - mx$mean

	#plot(x, rx^2)

	h2 <- deprintize(npregbw)(ydat=rx^2,xdat=xdat, bws=h1$bw,regtype='ll',bandwidth.compute = FALSE)
	sigma2x = deprintize(npreg)(h2, exdat=data.frame(newx))


	Delta2x2 = (sigma2x$mean - sigma02)/( pixnn*(1-pixnn) )

	hc = deprintize(npregbw)(ydat=yvec[z==0],xdat=xdat[z==0,],regtype='ll')#$bw
	ht = deprintize(npregbw)(ydat=yvec[z==1],xdat=xdat[z==1,],regtype='ll')#$bw

	mx0 = deprintize(npreg)(hc, exdat=data.frame(newx))
	mx1 = deprintize(npreg)(ht, exdat=data.frame(newx))


	# run constraint optimization
	Delta2x1 = sapply(1:nrow(newx), function(id) 
			deprintize(const_loclin2_d)(newx[id,], h=ht$bw, hs=hc$bw, x, d, const=Delta2x2[id], y))

	if(is.list(Delta2x1))
		Delta2x1 = sapply(Delta2x1, function(x) ifelse(length(x)==0, NA, x))
	Delta2x = pmax(Delta2x1^2 - Delta2x2, 0)


	
	res_mse_temp = c( (mean(rm_outliers(Delta2x1 - sqrt(Delta2x) - tau), na.rm = TRUE))^2 ,
		(mean(rm_outliers(Delta2x1 + sqrt(Delta2x) - tau), na.rm = TRUE))^2,
		(mean(rm_outliers(Delta2x1 - tau), na.rm = TRUE))^2 )
	
	res_mse_temp
	
	restemp = rbind(#newx,NA,
		c(tau, NA, NA),
		c(Delta2x1 - sqrt(Delta2x), NA, res_mse_temp[1]),
		c(Delta2x1 + sqrt(Delta2x), NA, res_mse_temp[2]),
		c(Delta2x1, NA, res_mse_temp[3]) )
		
	res = rbind(res, c(#newx,NA,
		Delta2x1 - sqrt(Delta2x), NA, 
		Delta2x1 + sqrt(Delta2x), NA, 
		Delta2x1, NA, NA, res_mse_temp ))
	
	
	res_mse = rbind(res_mse, res_mse_temp)
	
	print(c(itr, res_mse_temp ))
	
	#write.csv(cbind(thisseed,d,n,res_mse), file=paste0('linear_stresstest_d_mse.csv'), append=TRUE, row.names = FALSE, col.names = FALSE)
	
	write.csv(cbind(thisseed,d,n,restemp), file=paste0('linear_stresstest_d_',n,'_',d,'_',thisseed,'_',itr,'.csv'), row.names = FALSE, col.names = FALSE)

}

#write.csv(cbind(thisseed,d,n,res), file=paste0('linear_stresstest_d_',n,'_',d,'_',thisseed,'.csv'), row.names = FALSE, col.names = FALSE)


#write.csv(cbind(thisseed,d,n,res), file=paste0('linear_stresstest_d_',n,'_',d,'_',thisseed,'.csv'), row.names = FALSE, col.names = FALSE)
print('final')
print(cbind(thisseed,d,n,res_mse))
write.csv(cbind(thisseed,d,n,res), file=paste0('linear_stresstest_d_',n,'_',d,'_',thisseed,'.csv'), append=TRUE, row.names = FALSE, col.names = FALSE)
