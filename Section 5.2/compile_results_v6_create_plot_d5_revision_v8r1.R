## Compile results for the proposed method

dseq = 5 #c(5,5,5,5,5,5,5,5,5)
seedseq = 8041
nseq = c(2000, 4000, 6000, 8000, 10000, 12000)
cscaleseq1 = c(1,1.5,2,2.5,3, 3.5, 4)
cscaleseq2 = c(.1, .25, .5, 4.5, 5)
cscaleseq3 = c(.02, .05, .1, .25, .5, .75, 10, 25)/10
cscaleseq = sort(c(cscaleseq1, cscaleseq2, cscaleseq3))
#seedseq = c(1234, 8041, 1546, 1362, 2315, 2141, 5362, 6163, 2326, 6126)
#seedseq = c(seedseq, c(123, 804, 154))

foo <- function(x, y) c(x, y)
setup = outer(dseq, seedseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup, nseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup, cscaleseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)
#setup=setup[,setup[2,]<=10000]
#setup=setup[,setup[2,]>10000]

## Evaluation points
d= 5
	xiseq = seq(1,5,.1*d)
	if(d>3) xiseq = seq(1.2,3.8,.5) #xiseq = seq(1,5,.2*d)
	if(d>1){
		newx = outer(xiseq, xiseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
		if(d>2)
		for(i in 1:(d-2))
			newx = outer(newx, xiseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
	}
	if(d==1){ 
		newx = matrix(xiseq,ncol=1)
	} else newx = t(sapply(newx, cbind))
	
	set.seed(0841)
	newx = newx[seq(1,nrow(newx),4),]
	newx = newx + matrix(rnorm(nrow(newx)*d, 0,.2), ncol=d)

rm_outliers_mean <- function(xtemp){
	uq = quantile(xtemp, .97, na.rm = TRUE)
	lq = quantile(xtemp, .000, na.rm = TRUE)
	xtemp_out1 = which(xtemp<=lq | xtemp >=uq)
	xtemp_out2 = which(xtemp %in% boxplot.stats(xtemp)$out)
	
	xtemp_out = intersect(xtemp_out1, xtemp_out2)
	mean(xtemp[-xtemp_out], na.rm=TRUE)
	
}

thisseed = 8041
thisrun = 0

keepbestbw = list()
res = list()
len = list()
for(n in nseq){
	print(n)
	sigma02tr = var(runif(n*2,0,2))/4+.25*(newx[,1]-3)^2+1

	newx1 = newx[sigma02tr>1.3,]
	
	taufn <- function(x) .5*(x[1]-3)+mean(x[-1])/3
	tau = apply(newx1, 1, taufn)
	
	
	restemp = NULL
	lentemp = NULL
	bestbw = NULL
	
	##############################
	thisseed = 8041	
	##############################
	for(itr in 1:100){
	cat(itr,' '); if(itr%%20==0) cat('\n')
	
	res_acbw = NULL
	len_acbw = NULL
	for(cscale in cscaleseq){		
		# File name to read from
		if(cscale %in% cscaleseq1)
		f = paste0('simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw',cscale*10,'.csv')
		if(cscale %in% cscaleseq2)
		f = paste0('simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw',cscale*100,'.csv')
		if(cscale %in% cscaleseq3)
		f = paste0('simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw',cscale*1000+99,'.csv')
		if(file.exists(f)){
				temp = read.csv(f)
				# Coverage when the bias test is rejected
				cover1 = (temp$tau < temp$Delta2x10 - temp$absDeltax + 2*temp$setauminus) |
						(temp$tau < temp$Delta2x10 - temp$absDeltax + 2*temp$seminus) |
							(temp$tau < temp$betaU - temp$absDeltax + 2*temp$seminus) |
							(temp$tau < temp$betaU - temp$absDeltax + 2*temp$seminus)
				
				cover1 = cover1 & ((temp$tau > temp$Delta2x10 - temp$absDeltax + 2*temp$setauminus) |
									(temp$tau > temp$Delta2x10 - temp$absDeltax + 2*temp$seminus) |
									(temp$tau > temp$betaU - temp$absDeltax + 2*temp$seminus) |
									(temp$tau > temp$betaU - temp$absDeltax + 2*temp$seminus))
				
				# Coverage when the bias test is not rejected
				cover2 = (temp$tau < temp$Delta2x10 + 2*temp$sebetaU) |
						(temp$tau < temp$betaU + 2*temp$sebetaU)
							
				cover2 = cover2 & ((temp$tau > temp$Delta2x10 - 2*temp$sebetaU) |
						(temp$tau > temp$betaU - 2*temp$sebetaU))
						
				# Combine the two based on whether the test is rejected or not
				cover = (temp$absDeltax^2 > pmin(1, temp$cutoffn))*cover1 + (temp$absDeltax^2 < pmin(1, temp$cutoffn))*cover2
				res_acbw = cbind(res_acbw, cover[temp$tau>-.2 & temp$tau<0.1])
				
				len_acbw = cbind(len_acbw, ((temp$absDeltax^2 > pmin(1, temp$cutoffn))*pmin(temp$seminus, temp$setauminus) + 
							(temp$absDeltax^2 < pmin(1, temp$cutoffn))*temp$sebetaU)[temp$tau>-.2 & temp$tau<0.1])
		}
	}
	restemp = cbind(restemp, apply(res_acbw, 1, max))	## iterations in columns
	bestbw = cbind(bestbw, apply(res_acbw, 1, which.max))	## iterations in columns
	lentemp_j = NULL
	for(j in 1:nrow(len_acbw)){
		lentemp_j = c(lentemp_j, min(len_acbw[j,res_acbw[j,]==max(res_acbw[j,])]))
	}
	lentemp = cbind(lentemp, lentemp_j)
	}
	
	
	restemp1 = NULL
	lentemp1 = NULL
	bestbw1 = NULL
	############################## Add more iterations
	for(thisseed in c(8042, 8043, 8044, 8045)){
	##############################
		print(thisseed)
		for(itr in 1:20){
		cat(itr,' '); if(itr%%20==0) cat('\n')
		
		res_acbw = NULL
		len_acbw = NULL
		for(cscale in cscaleseq){		
			# if(cscale %in% cscaleseq1)
			# f = paste0('simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw',cscale*10,'.csv')
			# if(cscale %in% cscaleseq2)
			# f = paste0('simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw',cscale*100,'.csv')
			# if(cscale %in% cscaleseq3)
			f = paste0('simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw',cscale*1000+99,'.csv')
			if(file.exists(f)){
					temp = read.csv(f)
					cover1 = (temp$tau < temp$Delta2x10 - temp$absDeltax + 2*temp$setauminus) |
							(temp$tau < temp$Delta2x10 - temp$absDeltax + 2*temp$seminus) |
								(temp$tau < temp$betaU - temp$absDeltax + 2*temp$seminus) |
								(temp$tau < temp$betaU - temp$absDeltax + 2*temp$seminus)
					
					cover1 = cover1 & ((temp$tau > temp$Delta2x10 - temp$absDeltax + 2*temp$setauminus) |
										(temp$tau > temp$Delta2x10 - temp$absDeltax + 2*temp$seminus) |
										(temp$tau > temp$betaU - temp$absDeltax + 2*temp$seminus) |
										(temp$tau > temp$betaU - temp$absDeltax + 2*temp$seminus))
					
					cover2 = (temp$tau < temp$Delta2x10 + 2*temp$sebetaU) |
							(temp$tau < temp$betaU + 2*temp$sebetaU)
								
					cover2 = cover2 & ((temp$tau > temp$Delta2x10 - 2*temp$sebetaU) |
							(temp$tau > temp$betaU - 2*temp$sebetaU))
							
					cover = (temp$absDeltax^2 > pmin(1, temp$cutoffn))*cover1 + (temp$absDeltax^2 < pmin(1, temp$cutoffn))*cover2
					res_acbw = cbind(res_acbw, cover[temp$tau>-.2 & temp$tau<0.1])
					
					len_acbw = cbind(len_acbw, ((temp$absDeltax^2 > pmin(1, temp$cutoffn))*pmin(temp$seminus, temp$setauminus) + 
								(temp$absDeltax^2 < pmin(1, temp$cutoffn))*temp$sebetaU)[temp$tau>-.2 & temp$tau<0.1])
			}
		}
		restemp1 = cbind(restemp1, apply(res_acbw, 1, max))	## iterations in columns
		bestbw1 = cbind(bestbw1, apply(res_acbw, 1, which.max))	## iterations in columns
		lentemp_j = NULL
		for(j in 1:nrow(len_acbw)){
			lentemp_j = c(lentemp_j, min(len_acbw[j,res_acbw[j,]==max(res_acbw[j,])]))
		}
		lentemp1 = cbind(lentemp1, lentemp_j)
		}
	}
	
	keepbestbw[[which(n==nseq)]] = cbind(bestbw, bestbw1)
	res[[which(n==nseq)]] = cbind(restemp, restemp1)
	len[[which(n==nseq)]] = cbind(lentemp, lentemp1)
	print(summary(apply(res[[which(n==nseq)]],1,mean,na.rm=TRUE)))
	print(summary(apply(len[[which(n==nseq)]],1,rm_outliers_mean)))
	#print(summary(as.numeric(keepbestbw[[which(n==nseq)]])))
	
}

########### Results for the other methods

dseq = 5 #c(5,5,5,5,5,5,5,5,5)
seedseq = as.numeric(sapply((1:5), function(x) x*c(432, 140, 645, 236, 513, 142, 2635, 3616, 6231, 6216, 321, 408, 451)))
nseq = c(2000, 4000, 6000, 8000, 10000, 12000)


foo <- function(x, y) c(x, y)
setup = outer(dseq, seedseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup, nseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
#setup = outer(setup, cscaleseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)


## Evaluation points
d= 5
	xiseq = seq(1,5,.1*d)
	if(d>3) xiseq = seq(1.2,3.8,.5) #xiseq = seq(1,5,.2*d)
	if(d>1){
		newx = outer(xiseq, xiseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
		if(d>2)
		for(i in 1:(d-2))
			newx = outer(newx, xiseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
	}
	if(d==1){ 
		newx = matrix(xiseq,ncol=1)
	} else newx = t(sapply(newx, cbind))
	
	set.seed(0841)
	newx = newx[seq(1,nrow(newx),4),]
	newx = newx + matrix(rnorm(nrow(newx)*d, 0,.2), ncol=d)


f = paste0('simres_se_coverage_n_',4000,'_d_',5,'_seed_', 8041,'_itr_', 1,'_newmethod_v14_r1_',0,'_fixedbw',4*10,'.csv')
temp = read.csv(f)
tauorg = temp$tau

res_cf = res_xl = list()
len_cf = len_xl = list()

for(n in nseq){
	
	print(n)
	sigma02tr = var(runif(n*2,0,2))/4+.25*(newx[,1]-3)^2+1

	newx1 = newx[sigma02tr>1.3,]
	
	taufn <- function(x) .5*(x[1]-3)+mean(x[-1])/3
	tau = apply(newx1, 1, taufn)
	restemp_cf = restemp_xl = NULL
	lentemp_cf = lentemp_xl = NULL
	thisrun = 0
	for(thisseed in seedseq){
		print(thisseed)
		for(itr in 1:5){
			f = paste0('simres_se_coverage_n_',n,'_d_',d,'_seed_', thisseed,'_itr_', itr,'_newmethod_v14_r1_',thisrun,'_fixedbw_othermethods.csv')
			
			if(file.exists(f)){
				temp = read.csv(f)
				# Causal forest method
				cover1 = temp[tauorg>-.2 & tauorg<0.1, 'X.2']
				# X learner method
				cover2 = temp[tauorg>-.2 & tauorg<0.1, 'X.5']
				
				len1 = temp[tauorg>-.2 & tauorg<0.1, 'X.1']
				len2 = (temp[tauorg>-.2 & tauorg<.1, 'X.4'] - temp[tauorg>-.2 & tauorg<0.1, 'X.3'])/4
			
			restemp_cf = cbind(restemp_cf, cover1)
			restemp_xl = cbind(restemp_xl, cover2)
			
			lentemp_cf = cbind(lentemp_cf, len1)
			lentemp_xl = cbind(lentemp_xl, len2)
			}
		}
	}
	
	res_cf[[which(n==nseq)]] = restemp_cf
	len_cf[[which(n==nseq)]] = lentemp_cf
	res_xl[[which(n==nseq)]] = restemp_xl
	len_xl[[which(n==nseq)]] = lentemp_xl
	
	print(summary(apply(res_cf[[which(n==nseq)]],1,mean,na.rm=TRUE)))
	print(summary(apply(len_cf[[which(n==nseq)]],1,mean,na.rm=TRUE)))
	
	print(summary(apply(res_xl[[which(n==nseq)]],1,mean,na.rm=TRUE)))
	print(summary(apply(len_xl[[which(n==nseq)]],1,mean,na.rm=TRUE)))
}





## Plot results


layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2, byrow=TRUE)

layout(mat = layout.matrix,
       heights = c(1), # Heights of the two rows
       widths = c(7, 4)) # Widths of the two columns

resf = NULL
for(n in nseq[1:4]){
	nid = which(n==nseq)
	resf = rbind(resf, 
			cbind(apply(res[[nid]],1,mean,na.rm=TRUE), (nid-1)*4+1))
	resf = rbind(resf, 
			cbind(apply(res_cf[[nid]],1,mean,na.rm=TRUE), (nid-1)*4+2))
	resf = rbind(resf, 
			cbind(apply(res_xl[[nid]],1,mean,na.rm=TRUE), (nid-1)*4+3))	

}

lenf = NULL

for(n in nseq[1:4]){
	nid = which(n==nseq)
	lenf = rbind(lenf, 
			cbind(apply(len[[nid]],1,rm_outliers_mean), (nid-1)*4+1))
	lenf = rbind(lenf, 
			cbind(apply(len_cf[[nid]],1,rm_outliers_mean), (nid-1)*4+2))
	lenf = rbind(lenf, 
			cbind(apply(len_xl[[nid]],1,rm_outliers_mean), (nid-1)*4+3))	

}

#resf = resf[!is.na(lenf[,1]),]
lenf = lenf[!is.na(lenf[,1]),]

par(las=1, mar=c(5.1,4.1,1.1,3.5), oma=c(.2,0,0,1))

boxplot(resf[,1] ~ factor(resf[,2], levels=c(1:15)), xlab='', ylab='Coverage', 
			xaxt='n', yaxt='n', cex.lab=1.35, cex.axis=.95)
#axis(rep(c('Prop.', 'Causal F.', 'X-learn.', ''), 4)[c(1:3,5:7,9:11,13:15)], 
#		at = c(1:3,5:7,9:11,13:15), side=1) 

axis(c('0.0',0.2,0.4,0.6,0.8,0.95,1.0),
		at = c(0.0,0.2,0.4,0.6,0.8,0.95,1.0), side=2, cex.axis=.9) 

text(x = c(1:3,5:7,9:11,13:15)+.2,
     y = par("usr")[3],
     labels = rep(c(expression(Propos.(tau['-'])), 'Causal F.', 'X-learn.', ''), 4)[c(1:3,5:7,9:11,13:15)],
     ## Change the clipping region to fix label drawing.
     xpd = NA,
	srt=35, adj=1.2,
     cex = .9)


mtext('Sample size', side=1,at=8, line=4.2, cex=1.1)
nseqs = sort(nseq)
mtext(nseqs[1], side=1,at=2, line=2.8, cex=1.05)
mtext(nseqs[2], side=1,at=6, line=2.8, cex=1.05)
mtext(nseqs[3], side=1,at=10, line=2.8, cex=1.05)
mtext(nseqs[4], side=1,at=14, line=2.8, cex=1.05)

abline(h=.95)


boxplot(I(4*lenf[,1]) ~ factor(lenf[,2], levels=rev(c(1:15))), horizontal=TRUE,
		yaxt='n', xlab='CI length', ylab='', cex.lab=1.35)
par(las=0)
mtext(nseq[4], side=2,at=2, line=4.8, cex=.95)
mtext(nseq[3], side=2,at=6, line=4.8, cex=.95)
mtext(nseq[2], side=2,at=10, line=4.8, cex=.95)
mtext(nseq[1], side=2,at=14, line=4.8, cex=.95)
mtext("Sample size", side=2, at=8, line=6.1, cex=1.1)


text(x = par("usr")[1],
     y = c(1:3,5:7,9:11,13:15)+.2,
     labels = rev(rep(c(expression(Propos.(tau['-'])), 'Causal F.', 'X-learn.', ''), 4)[c(1:3,5:7,9:11,13:15)]),
     xpd = NA,
	 adj=1.2,
     cex = .9)

