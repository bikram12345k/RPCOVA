setwd('\\\\exasmb.rc.ufl.edu\\blue\\bkarmakar\\bkarmakar\\RPCOVA_ci\\new_ll_fixed')


dseq = c(1,2)
nseq = c(2000, 3000, 4000, 5000)
seedseq = c(1234, 8041, 1546, 1362, 2315, 2141, 5362, 6163, 2326, 6126)
seedseq = c(seedseq, c(123, 804, 154))
seedseq = c(seedseq, c(4321, 1408, 6451, 2362, 5132, 1421, 2635, 3616, 6231, 6216, 321, 408, 451))

seedseq1 = c(1234, 8041, 1546, 1362, 2315, 2141, 5362)
seedseq = c(seedseq1, c(4321, 1408, 6451, 2362, 5132, 1421, 2635, 3616, 6231, 6216, 321, 408, 451))


foo <- function(x, y) c(x, y)
setup = outer(dseq, nseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup, seedseq, FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)
dim(setup)

keepx = list()
cilength = list()
res = list()
for(d in dseq){
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

	keepx[[d]] = apply(newx, 1, function(x) prod(abs(x-3)<1.5))

	for(nid in 1:4){
		res[[d + max(dseq)*(nid-1)]] = matrix(0,nrow=4*sum(keepx[[d]]),ncol=0)	#ncol=23 in the revised code
		cilength[[d + max(dseq)*(nid-1)]] = matrix(0,nrow=3*sum(keepx[[d]]),ncol=0)
	}
	
}


for(id in 1:ncol(setup)){
	n = setup[2,id]; d = setup[1,id]; thisseed = setup[3,id]
	nid = which(nseq == n)
	if(d %in% 1:2 & !(thisseed %in% seedseq1)) next;

	for(itr in 1:ifelse(d %in% 1:2, 140, 40)){
		f = paste0('simres_se_coverage_n_',n,'_d_',d,
			'_seed_', thisseed,'_itr_', itr,'_newmethod_v6.csv')
		
		if(file.exists(f)){
			temp = read.csv(f)

			temp = temp[keepx[[d]]==1,]	# remove xs that are away from the center
			#if( (temp$nhd[1]/n)^(1/d) > .4 ) next;	# remove	cases where bandwidth was incorrect

			if(d<=1)
			cover_maxminus =  (temp$absDeltax^2 > pmin((temp$nhd)^(-.3), 10e-2) )*(( temp$tau > (temp$Delta2x10 - temp$absDeltax) - 2*pmax(temp$seminus, temp$seminus) & 
							temp$tau < (temp$Delta2x10 - temp$absDeltax) + 2*pmax(temp$seminus, temp$seminus) ) |
							( temp$tau > (temp$betaU - temp$absDeltax) - 2*pmax(temp$seminus, temp$seminus) & 
							temp$tau < (temp$betaU - temp$absDeltax) + 2*pmax(temp$seminus, temp$seminus) )
						) + 
				(temp$absDeltax^2 < pmin((temp$nhd)^(-.3), 10e-2))*( temp$tau > (temp$betaU) - 2*temp$sebetaU & 
							temp$tau < (temp$betaU) + 2*temp$sebetaU )
			if(d>=2)
			cover_maxminus =  (temp$absDeltax^2 > pmin((temp$nhd)^(-.3), 10e-2) )*(( temp$tau > (temp$Delta2x10 - temp$absDeltax) - 2*pmax(temp$seminus, temp$seplus) & 
							temp$tau < (temp$Delta2x10 - temp$absDeltax) + 2*pmax(temp$seminus, temp$seplus) ) |
							( temp$tau > (temp$betaU - temp$absDeltax) - 2*pmax(temp$seminus, temp$seplus) & 
							temp$tau < (temp$betaU - temp$absDeltax) + 2*pmax(temp$seminus, temp$seplus) )
						) + 
				(temp$absDeltax^2 < pmin((temp$nhd)^(-.3), 10e-2))*( temp$tau > (temp$betaU) - 2*temp$sebetaU & 
							temp$tau < (temp$betaU) + 2*temp$sebetaU )

		cover_maxplus =  (temp$absDeltax^2 > pmin((temp$nhd)^(-.3), 10e-2) )*(( temp$tau > (temp$Delta2x10 + temp$absDeltax) - 2*pmax(temp$seminus, temp$seplus) & 
							temp$tau < (temp$Delta2x10 + temp$absDeltax) + 2*pmax(temp$seminus, temp$seplus) ) |
							( temp$tau > (temp$betaU + temp$absDeltax) - 2*pmax(temp$seminus, temp$seplus) & 
							temp$tau < (temp$betaU + temp$absDeltax) + 2*pmax(temp$seminus, temp$seplus) )
						) + 
				(temp$absDeltax^2 < pmin((temp$nhd)^(-.3), 10e-2) )*( temp$tau > (temp$betaU) - 2*temp$sebetaU & 
							temp$tau < (temp$betaU) + 2*temp$sebetaU )

			se_maxminus = #pmax(temp$seminus, temp$seplus)
				(temp$absDeltax^2 > pmin((temp$nhd)^(-.3), 10e-2) )*pmax(temp$seminus, temp$seplus) + 
				(temp$absDeltax^2 < pmin((temp$nhd)^(-.3), 10e-2) )*temp$sebetaU

			cilength[[d + max(dseq)*(nid-1)]] = cbind(cilength[[d + max(dseq)*(nid-1)]],
					c(se_maxminus, temp$X.3, (temp$X.6-temp$X.5)/4))

			temp = cbind(cover_maxminus,
							cover_maxplus,
							#cover_maxminusF= temp$tau > (temp$Delta2x10 + temp$absDeltax) - 2*pmax(temp$seminus, temp$seplus) & 
							#temp$tau < (temp$Delta2x10 + temp$absDeltax) + 2*pmax(temp$seminus, temp$seplus),
							temp$X.4, temp$X.7)
			
			res[[d + max(dseq)*(nid-1)]] = cbind(res[[d + max(dseq)*(nid-1)]], as.numeric(temp))
		}
	}

}

#####################################################################
par(las=1, mar=c(5.1,4.1,1.1,4.1), mfrow=c(1,2))

for(d in dseq){

resf = NULL
for(nid in 1:4){
	temp= res[[d + max(dseq)*(nid-1)]]
	resf = rbind(resf, cbind(rowMeans(temp[-((1:(nrow(temp)/4))),]), 
				rep((nid-1)*4+1:3,each=nrow(temp)/4)))
	#temp = resf[,1]
	#temp[temp>.8] = temp[temp>.8]-.4
	#temp[resf[,1]>.4 & resf[,1]<.8] = .4
	#resf[,1] = temp
	print(dim(resf))
}
#print(c(d, nseq[nid]))
#print(aggregate(resf[,1], by =list(resf[,2]), mean))
#print(mean(aggregate(resf[,1], by =list(resf[,2]), mean)[c(1,4,7,10),'x'] ))

#boxplot(resf[,1] ~ factor(resf[,2], levels=c(1:15)), xlab='', ylab='Coverage', 
#			xaxt='n', cex.lab=1.35, cex.axis=.95)

}
#####################################################################
#windows()
layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow=TRUE)

layout(mat = layout.matrix,
       heights = c(1, 1), # Heights of the two rows
       widths = c(7, 4)) # Widths of the two columns


par(las=1, mar=c(5.1,4.1,1.1,3.5), oma=c(.2,0,0,1))

#par(las=1, mar=c(5.2,4.1,1.1,4.1), mfrow=c(2,2))

for(d in dseq){

resf = NULL
for(nid in 1:4){
	temp= res[[d + max(dseq)*(nid-1)]]
	resf = rbind(resf, cbind(rowMeans(temp[-(nrow(temp)/4 + (1:(nrow(temp)/4))),]), 
				rep((nid-1)*4+1:3,each=nrow(temp)/4)))
	#temp = resf[,1]
	#temp[temp>.8] = temp[temp>.8]-.4
	#temp[resf[,1]>.4 & resf[,1]<.8] = .4
	#resf[,1] = temp
	dim(resf)
}
par(las=1)

boxplot(resf[,1] ~ factor(resf[,2], levels=c(1:15)), xlab='', ylab='Coverage', 
			xaxt='n', yaxt='n', cex.lab=1.35, cex.axis=.95)
#axis(rep(c('Prop.', 'Causal F.', 'X-learn.', ''), 4)[c(1:3,5:7,9:11,13:15)], 
#		at = c(1:3,5:7,9:11,13:15), side=1) 

axis(c('0.0',0.2,0.4,0.6,0.8,0.95,1.0),
		at = c(0.0,0.2,0.4,0.6,0.8,0.95,1.0), side=2, cex.axis=.9) 


#axis(c('0.0',0.2,0.3,0.8,0.95,1.0),
#		at = c(0.0,0.2,0.3,0.4,0.55, 0.6), side=2, cex.axis=.9) 


text(x = c(1:3,5:7,9:11,13:15)+.2,
     y = par("usr")[3],
     labels = rep(c(expression(Propos.(tau['-'])), 'Causal F.', 'X-learn.', ''), 4)[c(1:3,5:7,9:11,13:15)],
     ## Change the clipping region to fix label drawing.
     xpd = NA,
	srt=35, adj=1.2,
     cex = .9)

mtext('Sample size', side=1,at=8, line=4.2, cex=1.1)

mtext(nseq[1], side=1,at=2, line=2.8, cex=1.05)
mtext(nseq[2], side=1,at=6, line=2.8, cex=1.05)
mtext(nseq[3], side=1,at=10, line=2.8, cex=1.05)
mtext(nseq[4], side=1,at=14, line=2.8, cex=1.05)

abline(h=.95)
#mtext(0.95, side=2, at=0.95, line = .5, cex=1)


#plot(NULL, xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', ylab="", xlab="", bty='n')

#par(las=1, mar=c(5.1,8.1,1.1,1.1), mfrow=c(2,1))


resf = NULL
for(nid in 1:4){
	temp = cilength[[d + max(dseq)*(nid-1)]]
	resf = rbind(resf, cbind(rowMeans(temp), 
				rep((nid-1)*4+1:3,each=nrow(temp)/3)))
}
boxplot(I(4*resf[,1]) ~ factor(resf[,2], levels=rev(c(1:15))), horizontal=TRUE,
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
     ## Change the clipping region to fix label drawing.
     xpd = NA,
	 adj=1.2,
     cex = .9)

}
par(las=3)
mtext(paste0("d = ",dseq[1]), side=4,at=.8, line=-2, cex=1.25, outer=TRUE)
mtext(paste0("d = ",dseq[2]), side=4,at=.3, line=-2, cex=1.25, outer=TRUE)

