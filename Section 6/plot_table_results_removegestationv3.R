## Load Est_INFERENCE_032323_fullv3.RData created by est_inference_fullv3.R
newxold = newx
#newx[newx[,3]>3 & newx[,3]<4,3] = 3
#newx[newx[,2]>11 & newx[,2]<12,3] = 11



#Only use number of visits to 30/35
ids = newx[,1] <= 43 & newx[,2] >= 0 & newx[,2] <= 25 & newx[,3] <= 8
#rep(TRUE, nrow(newx)) #(newx[,1] >= 15 & newx[,1] <= 43 & newx[,2] <= 24 & newx[,3] <= 9 & newx[,3] >= 0 & newx[,3] <= 8)

#windows()
dev.set(2)

#plot(Delta2x10.1, mx1$mean - mx0$mean, xlim=c(-2000,2000))

Delta2x10.1 = Delta2x10
Delta2x10.1[Delta2x10==0] = (mx1$mean - mx0$mean)[Delta2x10==0]

tauminus0 = Delta2x10.1 - sqrt(pmax(Delta2x10.1^2 - Delta2x2, 0))
#tauminus0 = pmin(tauminus0, mx1$mean - mx0$mean)

#tauminus0[which(is.na(tauminus0))] = (mx1$mean - mx0$mean)[which(is.na(tauminus0))]
tauplus0 = Delta2x10.1 + sqrt(pmax(Delta2x10.1^2 - Delta2x2, 0))

Delta2x2.sen = (sigma2x$mean - 9*sigma02x$mean)/( pixnn*(1-pixnn) )
Delta2x.sen = pmax(Delta2x10.1^2 - Delta2x2.sen, 0)
tauminus.sen0 =  Delta2x10.1 - sqrt(Delta2x.sen)


Delta2xeps.1 = pmax(Delta2x10.1^2 - Delta2x2, 0) + .00001


term1 = tauminus0^2*(v12$mean/(pixnn*alphad[3]) + v02$mean/((1-pixnn)*alphad[4]))/Delta2xeps.1
term2 = sigma4lambda2/(4*alphad[2]*Delta2xeps.1*pixnn^2*(1-pixnn)^2)
term3 = tauminus0*(sigma2nu0eta0/sqrt(alphad[2]*alphad[4]) - sigma2nu1eta1/sqrt(alphad[2]*alphad[3]))/(Delta2xeps.1*pixnn*(1-pixnn))

negterm3 = which(term1+term2+term3 < 0)
term3[negterm3] = 2*sqrt(term1*term2)[negterm3]

seminus = sqrt( (term1+term2+term3)*thetaK/fx )/sqrt(n*hprime) # ignore nans


## Create confidence intervals
# For plus
#tauplusest = NULL
#tauminusest = NULL
#setauplus = NULL
setauminus = NULL
cutoffn = pmin(seDelta2x/2, 1)*max(3, (sqrt(n*hprime))^(max(1/3,d/6)))
#cutoffn = (2/sqrt(n*hprime))*max(3, (sqrt(n*hprime))^(1/3))
#deltaused = 0
for(i in 1:nrow(newxold)){
	if(is.nan(seDelta2x[i]) | is.na(pmax(Delta2x10.1^2 - Delta2x2, 0)[i])){
		#setauplus = c(setauplus, sebetaU[i])
		#tauplusest = c(tauplusest, mx1$mean[i] - mx0$mean[i])
		setauminus = c(setauminus, sebetaU[i])
		#tauminusest = c(tauminusest, mx1$mean[i] - mx0$mean[i])
	} else if(pmax(Delta2x10.1^2 - Delta2x2, 0)[i] < cutoffn[i]){
		#setauplus = c(setauplus, sebetaU[i])
		#tauplusest = c(tauplusest, mx1$mean[i] - mx0$mean[i])
		setauminus = c(setauminus, sebetaU[i])
		#tauminusest = c(tauminusest, mx1$mean[i] - mx0$mean[i])
	} else {
		#deltaused = deltaused + 1
		#setauplus = c(setauplus, seplus[i])
		#tauplusest = c(tauplusest, tauplus[i])
		setauminus = c(setauminus, seminus[i])
		#tauminusest = c(tauminusest, tauminus[i])			
	}
}




rm_outliers <- function(xtemp){
	uq = quantile(xtemp, .995, na.rm = TRUE)
	lq = quantile(xtemp, .005, na.rm = TRUE)
	xtemp_out1 = which(xtemp<=lq | xtemp >=uq)
	xtemp_out2 = which(xtemp %in% boxplot.stats(xtemp)$out)
	
	xtemp_out = intersect(xtemp_out1, xtemp_out2)
	xtemp[xtemp_out]=NA
	xtemp
}

tauminus0 = rm_outliers(tauminus0)

length(tauminus)
dim(newx)
#windows()
par(mfrow=c(1,3), mar=c(2,2.5,2,2), oma = c(2,2,0,0), mgp=c(1,.75,0))

for(i in 1:3){
temp = aggregate(tauminus0[ids], by=list(newx[ids,i]), mean, na.rm=TRUE)
temp1 = aggregate(tauminus.sen0[ids], by=list(newx[ids,i]), mean, na.rm=TRUE)
temp2 = aggregate(tauplus0[ids], by=list(newx[ids,i]), mean, na.rm=TRUE)


sigfrac.temp = aggregate(tauminus0[ids]+1.645*pmax(setauminus, setauplus)[ids], by=list(newx[ids,i]), FUN=function(x) mean(x<0, na.rm=TRUE))
sigfrac.temppos = aggregate(tauminus0[ids]-1.645*pmax(setauminus, setauplus)[ids], by=list(newx[ids,i]), FUN=function(x) mean(x>0, na.rm=TRUE))

sigfrac.temp = aggregate(tauminus0[ids]/setauminus[ids], by=list(newx[ids,i]), FUN=function(x) mean(x < -1.645, na.rm=TRUE))
sigfrac.tempU = aggregate((mx1$mean[ids] - mx0$mean[ids])/sebetaU[ids], by=list(newx[ids,i]), FUN=function(x) mean(x < -1.645, na.rm=TRUE))





setemp = aggregate(setauminus[ids]^2, by=list(newx[ids,i]), sum, na.rm=TRUE)
setemp[,2] = sqrt(setemp[,2])
setemp1 = aggregate(setauplus[ids]^2, by=list(newx[ids,i]), sum, na.rm=TRUE)
setemp1[,2] = sqrt(setemp1[,2])

setemp[,2] = pmax(setemp[,2],setemp1[,2])

tempbase = aggregate(mx1$mean[ids] - mx0$mean[ids], by=list(newx[ids,i]), mean, na.rm=TRUE)

head(temp)

if(i==1)  xlimits = c(15, 43)
if(i==2)  xlimits = c(1, 24)
if(i==3)  xlimits = c(0, 8)


plot(temp[temp[,1]<=xlimits[2] & temp[,1]>=xlimits[1],], type='l', lwd=1.45, col='black',ylim=c(-1500, 850), yaxt='n', xlab='', ylab='', xlim=xlimits)
points(tempbase[temp[,1]<=xlimits[2] & temp[,1]>=xlimits[1],], type='l', lwd=1.45, lty=2, col='black',ylim=c(-1500, 850), yaxt='n', xlab='', ylab='', xlim=xlimits)

par(las=1)
axis(label=c(0,1,-750, -500,-250,0,250,500,750), at = c(-1500,-1100,-750,-500,-250,0,250,500,750), side=2)
#mtext('Frac.', side = 2, line = 1.2, cex=1, at=-850, adj=1)
par(las=0)
mtext('Average effect on birthweight', side = 2, line = 3, cex=1.15)
mtext('Prop. sig. neg.', side = 2, line = 1.75, cex=1, at=-1300)


if(i==1)
	mtext("Mother's age at birth (yr)", side = 1, line = 2.22, cex=1.15)
if(i==2)
	mtext("Number of prenatal visits", side = 1, line = 2.22, cex=1.15)
if(i==4)
	mtext('Number of weeks gestation', side = 1, line = 2.22, cex=1.15)
if(i==3)
	mtext('Number of months of prenatal care', side = 1, line = 2.22, cex=1.15)

#points(temp1[temp[,1]<99,], type='l', lwd=1.35, col='blue',ylim=c(-2000, 800))

segments(x0=min(temp[temp[,1]<99,1]), x1=max(temp[temp[,1]<99,1]), y0=-1500+400, y1=-1500+400, col='grey30', lwd=1.45)
segments(x0=min(temp[temp[,1]<99,1]), x1=max(temp[temp[,1]<99,1]), y0=-1500+200, y1=-1500+200, lty=3, col='grey65', lwd=1.45)
segments(x0=min(temp[temp[,1]<99,1]), x1=max(temp[temp[,1]<99,1]), y0=-1500, y1=-1500, lwd=1.45)

gap = diff(range(temp[temp[,1]<99,1]))/350
if(i==3)	gap = .05

for(j in 1:nrow(sigfrac.temp[temp[,1]<99,])){
	segments(y0=-1500,y1=-1500+400*sigfrac.temp[temp[,1]<99,][j,2], x0=sigfrac.temp[temp[,1]<99,][j,1]-gap,
			x1 = sigfrac.temp[temp[,1]<99,][j,1]-gap, lwd=1.45)
	segments(y0=-1500,y1=-1500+400*sigfrac.tempU[temp[,1]<99,][j,2], x0=sigfrac.tempU[temp[,1]<99,][j,1]+gap,
			x1 = sigfrac.tempU[temp[,1]<99,][j,1]+gap, lty=2, lwd=1.45)
	#segments(y0=-1500+200*sigfrac.temp[temp[,1]<99,][j,2],y1=-1500+400*sigfrac.temp[temp[,1]<99,][j,2]+sigfrac.temppos[temp[,1]<99,][j,2], x0=sigfrac.temp[temp[,1]<99,][j,1],
	#		x1 = sigfrac.temp[temp[,1]<99,][j,1], col='red')
}

#points(temp2[temp[,1]<99,], type='l', lwd=1.25, col='blue',ylim=c(-2000, 1000))


x.temp = temp[temp[,1]<99,1]
y2.temp = temp[temp[,1]<99,2]
y1.temp = temp1[temp[,1]<99,2]

#polygon(c(x.temp, rev(x.temp)), c(y2.temp, rev(y1.temp)), col='grey', border='NA')

#points(temp[temp[,1]<99,1], temp[temp[,1]<99,2]-1.96*setemp[temp[,1]<99,2], type='l', lwd=1.25, col='blue',ylim=c(-100, 0))
#points(temp[temp[,1]<99,1], temp[temp[,1]<99,2]+1.96*setemp[temp[,1]<99,2], type='l', lwd=1.25, col='blue',ylim=c(-100, 0))

#points(temp1[temp[,1]<99,], type='l', lwd=1.25, ylim=c(-1000, 1400))
#points(tempbase[temp[,1]<99,], type='l', lwd=1.25, ylim=c(-1000, 1400))

abline(h=0)
#mtext(colnames(x)[i], side=3, line=-3, cex=2)


}

#temp = aggregate(sqrt(Delta2x), by=list(newx[,4]), mean, na.rm=TRUE)
#plot(temp[temp[,1]<99,], type='l', lwd=1.25, ylim=c(-1000, 1000))
#abline(h=0)

plot(x=rep(.75, sum(tauminus0[ids]/setauminus[ids] < - 1.645)),y=tauminus0[ids][tauminus0[ids]/setauminus[ids] < - 1.645], col= 1, xlim=c(.5, 1.5),
	ylim=c(-2000,2000))
points(x=rep(1.25, sum(tauminus0[ids]/setauminus[ids] >= - 1.645)),y=tauminus0[ids][tauminus0[ids]/setauminus[ids] >= - 1.645], col= 2, xlim=c(.5, 1.5),
	ylim=c(-2000,2000))

}


###################################
#windows()
dev.set(3)

tauminus0 = Delta2x10 - sqrt(Delta2x)
tauplus0 = Delta2x10 + sqrt(Delta2x)

Delta2x.sen = Delta2x + 8*sigma02x$mean/( pixnn*(1-pixnn) )
 #(sigma2x$mean - 9*sigma02x$mean)/( pixnn*(1-pixnn) )
#Delta2x.sen = pmax(Delta2x1^2 - Delta2x2.sen, 0)
tauminus.sen0 =  Delta2x10 - sqrt(Delta2x.sen)
tauminus.sen0[is.na(tauminus0)] = NA


length(tauminus)
dim(newx)
#windows()
par(mfrow=c(2,2), mar=c(2,2.5,2,2), oma = c(2,2,0,0), mgp=c(1,.75,0))

for(i in 1:4){
temp = aggregate(tauminus0, by=list(newx[,i]), mean, na.rm=TRUE)
temp1 = aggregate(tauminus.sen0, by=list(newx[,i]), mean, na.rm=TRUE)
temp2 = aggregate(tauplus0, by=list(newx[,i]), mean, na.rm=TRUE)

sigfrac.temp = aggregate(tauplus+1.645*pmax(setauminus, setauplus), by=list(newx[,i]), FUN=function(x) mean(x<0, na.rm=TRUE))


setemp = aggregate(setauminus^2, by=list(newx[,i]), sum, na.rm=TRUE)
setemp[,2] = sqrt(setemp[,2])
setemp1 = aggregate(setauplus^2, by=list(newx[,i]), sum, na.rm=TRUE)
setemp1[,2] = sqrt(setemp1[,2])

setemp[,2] = pmax(setemp[,2],setemp1[,2])

tempbase = aggregate(mx1$mean - mx0$mean, by=list(newx[,i]), mean, na.rm=TRUE)

head(temp)


par(las=1)

plot(temp[temp[,1]<99,], type='l', lwd=1.35, col='black',ylim=c(-4500, 1100), xlab='', ylab='')

#axis(label=c(0,100, -500,-250,0,250,500,750), at = c(-1000,-800,-500,-250,0,250,500,750), side=2)
#mtext('Frac.', side = 2, line = 1.2, cex=1, at=-850, adj=1)
par(las=0)
mtext('Average effect on birthweight', side = 2, line = 3.2, cex=1.15)
#mtext('% sig. neg', side = 2, line = 2.3, cex=1, at=-950)


if(i==1)
	mtext("Mother's age at birth (yr)", side = 1, line = 2.2, cex=1.15)
if(i==2)
	mtext("Number of prenatal visits", side = 1, line = 2.2, cex=1.15)
if(i==3)
	mtext('Number of weeks gestation', side = 1, line = 2.2, cex=1.15)
if(i==4)
	mtext('Month prenatal care began', side = 1, line = 2.2, cex=1.15)


#points(temp1[temp[,1]<99,], type='l', lwd=1.25, col='blue',ylim=c(-2000, 1000))

#segments(x0=min(temp[temp[,1]<99,1]), x1=max(temp[temp[,1]<99,1]), y0=-950+200, y1=-950+200, col='grey30')
#segments(x0=min(temp[temp[,1]<99,1]), x1=max(temp[temp[,1]<99,1]), y0=-950+100, y1=-950+100, lty=3, col='grey')
#segments(x0=min(temp[temp[,1]<99,1]), x1=max(temp[temp[,1]<99,1]), y0=-950, y1=-950)

#for(j in 1:nrow(sigfrac.temp[temp[,1]<99,]))
#	segments(y0=-950,y1=-950+200*sigfrac.temp[temp[,1]<99,][j,2], x0=sigfrac.temp[temp[,1]<99,][j,1],
#			x1 = sigfrac.temp[temp[,1]<99,][j,1])

#points(temp2[temp[,1]<99,], type='l', lwd=1.25, col='blue',ylim=c(-2000, 1000))


x.temp = temp[temp[,1]<99,1]
y2.temp = temp[temp[,1]<99,2]
y1.temp = temp1[temp[,1]<99,2]

polygon(c(x.temp, rev(x.temp)), c(y2.temp, rev(y1.temp)), col='grey', border='NA')

#points(temp[temp[,1]<99,1], temp[temp[,1]<99,2]-1.96*setemp[temp[,1]<99,2], type='l', lwd=1.25, col='blue',ylim=c(-100, 0))
#points(temp[temp[,1]<99,1], temp[temp[,1]<99,2]+1.96*setemp[temp[,1]<99,2], type='l', lwd=1.25, col='blue',ylim=c(-100, 0))

#points(temp1[temp[,1]<99,], type='l', lwd=1.25, ylim=c(-1000, 1400))
#points(tempbase[temp[,1]<99,], type='l', lwd=1.25, ylim=c(-1000, 1400))

abline(h=0)
#mtext(colnames(x)[i], side=3, line=-3, cex=2)


}
