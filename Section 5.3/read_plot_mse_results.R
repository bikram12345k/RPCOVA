mse = NULL

fnames = paste0('slurm-54246442_',1:240,'.out')

for(f in fnames){
con <- file(f, "r", blocking = FALSE)
res = readLines(con)

endf = which(res=="[1] \"final\"")
endfe = which(sapply(strsplit(res,':'), FUN=function(x) x[1]=='Error'))
endftm = which(sapply(strsplit(res,' '), FUN=function(x) '***' %in% x))

if(length(endftm) > 0)
	if(endftm <= 13) next;
if(length(endfe) > 0)
	if(endfe <= 13) next;

#head(res)

setup = res[7]
setup = as.numeric(strsplit(setup,' ')[[1]][-1])
setup = setup[!is.na(setup)]

if(length(endfe) > 0) endf = endfe+1
if(length(endftm) > 0) endf = endftm+1

if(length(endf) >0){
	values = res[14:(endf-2)]
} else  values = res[-(1:10)]

values = sapply(strsplit(values,' '), function(x){
			y = as.numeric(x[-1])
			y[(!is.na(y)) | is.nan(y)]
			})

mse = rbind(mse, cbind(setup[1],setup[2],setup[3],t(values)))
close(con)
}


rm_outliers <- function(xtemp){
	if(mean(abs(xtemp), na.rm=TRUE) < 100 ) {
	uq = quantile(xtemp, .99, na.rm = TRUE)
	lq = quantile(xtemp, .01, na.rm = TRUE)
	xtemp_out1 = which(xtemp<=lq | xtemp >=uq)
	xtemp_out2 = which(xtemp %in% boxplot.stats(xtemp)$out)
	
	xtemp_out = intersect(xtemp_out1, xtemp_out2)
	xtemp=xtemp[-xtemp_out]
	return(xtemp[!is.nan(xtemp)])

	}
	uq = quantile(xtemp, .95, na.rm = TRUE)
	lq = quantile(xtemp, .15, na.rm = TRUE)
	xtemp_out1 = which(xtemp<=lq*10 | xtemp >=uq/10)
	xtemp_out2 = which(xtemp %in% boxplot.stats(xtemp)$out)
	
	xtemp_out = intersect(xtemp_out1, xtemp_out2)
	xtemp=xtemp[-xtemp_out]
	xtemp[!is.nan(xtemp)]
}


mse_mean = aggregate(mse[,-(1:4)], by = list(mse[,1], mse[,2]), mean, na.rm=TRUE)

mse_mean = aggregate(mse[,-(1:4)], by = list(mse[,1], mse[,2]), function(x)
			mean(rm_outliers(x), na.rm=TRUE))
#round(mse_mean,2)[round(mse_mean,2)[,2] %in% c(2,3,4,5),]
tab = round(mse_mean[mse_mean[,1]%in%c(500,1000,1500,2000) & mse_mean[,2]%in% c(4:7),1:5],2)

tab = round(mse_mean[mse_mean[,1]%in%c(500,1000,1500,2000) & mse_mean[,2]%in% c(5:8),1:5],2)


#mean, na.rm=TRUE)
mse_sd = aggregate(mse[,-(1:4)], by = list(mse[,1], mse[,2]), function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))

mse_median = aggregate(mse[,-(1:4)], by = list(mse[,1], mse[,2]), median, na.rm=TRUE)


table(mse[,1], mse[,2])


##############

#par()$mar
#windows()
par(mfrow=c(1,3), mar=c(5.1,4.5,1.5,1.1), oma=c(0,1,0,0), las=1)
for(d in c(1,5,6)){
	par(las=1)
	mse_mean_sub = mse_mean[mse_mean[,2]==d & mse_mean[,1]>0,]


	plot(1000/mse_mean_sub[,1], mse_mean_sub[,3], type='l', lwd=2, 
			xlab=substitute(a%*%n^-1,list(a=1000)), ylab='', 
			lty=1, cex.lab=1.5, 
			xlim=c(0,max(1000/mse_mean_sub[,1])), 
			ylim = c(0, min(25,max(mse_mean_sub[,3])*1.05)),
			cex.axis=1.25, cex.lab=1.65)
	
	#slope = lm(mse_mean_sub[-(1:2),3]~ I(1/mse_mean_sub[-(1:2),1])-1)$coef

	segments(x0=0, y0=0, x1=min(1000/mse_mean_sub[,1]), y1 = min(mse_mean_sub[,3]), 
			col='grey90', lwd=2)


	mtext(bquote(d == .(d)), side=3, line=-3)
	par(las=0)
	#if(d==1) 
	mtext( substitute(paste("MSE(", widehat(ATE)[a],")"), list(a='–')), side=2, line=3.2)

}

###############################
## comparison methods

mse = NULL

fnames = paste0('slurm-54969968_',1:400,'.out')
fnames = c(fnames, paste0('slurm-55035271_',1:400,'.out'))	# for dimension 8
fnames = c(fnames, paste0('slurm-55041005_',1:400,'.out'))	# for dimension 8



for(f in fnames){
con <- file(f, "r", blocking = FALSE)
res = readLines(con)

endf = which(res=="[1] \"final\"")
endfe = which(sapply(strsplit(res,':'), FUN=function(x) x[1]=='Error'))
endftm = which(sapply(strsplit(res,' '), FUN=function(x) '***' %in% x))

if(length(endftm) > 0)
	if(endftm <= 18) next;
if(length(endfe) > 0)
	if(endfe <= 18) next;

#head(res)

#setup = res[7]
#setup = as.numeric(strsplit(setup,' ')[[1]][-1])
#setup = setup[!is.na(setup)]

if(length(endfe) > 0) endf = endfe+1
if(length(endftm) > 0) endf = endftm+1

if(length(endf) >0){
	values = res[seq(19,(endf-2),2)]
} else {
	 values = res[-(1:18)]
	values = values[seq(1, length(values), 2)]	
}

values = sapply(strsplit(values,' '), function(x){
			y = as.numeric(x)
			y[(!is.na(y)) | is.nan(y)]
			})

mse = rbind(mse, t(values[-1,]))
#mse = rbind(mse, cbind(setup[1],setup[2],setup[3],t(values)))
close(con)
}


#mse_all=mse
#mse = mse_all[complete.cases(mse_all),]

mse_mean = aggregate(mse[,-(1:2)], by = list(mse[,1], mse[,2]), mean, na.rm=TRUE)

mse_mean = aggregate(mse[,-(1:2)], by = list(mse[,1], mse[,2]), function(x)
			mean(rm_outliers(x), na.rm=TRUE))
#round(mse_mean,2)[round(mse_mean,2)[,2] %in% c(2,3,4,5),]
tab2 = round(mse_mean[mse_mean[,1]%in%c(500,1000,1500,2000) & mse_mean[,2]%in% c(5:8),1:4],2)

cbind(tab, tab2)[,-c(5:7)]



