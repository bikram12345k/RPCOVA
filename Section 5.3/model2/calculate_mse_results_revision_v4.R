## CALCULATE MSE results for model 2 in Section 5.3

seedseq <- 2*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721) #, 4241, 1921, 5693, 1218, 8216) 
seedseq <- c(seedseq, 3*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))
seedseq <- c(seedseq, 4*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))
seedseq <- c(seedseq, 5*c(4241, 1921, 5693, 1218, 8216, 1237, 4384, 4191, 2688, 6721))

dseq = 2:5
nseq = c(500, 1000, 1500, 2000)

foo <- function(x, y) c(x, y)
setup = outer(seedseq , dseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = outer(setup , nseq , FUN = Vectorize(foo, SIMPLIFY = FALSE))
setup = sapply(setup, cbind)

# CALCULATE ATT
att <- rep(NA,  8)
for(d in 2:8){
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
		
		att[d] = mean( (y1-y0)[z==1])
}

tauminus_list = list()
for(nid in 1:length(nseq))
	tauminus_list[[nid]] = matrix(0, 0, nseq[nid]+1)
	
mse = NULL

for(id in 1:ncol(setup)){
#for(itr in 1:10){
	print(id)

	thisseed = setup[1,id]
	d = setup[2,id]
	n = setup[3,id]
	nid = which(nseq==n)
	
	## Results file
	f = paste0('linear_stresstest_d_',n,'_',d,'_',thisseed,'_revision_v4_bartps.csv')
	if(file.exists(f)){
		res1 = read.csv(f)

	     res1 = as.matrix(res1)


		tauminus = res1[, 3+(d*n) + 1 + (1:n)]
		tauplus = res1[, 3+(d*n) + 1 + n + 1 + (1:n)]
		
		values = cbind( (att[d] - rowMeans(tauminus))^2,
					(att[d] - rowMeans(tauplus))^2
					)
		

		tauminus_list[[nid]] = rbind(tauminus_list[[nid]], cbind(d, tauminus))

		mse = rbind(mse, cbind(d, n, values)) 
	}
#	}
}


mse_mean = aggregate(mse[,-(1:2)], by = list(mse[,2], mse[,1]), function(x)
			mean(rm_outliers(x), na.rm=TRUE))

mse_mean

##################################################################
##################################################################
## MSE OF THE COMPETING METHODS ##

## These files are created by the slurm call to run 'bashrun_comparison_methods_ate_changesettings_revision_v4.sbatch'
#	The number 9998176 might be different.
fnames = paste0('slurm-9998176_',1:640,'.out')
#fnames = c(fnames, paste0('slurm-55035271_',1:400,'.out'))	# for dimension 8
#fnames = c(fnames, paste0('slurm-55041005_',1:400,'.out'))	# for dimension 8


mseothers = NULL
for(f in fnames){
	# Read from R logs.
	con <- file(f, "r", blocking = FALSE)
	res = readLines(con)

	endf = which(res=="[1] \"final\"")
	endfe = which(sapply(strsplit(res,':'), FUN=function(x) x[1]=='Error'))
	endftm = which(sapply(strsplit(res,' '), FUN=function(x) '***' %in% x))

	if(length(endftm) > 0)
		if(endftm <= 21) next;
	if(length(endfe) > 0)
		if(endfe <= 21) next;

	if(length(endfe) > 0) endf = endfe+1
	if(length(endftm) > 0) endf = endftm+1

	if(length(endf) >0){
		values = res[seq(22,(endf-2),2)]
	} else {
		 values = res[-(1:21)]
		values = values[seq(1, length(values), 2)]	
	}

	values = sapply(strsplit(values,' '), function(x){
				y = as.numeric(x)
				y[(!is.na(y)) | is.nan(y)]
				})
	valuestemp = NULL
	for(i in seq(1, length(values), 2)){
	valuestemp = rbind(valuestemp,
			c(values[[i]], values[[i+1]]))
	}
	mseothers = rbind(mseothers, valuestemp[,-1])
	#mse = rbind(mse, cbind(setup[1],setup[2],setup[3],t(values)))
	close(con)
}


#mse_all=mse
#mse = mse_all[complete.cases(mse_all),]

mse_meanothers = aggregate(mseothers[,-(1:2)], by = list(mseothers[,1], mseothers[,2]), mean, na.rm=TRUE)

mse_meanothers = aggregate(mseothers[,-(1:2)], by = list(mseothers[,1], mseothers[,2]), function(x)
			mean((x), na.rm=TRUE))
			
mse_meanothers


