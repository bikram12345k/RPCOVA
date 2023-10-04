##install relax iv package
drat::addRepo("rrelaxiv", "https://errickson.net/rrelaxiv")
install.packages("rrelaxiv")
#install.packages('approxmatch')
install.packages('optmatch')

#library(approxmatch)
library(optmatch)



source('C:/Users/bkarmakar/Dropbox (UFL)/rt_os_bridging/code/Multi Grp matching samle code/approxmatch/R/covbalance.R')
source('C:/Users/bkarmakar/Dropbox (UFL)/rt_os_bridging/code/Multi Grp matching samle code/approxmatch/R/kwaymatching.R')
source('C:/Users/bkarmakar/Dropbox (UFL)/rt_os_bridging/code/Multi Grp matching samle code/approxmatch/R/multigrp_dist_struc.R')
source('C:/Users/bkarmakar/Dropbox (UFL)/rt_os_bridging/code/Multi Grp matching samle code/approxmatch/R/nrbalancematch.R')
source('C:/Users/bkarmakar/Dropbox (UFL)/rt_os_bridging/code/Multi Grp matching samle code/approxmatch/R/tripletmatching.R')


table(countydata$Mining, countydata$alapacian)
countydata1 = countydata[#countydata$Mining!='Adjacent' & 
				!(countydata$ADMIN_FIPS %in% trt_counties$ADMIN_FIPS[trt_counties$InDSpaper==0]),]

countydata1$Mining[countydata1$Mining=='Adjacent'] = 'No Mining'
countydata1 = countydata1[!(countydata1$Mining=='Mining' & !countydata1$alapacian),]

countydata1$poverty.rate.1990.q = datawizard::categorize(countydata1$poverty.rate.1990, split = "quantile", n_groups = 10)
countydata1$poverty.rate.1990.q2 = datawizard::categorize(countydata1$poverty.rate.1990, split = c(8, 14.10, 17.40, 19, 20))
table(countydata1$poverty.rate.1990.q2, countydata1$Mining)
countydata1$pr_white.q = datawizard::categorize(countydata1$pr_white, split = c(0.9052825,0.9252825, 0.9552825, 0.9750257, 0.9865438, 0.9928119))
countydata1$pr_black.q = datawizard::categorize(countydata1$pr_black, split = "quantile", n_groups = 4)
countydata1$income.1990.q = datawizard::categorize(countydata1$income.1990, split = "quantile", n_groups = 10)


## specify the confounders
vars = c('poverty.rate.1990', 'income.1990', 'noths.1990',  
		'pr_white', 'pr_black', 'smoking_rate')#, 
		#'uninsured.adults', 'primary.care.provider.rate',
		#'single.parent.households')

dist_str <- multigrp_dist_struc(countydata1, 'Mining', 
				components=list(prop=vars,smahal=vars[-c(1,2,3,5)],
					smahal=vars[c(1,2,3)], smahal=c('pr_white')), wgts=c(4,1,4,1))

res = kwaymatching(dist_str, 'Mining', .data=countydata1,
				finebalanceVars = c('poverty.rate.1990.q','pr_white.q','pr_black.q'),
				exactmatchon='poverty.rate.1990.q2') 
			# can include e.g., finebalanceVars='race' or 

countydata1$matches = NA
countydata1[res$matches[,1],'matches'] = 1:nrow(res$matches)
countydata1[res$matches[,2],'matches'] = 1:nrow(res$matches)



details = c('std_diff', 'function(x) round(mean(x,na.rm=TRUE),3)', 
				'function(x) round(median(x,na.rm=TRUE),3)')
names(details) <- c('std_diff', 'mean','median')
covbalance(.data=countydata1, grouplabel='Mining', 
	                 matches = res$matches, vars = vars, details)
					 
View(countydata1[order(countydata1$matches), c('matches',vars)])

sapply(vars, function(v) summary(lm(formula(paste0(v,'~Mining+matches')), data=countydata1))$coef[2,4])

