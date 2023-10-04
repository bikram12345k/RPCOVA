setwd("C:/Natality data/NATL1995USPS.AllCnty")

con = file("NATL2010us.AllCnty.txt", "r")


res <- c()
while(TRUE) {
  line = readLines(con, 1)
  if(length(line) == 0) break
  else if(substr(line,107,108) %in% c('TN', 'VA', 'KY', 'WV'))
	res = c(res, line)
}

close(con)



d <- data.frame(state = sapply(res, substr, start=107,stop=108),
fipsco = sapply(res, substr, start=114,stop=116),

#bday = sapply(res, substr, start=22,stop=23),
bmonth = sapply(res, substr, start=19,stop=20),
byear = sapply(res, substr, start=15,stop=18),

#momday = sapply(res, substr, start=67,stop=68),
#mommonth = sapply(res, substr, start=65,stop=66),
#momyear = sapply(res, substr, start=70,stop=73),

age = sapply(res, substr, start=89,stop=90),
 
#adequacy = sapply(res, substr, start=275,stop=275),

nprenvisits =  sapply(res, substr, start=270,stop=271),

gestation = sapply(res, substr, start=451,stop=452),
moncare = sapply(res, substr, start=245,stop=246),	## we can get day count

#fpcareday= sapply(res, substr, start=250,stop=251),
#fpcaremonth = sapply(res, substr, start=248,stop=249),
#fpcareyear = sapply(res, substr, start=252,stop=255),


sex = sapply(res, substr, start=436,stop=436),

birthweight = sapply(res, substr, start=463,stop=466),

cigar = sapply(res, substr, start=294,stop=294),	#?

prev_pret_birth = sapply(res, substr, start= 318, stop=318))

rownames(d) = NULL
d$state[d$state=='KY'] = '21'
d$state[d$state=='WV'] = '54'
d$state[d$state=='VA'] = '51'
d$state[d$state=='TN'] = '47'

d$fipsco = paste0(d$state, d$fipsco)

d$sex[d$sex=='M'] = '0'
d$sex[d$sex=='F'] = '1'

d$cigar[d$cigar=='Y'] = 1
d$cigar[d$cigar=='N'] = 0
d$cigar[d$cigar=='U'] = 9

d$prev_pret_birth[d$prev_pret_birth=='Y'] = 1
d$prev_pret_birth[d$prev_pret_birth=='N'] = 0
d$prev_pret_birth[d$prev_pret_birth=='U'] = 9 

	#drink = sapply(res, substr, start=247,stop=248))	# not available

d = data.frame(sapply(d, as.numeric))


setwd('C:/Natality data/NATL1995USPS.AllCnty')
trt_counties = read.csv('mountaintop_mining_counties.csv')
#trt_counties = trt_counties[trt_counties$InDSpaper==1,]
names(trt_counties)
mining_permit = read.csv('mining_permit_data.csv')
names(mining_permit)
mining_permit$Extended = paste0(mining_permit$ADMIN_NAME,', ',mining_permit$STATE)

## Get fips codes of the treated counties from mining_permit data set
county_fips = mining_permit[,c('Extended', 'ADMIN_FIPS')]
trt_counties = plyr::join(trt_counties, county_fips)

## Identify treated counties
mining_permit$treated = (mining_permit$ADMIN_FIPS %in% trt_counties$ADMIN_FIPS)

table(mining_permit$treated, mining_permit$Mining)	# IMPORTANT
#View(mining_permit[order(mining_permit$treated, mining_permit$Mining, decreasing=TRUE),])


## Now which are the alapacian counties?
alapacian_counties = read.csv('alapacian_counties.csv')

mining_permit$alapacian = mining_permit$ADMIN_FIPS %in% alapacian_counties$FIPS


table(mining_permit$treated, mining_permit$Mining, mining_permit$alapacian)	# IMPORTANT





#####################
poverty = read.csv('Poverty_Rates.csv')
poverty = poverty[,c('FIPS', "Poverty.Rate.1990")]
names(poverty) = c('ADMIN_FIPS', 'poverty.rate.1990')
names(poverty)

countydata = mining_permit
countydata = plyr::join(countydata, poverty)

income = read.csv('income1990short.csv')
names(income) = c('Extended', 'income.1990')

countydata = plyr::join(countydata, income)

education = read.csv('Education.csv')
education = education[,c('FIPS.Code',"Percent.of.adults.with.less.than.a.high.school.diploma..1990")]
names(education) = c('ADMIN_FIPS', 'noths.1990')

countydata = plyr::join(countydata, education)

race = read.csv('race_county.csv')
race$total = rowSums(race[,-(1:2)])
race = race[race$Year==1990,c('ADMIN_FIPS', 'White', 'Black', 'total')]
race$pr_white = race$White/race$total
race$pr_black = race$Black/race$total

race = race[,c('ADMIN_FIPS', 'pr_white', 'pr_black')]

countydata = plyr::join(countydata, race)

head(countydata)
head(d1)
table(d1$cigar)
table(d1$prev_pret_birth)

msmoking = read.csv('county_maternal_smoking_1989_2003.csv')
countydata = plyr::join(countydata, msmoking)

health = read.csv('county_health_data_2010.csv')
health = health[,c('ADMIN_FIPS', 'Uninsured.adults.raw.value',
				'Primary.care.provider.rate.per.100000.population',
				'Single.parent.households.raw.value')]
colnames(health) = c('ADMIN_FIPS', 'uninsured.adults',
					'primary.care.provider.rate',
					'single.parent.households')


countydata = plyr::join(countydata, health)

############
## Use matching code to match counties
#created countydata1
########################

## Births in treated counties
d$treated = d$fipsco %in% trt_counties$ADMIN_FIPS[trt_counties$InDSpaper==1]
d$alapacian = d$fipsco %in% alapacian_counties$FIPS

## Determine control and resistant populations. So far, not using information of whether they are alapacian
#d$control = d$fipsco %in% mining_permit$ADMIN_FIPS[mining_permit$Mining == 'Mining' & !mining_permit$treated] #& mining_permit$alapacian]
#d$resistant = d$fipsco %in% countydata1$ADMIN_FIPS[countydata1$Mining == 'No Mining' & 
#							!countydata1$treated & 
#							!is.na(countydata1$matches)] # & mining_permit$alapacian]

d$control = d$fipsco %in% setdiff(mining_permit$ADMIN_FIPS[mining_permit$Mining == 'Mining' & 
				!mining_permit$treated & mining_permit$alapacian],
				trt_counties$ADMIN_FIPS)

#d$resistant = d$fipsco %in% setdiff(countydata1$ADMIN_FIPS[countydata1$Mining == 'No Mining' & 
#							!countydata1$treated & 
#							countydata1$alapacian],
#				trt_counties$ADMIN_FIPS)
					
d$resistant = d$fipsco %in% setdiff(countydata1$ADMIN_FIPS[countydata1$Mining == 'No Mining' & 
							!countydata1$treated & 
							!is.na(countydata1$matches)],
				trt_counties$ADMIN_FIPS)



d$poptype = 1*d$treated + 2*d$control + 3*d$resistant
#boxplot(d$birthweight ~ d$poptype)
aggregate(d$birthweight, by = list(d$poptype), summary)

d1 = d[d$birthweight<9999.000,]

aggregate(d1$birthweight, by = list(d1$poptype), sd)
table(d1$poptype)

aggregate(d1$birthweight, by = list(d1$poptype), summary)

boxplot(d1$birthweight ~ d1$poptype, ylim=c(2500, 3800))

write.csv(d1, file = 'compiled_data_2010.csv')

