covbalance <-
function(.data, grouplabel, matches, vars, details){

stopifnot(all(vars %in% colnames(.data)))
stopifnot(length(vars)>0)
if(!missing(details)){
stopifnot(is.character(details))
if(is.null(names(details))){
count = 1; noname = TRUE
names(details) <- 1:length(details)
} else noname = FALSE
for(i in seq(from = 1, length.out = length(details))){
if(details[i] %in% c('std.diff', 'std_diff', 'standarized difference')){
 next; } else  if(agrepl('function', details[i])){
if(noname){names(details)[i] <- paste0('genFUN.',count)
count = count + 1 }
 next;
} else if(!exists(details[i], mode='function')){
 stop('Function "',details[i], '" not found.\n')
} else names(details)[i] = details[i]
}
}
if(missing(details))
details <- 'std_diff'

if(is.null(rownames(.data))){
warning('Missing rownames of ".data", assuming 1:nrow(.data)')
rownames(.data) <- 1:nrow(.data)
}

# This part is common to all the functions ##
# parsing 'grouplabel' for convenience  ##

stopifnot(!missing(grouplabel))
stopifnot(is.vector(grouplabel) | is.matrix(grouplabel) | is.data.frame(grouplabel))

gpname <- ''
if(!is.vector(grouplabel)){
if(is.null(colnames(grouplabel)))
grouplabel = apply(grouplabel, 2, function(i) i*grouplabel[,i])
if(!is.null(colnames(grouplabel))){
grouplabeltemp = rep('dummy', nrow(grouplabel))
colLevels <- colnames(grouplabel)
for(i in 1:length(colLevels))
grouplabeltemp[grouplabel[,i]==1] = colLevels[i]
grouplabel = grouplabeltemp
rm(grouplabeltemp)
}
} else {
if(!missing(.data) & all(grouplabel %in% colnames(.data))){
if(length(grouplabel)==1){
grouplabeltemp = .data[,grouplabel]
if(is.numeric(grouplabeltemp))
gpname <- grouplabel
grouplabel = grouplabeltemp
rm(grouplabeltemp)
} else {
grouplabeltemp = rep('dummy', nrow(.data))
colLevels <- grouplabel
for(i in 1:length(colLevels))
grouplabeltemp[.data[,colLevels[i]]==1] = colLevels[i]
grouplabel = grouplabeltemp
rm(grouplabeltemp)
}
}
}

if(is.null(names(grouplabel))){
names(grouplabel) = rownames(.data)
} else grouplabel = grouplabel[rownames(.data)]
stopifnot(all.equal(names(grouplabel), rownames(.data)))

## end of parsing grouplabel ##

.data <- cbind(.data, .grouplabel = grouplabel)
grpnames <- unique(grouplabel)
k <- dim(table(grouplabel))

## create the after match data structure

matches <- as.matrix(matches)
.data.after <- .data[rownames(.data)%in%as.character(matches),]
stopifnot(dim(.data.after)[1]>0)
grouplabel.after <- grouplabel[rownames(.data)%in%as.character(matches)]
grouplabel.after <- grouplabel.after[rownames(.data.after)]

if(is.null(rownames(matches)))
rownames(matches) <- paste('strata', 1:nrow(matches), sep='_')

.strata <- rep(NA, nrow(.data.after))
for(i in 1:nrow(matches))
.strata[match(as.character(matches[i,]), rownames(.data.after), nomatch=0)] <- rownames(matches)[i]

## take within strata averages
.data.after <- aggregate(.data.after[, vars], list(strata=.strata, grp = grouplabel.after), mean, na.rm=TRUE)


## Create the structure of the result
details <- details[which(!(details %in% c('std.diff', 'std_diff', 'standarized difference')))]

if(length(details)>0){
result <- list()

for(i in seq(from = 1, length = length(details))){

variable <- paste0(gsub("[[:space:]]", ".", names(details)[i]),'_before')

evalfunc = paste0('aggregate(.data[,vars], by = list(grp=grouplabel), FUN = ',details[i],')')
temp <- eval(parse(text=evalfunc))
result[[variable ]] <- t(temp[,vars])
colnames(result[[variable ]]) <- temp[,'grp']
result[[variable]] <- result[[variable]][vars, grpnames]
rm(temp)

variable <- paste0(gsub("[[:space:]]", ".", names(details)[i]),'_after')

evalfunc = paste0('aggregate(.data.after[,vars], by = list(grp=.data.after[,"grp"]), FUN = ',details[i],')')
temp <- eval(parse(text=evalfunc))
result[[variable]] <- t(temp[,vars])
colnames(result[[variable]]) <- temp[,'grp']
result[[variable]] <- result[[variable]][vars, grpnames]
rm(temp)
}
}

if('mean' %in% details){
meanbefore <- result$mean_before
meanafter <- result$mean_after
} else {
temp <- aggregate(.data[,vars], by = list(grp=grouplabel), FUN = mean, na.rm=TRUE)
meanbefore <- t(temp[,vars])
colnames(meanbefore) <-  temp[,'grp']
meanbefore <- meanbefore[vars, grpnames]
rm(temp)

temp <- aggregate(.data.after[,vars], by = list(grp=.data.after[,"grp"]), FUN = mean, na.rm=TRUE)
meanafter <- t(temp[,vars])
colnames(meanafter) <- temp[,'grp']
meanafter <- meanafter[vars, grpnames]
rm(temp)
}

if('var' %in% details){
varbefore <- result$var_before
} else { 
temp <- aggregate(.data[,vars], by = list(grp=grouplabel), FUN = var, na.rm=TRUE)
varbefore <- t(temp[,vars])
colnames(varbefore) <-  temp[,'grp']
varbefore <- varbefore[vars, grpnames]
rm(temp)
}

## calculate standized differences

std_diff <- list()

for(i in 1:(k-1))
for(j in (i+1):k){
i1 <- grpnames[i]
j1 <- grpnames[j]

std_var <- sqrt((varbefore[,i1]+varbefore[,j1])/2)
std_diff[[paste0(i1,'-',j1)]] <- 
cbind((meanbefore[,i1]-meanbefore[,j1])/std_var,
(meanafter[,i1]-meanafter[,j1])/std_var)
colnames(std_diff[[paste0(i1,'-',j1)]]) <- c('std_diff_before', 'std_diff_after')
}


if(length(details)>0)
return(list(std_diff=std_diff, details=result))
if(length(details)==0)
return(list(std_diff=std_diff))
}
