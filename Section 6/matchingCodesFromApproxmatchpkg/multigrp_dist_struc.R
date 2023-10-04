multigrp_dist_struc <-
function(.data, grouplabel, components, wgts){
stopifnot(length(wgts) == length(components))
stopifnot(is.numeric(wgts))
stopifnot(all(wgts>=0))

smahal <-  function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,,drop=FALSE]
    Xt<-X[z==1,,drop=FALSE]
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
	if (!requireNamespace("MASS", quietly = TRUE)) {
		stop("Please load package 'MASS'")
	}
    icov<-MASS::ginv(cv)
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
}

stopifnot(!is.null(names(components)))
compt <- names(components)
compt[tolower(compt) %in% c("prop", "propensity", "mahal", "mahalanobis")] <- ''
compt <- setdiff(unique(compt), '')

for(i in seq(from = 1, length.out = length(compt))){
if(!exists(compt[i], mode='function')) stop('Function ',compt[i], ' not found.\n')
}

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

.data <- cbind(.data, .grouplabel = grouplabel)
grpnames <- unique(grouplabel)
k <- dim(table(grouplabel))
stopifnot(k>1)

calculatedist <- function(id, components, .data, i1, j1, wgts){

.data <- .data[.data[,'.grouplabel'] %in% c(i1, j1),]
grpi1 <- .data[,'.grouplabel'] == i1
names(grpi1) <- rownames(.data)
desc <- names(components)[id]
vars <- components[[id]]

if((!is.null(names(wgts))) & all(names(wgts) %in% names(components))){
wgt <- wgts[desc]
} else wgt <- wgts[id]
# print(c(wgt, dim(.data), sum(grpi1)))

if(tolower(desc) %in% c("prop", "propensity", "logprop", "logpropensity")){
.data <- cbind(.data, grpi1 = grpi1)
proepnsity.model.formula <- as.formula(paste('grpi1 ','~', paste(vars, collapse = " + "), collapse= ' '))

propscore.model = glm(proepnsity.model.formula, data = .data, family = 'binomial')

d <- t(sapply((propscore.model$fitted.values[names(which(grpi1))]),
function(x) abs(x - (propscore.model$fitted.values[names(which(!grpi1))]))))
} else {
if(tolower(desc) %in% c("mahalanobis", "mahal")){
d <- sqrt(smahal(grpi1, .data[, vars]))
  } else {
desc <- eval(parse(text=desc))
#if(exists(desc, mode = "function")){
d <- desc(grpi1, .data[,vars])
#d <- eval(parse(text = paste0(desc, "(grpi1, .data[,vars])" )))
#}
}
}

d <- wgt*d
rownames(d) <- rownames(.data)[grpi1]
colnames(d) <- rownames(.data)[!grpi1]
return(d)
}

distmat <- list()
nparts <- length(components)

for(i in 1:(k-1)){
i1 <- grpnames[i]
for(j in (i+1):k) {
j1 <- grpnames[j]

dlist <- lapply(1:nparts, calculatedist, .data = .data, 
components = components, i1 = i1, j1 = j1, wgts = wgts)
distmat[[paste0(i1, '-', j1)]] <- Reduce("+", dlist)

}
}

return(distmat)
}
