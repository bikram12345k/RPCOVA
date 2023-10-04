tripletmatching <-
function(distmat, grouplabel, design, indexgroup = 1, .data,
finebalanceVars, exactmatchon){

res1 <- kwaymatching(distmat, grouplabel, design, indexgroup, .data,
finebalanceVars, exactmatchon, reorder=TRUE)
res2 <- kwaymatching(distmat, grouplabel, design, indexgroup = indexgroup, .data,
finebalanceVars, exactmatchon, reorder=FALSE, verbose = FALSE)

if(res1$cost < res2$cost){
cat("Cost: ", res1$cost,"\n")
return(matches = res1$matches)
} else {
cat("Cost: ", res2$cost,"\n")
return(matches = res2$matches)
}
}
