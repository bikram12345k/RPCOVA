kwaymatching <- function (distmat, grouplabel, design, indexgroup = 1, .data, 
    finebalanceVars, exactmatchon, ordering, reorder = FALSE, 
    verbose = TRUE) 
{
    stopifnot(!missing(grouplabel))
    stopifnot(is.vector(grouplabel) | is.matrix(grouplabel) | 
        is.data.frame(grouplabel))
    gpname <- ""
    if (!is.vector(grouplabel)) {
        if (is.null(colnames(grouplabel))) 
            grouplabel = apply(grouplabel, 2, function(i) i * 
                grouplabel[, i])
        if (!is.null(colnames(grouplabel))) {
            grouplabeltemp = rep("dummy", nrow(grouplabel))
            colLevels <- colnames(grouplabel)
            for (i in 1:length(colLevels)) grouplabeltemp[grouplabel[, 
                i] == 1] = colLevels[i]
            grouplabel = grouplabeltemp
            rm(grouplabeltemp)
        }
    }
    else {
        if (!missing(.data)){
            if (all(grouplabel %in% colnames(.data))) {
                if (length(grouplabel) == 1) {
                  grouplabeltemp = .data[, grouplabel]
                  if (is.numeric(grouplabeltemp)) 
                    gpname <- grouplabel
                  grouplabel = grouplabeltemp
                  rm(grouplabeltemp)
                }
                else {
                  grouplabeltemp = rep("dummy", nrow(.data))
                  colLevels <- grouplabel
                  for (i in 1:length(colLevels)) grouplabeltemp[.data[, 
                    colLevels[i]] == 1] = colLevels[i]
                  grouplabel = grouplabeltemp
                  rm(grouplabeltemp)
                }
            }
        }
    }
    if (missing(finebalanceVars) & missing(exactmatchon)) {
        if (missing(.data)) {
            warning("Assuming units labels 1:N\nplease make sure that rownames and column names of distance matrices are names so.")
            .data = data.frame(spholder = rep(NA, length(grouplabel)))
            rownames(.data) = 1:length(grouplabel)
        }
    }
    if (!missing(finebalanceVars) | !missing(exactmatchon)) {
        stopifnot(!missing(.data))
        if (!missing(finebalanceVars)) 
            stopifnot(all(finebalanceVars %in% colnames(.data)))
        if (!missing(exactmatchon)) 
            stopifnot(all(exactmatchon %in% colnames(.data)))
    }
    stopifnot(dim(.data)[1] == length(grouplabel))
    if (!missing(ordering)) {
        orgplevels = ordering
        if (reorder) {
            cat("NOTE: Since 'ordering' is given, ignoring 'reorder=TURE'\n")
            reorder = FALSE
        }
    }
    else orgplevels <- levels(factor(grouplabel))
    if (reorder) {
        orgplevels0 = orgplevels
        stopifnot(!missing(indexgroup))
        while (1) {
            orgplevelstemp = sample(orgplevels)
            if (!all((setdiff(orgplevelstemp, indexgroup) == 
                setdiff(orgplevels, indexgroup)))) {
                orgplevels = orgplevelstemp
                break
            }
        }
    }
    grouplabel <- as.numeric(factor(grouplabel, levels = orgplevels, 
        ordered = T))
    if (is.null(names(grouplabel))) {
        names(grouplabel) = rownames(.data)
    }
    else grouplabel = grouplabel[rownames(.data)]
    grpsizes <- table(grouplabel)
    names(grpsizes) <- orgplevels
    k <- dim(grpsizes)
    #stopifnot(k >= 3)
    stopifnot(all.equal(names(grouplabel), rownames(.data)))
    if (missing(design)) 
        design = rep(1, k)
    stopifnot(length(design) == k)
    stopifnot(all(design%%1 == 0))
    if (is.null(names(design))) {
        if (exists("colLevels")) 
            names(design) = colLevels
        if (!exists("colLevels")) {
            names(design) = orgplevels
            if (reorder) 
                names(design) = orgplevels0
        }
        if (verbose) {
            cat("Using the design structure: \n")
            print(design)
            cat("If this is incorrect please provide appropriate named vector for 'design'\n")
        }
    }
    else if (verbose) {
        cat("Using the design structure: \n")
        print(design)
    }
    if (missing(indexgroup)) 
        indexgroup <- orgplevels[which.min(grpsizes[orgplevels])]
    if (!is.character(indexgroup)) 
        indexgroup <- as.character(indexgroup)
    stopifnot(design[indexgroup] == 1)
    notindexset = orgplevels[which(orgplevels != indexgroup)]
    stopifnot(all(grpsizes[notindexset] >= (design[notindexset] * 
        grpsizes[indexgroup])))
    for (i in 1:(k - 1)){
      for (j in (i + 1):k) {
        i1 = orgplevels[i]
        j1 = orgplevels[j]
        stopifnot(exists(paste0(i1, "-", j1), where = distmat) | 
            exists(paste0(j1, "-", i1), where = distmat))
        if (!exists(paste0(i1, "-", j1), where = distmat) & j1 != 
            indexgroup) {
            distmat[[paste0(i1, "-", j1)]] = t(distmat[[paste0(j1, 
                "-", i1)]])
            distmat[[paste0(j1, "-", i1)]] = NULL
        }
        if (!exists(paste0(j1, "-", i1), where = distmat) & j1 == 
            indexgroup) {
            distmat[[paste0(j1, "-", i1)]] = t(distmat[[paste0(i1, 
                "-", j1)]])
            distmat[[paste0(i1, "-", j1)]] = NULL
        }
      }
    }
    if (missing(finebalanceVars)) 
        stratify = rep("-0-", nrow(.data))
    if (!missing(finebalanceVars)) 
        stratify = paste0("-", apply(.data[, finebalanceVars, 
            drop = FALSE], 1, paste, collapse = "-"))
    names(stratify) = rownames(.data)
    balanceValues = floor(apply(table(stratify, grouplabel), 
        1, function(x) min(x/design[orgplevels])))
    if (missing(exactmatchon)) 
        exactmatchon = NA
    disttemp = distmat[[paste0(indexgroup, "-", notindexset[1])]]
    interMatch <- nrbalancematch(.data, grouplabel == which(indexgroup == 
        orgplevels), grouplabel == which(notindexset[1] == orgplevels), 
        stratify, extmatch = exactmatchon, disttemp, balanceValues, 
        nmatch = design[notindexset[1]])
    interMatch <- interMatch[order(interMatch[, 1]), ,drop=FALSE]
    cost <- sum(apply(interMatch, 1, function(x) disttemp[x[1], 
        x[2]]))
    interMatch <- as.matrix(aggregate(interMatch[, 2], by = list(interMatch[, 
        1]), function(x) x))
    colnames(interMatch) <- paste(gpname, c(indexgroup, paste(notindexset[1], 
        1:design[notindexset[1]], sep = ".")), sep = "_")
    matchedsofar = c(indexgroup, notindexset[1])
    matchedsofar = rep(matchedsofar, design[matchedsofar])
    for (j in seq(from = 2, length = (k - 2))) {
        j1 = notindexset[j]
        disttemp = lapply(1:length(matchedsofar), function(i) distmat[[paste0(matchedsofar[i], 
            "-", j1)]][interMatch[, i], ]/design[matchedsofar[i]])
        disttemp = Reduce("+", disttemp)
        rownames(disttemp) = as.character(interMatch[, 1])
        grouplabeltemp = grouplabel
        grouplabeltemp[grouplabeltemp != which(j1 == orgplevels)] <- 0
        grouplabeltemp[as.character(interMatch[, 1])] <- 1
        secondMatch <- nrbalancematch(.data, grouplabeltemp == 
            1, grouplabeltemp == which(j1 == orgplevels), stratify, 
            extmatch = exactmatchon, disttemp, balanceValues, 
            nmatch = design[j1])
        rm(grouplabeltemp)
        cost <- cost + sum(apply(secondMatch, 1, function(x) disttemp[x[1], 
            x[2]]))
        secondMatch <- as.matrix(aggregate(secondMatch[, 2], 
            by = list(secondMatch[, 1]), function(x) x))
        colnames(secondMatch) <- paste(gpname, c(indexgroup, 
            paste(j1, 1:design[j1], sep = ".")), sep = "_")
        interMatch <- merge(data.frame(interMatch), data.frame(secondMatch), 
            by = colnames(data.frame(interMatch))[1])
        matchedsofar = c(matchedsofar, rep(j1, design[j1]))
    }
    return(list(matches = as.matrix(interMatch), cost = cost))
}
