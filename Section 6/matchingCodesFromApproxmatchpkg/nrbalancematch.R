nrbalancematch <- function (cardata.fil, trt_labs, ctrl_labs, stratify, extmatch = NA, 
    distmat, balanceValues, nmatch = 1, solver = 'rlemon') 
{
    if (!requireNamespace("optmatch", quietly = TRUE)) {
        stop("Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.")
    }
	#print(balanceValues)
	#print(which(balanceValues>0))
	#print(stratify)
	#balanceValues <- balanceValues[balanceValues>0]
	#stratify <- stratify[stratify %in% names(balanceValues)]
	#print(balanceValues)
	
	
    treated = seq(from = 1, length = sum(trt_labs))
    treated.names = rownames(cardata.fil)[trt_labs]
    control = seq(from = length(treated) + 1, length = sum(ctrl_labs))
    control.names = rownames(cardata.fil)[ctrl_labs]
    controlLabs = seq(from = length(treated) + length(control) + 
        1, length = length(table(stratify)))
    controlLabs.names = names(balanceValues)
    nnodes = c(treated.names, control.names, controlLabs.names)
    controlLab_extended = unlist(sapply(control.names, function(x){
				tag <- which(nnodes == stratify[x])
				ifelse(length(tag)>0, tag, -9999)
			}))
    names(controlLab_extended) = c()
	
	
	startn = c(rep(treated, length(control)), control, controlLabs)
	endn = c(rep(control, rep(length(treated), length(control))), 
			controlLab_extended, rep(length(nnodes) + 1, length(controlLabs)))
		
	startn = startn[endn!=-9999]
	endn = endn[endn!=-9999]
	controlLab_extended1 <- controlLab_extended
	controlLab_extended <- controlLab_extended[controlLab_extended!=-999]
	
	#print(balanceValues)
    ucap = rep(1, length(startn)-length(controlLabs))
	ucap = c(ucap, rep(length(treated), length(controlLabs)))
	
	b = c(rep(nmatch, length(treated)), rep(0, length(control)), 
			-nmatch * balanceValues, -(nmatch * (length(treated) - 
				sum(balanceValues))))
	
	#print(table(b))
	#print(sum(balanceValues))
	
	
	arc.names.start = c(rep(treated.names, length(control)), 
		control.names[controlLab_extended1!=-999], controlLabs.names)
	arc.names.end = c(rep(control.names, rep(length(treated), 
		length(control))), rep("balance", length(controlLab_extended)), 
		rep("outflow", length(controlLabs)))
	cost = c()
	cost = c(rep(0, (length(startn) - length(controlLabs))), rep(10 * 
			max(distmat, na.rm = T), length(controlLabs)))
	
	distmat = distmat[treated.names, control.names]
    cost[1:(length(treated) * length(control))] = as.numeric(distmat)
    if (!is.na(extmatch)) {
        extmatch <- outer(treated.names, control.names, FUN = function(x, 
            y) cardata.fil[x, extmatch] == cardata.fil[y, extmatch])
        ucap[1:(length(treated) * length(control))] = as.numeric(extmatch)
		cost[1:(length(treated) * length(control))] = cost[1:(length(treated) * length(control))] + 1000 * 
        max(distmat, na.rm = T)*as.numeric(extmatch)
    }
	cost <- cost/max(cost)*1e8
	#print(c("Maximum capacity: ", max(ucap)))
	#print(table(ucap))
#    my.expr <- parse(text = ".Fortran(\"relaxalg\", length(b), as.integer(length(startn)), \n    \t    as.integer(startn), as.integer(endn), as.integer(cost), \n    \t    as.integer(ucap), as.integer(b), x1 = integer(length(startn)), \n    \t    crash1 = as.integer(0), large1 = as.integer(.Machine$integer.max/4), \n    \t    feasible1 = integer(1), NAOK = FALSE, DUP = TRUE, PACKAGE = \"optmatch\")")
#    res <- eval(my.expr)
    
    if(solver == 'rrelaxiv'){
      if(requireNamespace('rrelaxiv', quietly = TRUE)) {
        rout <- rrelaxiv::RELAX_IV(startnodes = as.integer(startn),
                                   endnodes = as.integer(endn),
                                   arccosts = as.integer(cost),
                                   arccapacity = as.integer(ucap),
                                   supply = as.integer(b))
        res <- list(crash = 0,
                           feasible = !all(rout == 0),
                           x = rout)
      } else {
        solver = 'rlemon'
        warning('Package rrelaxiv not available, using rlemon instead.')
      }
    }
    if(solver == 'rlemon'){    
      lout <- rlemon::MinCostFlow(arcSources = as.integer(startn),
                                  arcTargets = as.integer(endn),
                                  arcCapacities = as.integer(ucap),
                                  arcCosts = as.integer(cost),
                                  nodeSupplies = as.integer(b),
                                  numNodes = max(c(startn, endn)),
                                  algorithm = 'CycleCancelling')
      res <- list(crash = 0,
                         feasible = !all(lout[[1]] == 0),
                         x = lout[[1]])
    }else{
      stop(
        'Argument to solver not recognized: please use one of rlemon and rrelaxiv')
    }
    
    if (res$feasible == 0) {
        print("NOT FEASIBLE")
    }
    
	#print(names(res))
    res = res$x
    arc.names.start.res = arc.names.start[res == 1]
    arc.names.end.res = arc.names.end[res == 1]
	match.treated = arc.names.start.res[!(arc.names.end.res %in% c("balance", "outflow"))]
	match.control = arc.names.end.res[!(arc.names.end.res %in% c("balance", "outflow"))]
	
	temp <- cbind(arc.names.start, arc.names.end, startn, endn, cost, ucap, res)
	
	# print(temp[order(as.numeric(temp[,1])),])
	
	# print(c("len of start arc", length(arc.names.start.res)))
	# print(c("unique len of start arc", length(unique(arc.names.start.res))))
	
	# print(c("len of end arc", length(arc.names.end.res)))
	# print(c("unique len of end arc", length(unique(arc.names.end.res))))
	
	# print(length(match.control))
	# print(length(unique(match.control)))
    #match.control = arc.names.start.res[arc.names.end.res == 
    #    "balance"]
    #match.treated = sapply(match.control, function(x) {
    #    trts = arc.names.start.res[arc.names.end.res == x]
    #    trts = trts[which.min(distmat[trts, x])]
    #    ifelse(length(trts) > 0, trts, NA)
    #})
    triplets = cbind(match.treated, match.control)
	#print(triplets)
    triplets
}
