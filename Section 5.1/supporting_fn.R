## Func to calculate sum_i a_i t(b_i)

cust_outer <- function(a, b){
	if(!is.matrix(a)) a = as.matrix(a)
	if(!is.matrix(b)) b = as.matrix(b)
	n = nrow(a)
	stopifnot(n==nrow(b) | n>=1)
	res = 0
	for(i in 1:n)
		res = res + a[i,]%*%t(b[i,])
	res
}

## Func for constrained local linear estimation for step 3
# const  is  (sigma2x - sigma02)/( pixnn*(1-pixnn) ) 	at 	arg

const_loclin <- function(arg, h, x, d, const, kernel = "gauss"){

	if(!is.matrix(x)) x<-matrix(x,ncol=d)

	if (kernel=="bart") 
	ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
	if (kernel=="gauss") 
	ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
	if (kernel=="uniform") 
	ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
					  return( ans ) }

	argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
	w<-ker((x-argu)/h)/h^d  
	
	## Quadratic optimization function for z=1 group,
	#  written as  x'P1x + q1'x, where x = (intercept1, slope1)

	w1 = w*(z==1)

	P1 = matrix(0, d+1, d+1)
	P1[1,1] = sum(w1)
	P1[-1,-1] = cust_outer(w1*(x-argu), x-argu)
	P1[1,-1] = cust_outer((x-argu), w1)
	P1[-1,1] = t(P1[1,-1])

	q1 = c(-sum(w1*y), -cust_outer(w1*(x-argu), y))

	## Quadratic optimization function for z=0 group,
	#  written as  x'P2x + q2'x, where x = (intercept2, slope2)

	w2 = w*(z==0) 

	P2 = matrix(0, d+1, d+1)
	P2[1,1] = sum(w2)
	P2[-1,-1] = cust_outer(w2*(x-argu), x-argu)
	P2[1,-1] = cust_outer((x-argu), w2)
	P2[-1,1] = t(P2[1,-1])
	
	q2 = c(-sum(w2*y), -cust_outer(w2*(x-argu), y))

	## Combine them into one objective function, for (intercept1, slope1, intercept2, slope2),
	# written as x'Px + q'x

	P = matrix(0, nrow(P1)+nrow(P2), ncol(P1)+ncol(P2))
	P[1:nrow(P1),1:ncol(P1)] = P1
	P[-(1:nrow(P1)),-(1:ncol(P1))] = P2
	P = (P + t(P))/2

	q = c(q1, q2)

	### Set up gurobi optimization problem

	opt_prob <- list()
	opt_prob$Q = rbind(cbind(rbind(cbind(P,0),0),0),0)
	opt_prob$obj  = 2*c(q, 0, 0)	# Two pseudo variables are used to specify the constraint, see below
	opt_prob$modelsense = "min"
	
	## Set the first pseudo variable as (intercept1 - intercept2)
	## Put constraint on the second pseudo to be at least sqrt(const)

	opt_prob$A = matrix( c(1,rep(0,nrow(P1)-1),-1,rep(0,nrow(P2)-1), -1, 0), nrow=1)
	opt_prob$A = rbind(opt_prob$A, c(rep(0, nrow(P)), 0, 1))
	opt_prob$rhs = c(0, sqrt(max(const,0)))
	opt_prob$sense = c('=','>')
	
	## Set the second pseudo variable to abs of the first pseudo variable

	opt_prob$genconabs = list()
	opt_prob$genconabs[[1]] = list()
	opt_prob$genconabs[[1]]$resvar = nrow(P)+2
	opt_prob$genconabs[[1]]$argvar = nrow(P)+1
	
	## Solve optimization with gurobi solver

	result <- gurobi::gurobi(opt_prob)#$x
	
	result$x[1] - result$x[3]
}

const_loclin2 <- function(arg, h, hs, x, d, const, kernel = "gauss"){

	if(!is.matrix(x)) x<-matrix(x,ncol=d)
	if(missing(h)) h = .5
	if(missing(hs)) hs = h
	
	if (kernel=="bart") 
	ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
	if (kernel=="gauss") 
	ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
	if (kernel=="uniform") 
	ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
					  return( ans ) }

	argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
	w<-ker((x-argu)/h)/h^d  
	
	## Quadratic optimization function for z=1 group,
	#  written as  x'P1x + q1'x, where x = (intercept1, slope1)

	w1 = w*(z==1)

	P1 = matrix(0, d+1, d+1)
	P1[1,1] = sum(w1)
	P1[-1,-1] = cust_outer(w1*(x-argu), x-argu)
	P1[1,-1] = cust_outer((x-argu), w1)
	P1[-1,1] = t(P1[1,-1])

	q1 = c(-sum(w1*y), -cust_outer(w1*(x-argu), y))

	## Quadratic optimization function for z=0 group,
	#  written as  x'P2x + q2'x, where x = (intercept2, slope2)

	w2 = ker((x-argu)/hs)/hs^d*(z==0) 

	P2 = matrix(0, d+1, d+1)
	P2[1,1] = sum(w2)
	P2[-1,-1] = cust_outer(w2*(x-argu), x-argu)
	P2[1,-1] = cust_outer((x-argu), w2)
	P2[-1,1] = t(P2[1,-1])
	
	q2 = c(-sum(w2*y), -cust_outer(w2*(x-argu), y))

	## Combine them into one objective function, for (intercept1, slope1, intercept2, slope2),
	# written as x'Px + q'x

	P = matrix(0, nrow(P1)+nrow(P2), ncol(P1)+ncol(P2))
	P[1:nrow(P1),1:ncol(P1)] = P1
	P[-(1:nrow(P1)),-(1:ncol(P1))] = P2
	P = (P + t(P))/2

	q = c(q1, q2)

	### Set up gurobi optimization problem

	opt_prob <- list()
	opt_prob$Q = rbind(cbind(rbind(cbind(P,0),0),0),0)
	opt_prob$obj  = 2*c(q, 0, 0)	# Two pseudo variables are used to specify the constraint, see below
	opt_prob$modelsense = "min"
	
	## Set the first pseudo variable as (intercept1 - intercept2)
	## Put constraint on the second pseudo to be at least sqrt(const)

	opt_prob$A = matrix( c(1,rep(0,nrow(P1)-1),-1,rep(0,nrow(P2)-1), -1, 0), nrow=1)
	opt_prob$A = rbind(opt_prob$A, c(rep(0, nrow(P)), 0, 1))
	opt_prob$rhs = c(0, sqrt(max(const,0)))
	opt_prob$sense = c('=','>')
	
	## Set the second pseudo variable to abs of the first pseudo variable

	opt_prob$genconabs = list()
	opt_prob$genconabs[[1]] = list()
	opt_prob$genconabs[[1]]$resvar = nrow(P)+2
	opt_prob$genconabs[[1]]$argvar = nrow(P)+1
	
	## Solve optimization with gurobi solver

	result <- gurobi::gurobi(opt_prob)#$x
	
	result$x[1] - result$x[3]
}

## Constrained problem solution (fixed for d>1)

const_loclin2_d <- function(arg, h, hs, x, d, const, y, kernel = "gauss"){

	if(!is.matrix(x)) x<-matrix(x,ncol=d)
	if(missing(h)) h = .5
	if(missing(hs)) hs = h
	
	if (kernel=="bart") 
	ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
	if (kernel=="gauss") 
	ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
	if (kernel=="uniform") 
	ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
					  return( ans ) }

	argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
	w<-ker(t(t(x-argu)/h))/prod(h)
	
	## Quadratic optimization function for z=1 group,
	#  written as  x'P1x + q1'x, where x = (intercept1, slope1)

	w1 = w*(z==1)

	P1 = matrix(0, d+1, d+1)
	P1[1,1] = sum(w1)
	P1[-1,-1] = cust_outer(w1*(x-argu), x-argu)
	P1[1,-1] = cust_outer((x-argu), w1)
	P1[-1,1] = t(P1[1,-1])

	q1 = c(-sum(w1*y), -cust_outer(w1*(x-argu), y))

	## Quadratic optimization function for z=0 group,
	#  written as  x'P2x + q2'x, where x = (intercept2, slope2)

	w2 = ker(t(t(x-argu)/hs))/prod(hs)*(z==0) 

	P2 = matrix(0, d+1, d+1)
	P2[1,1] = sum(w2)
	P2[-1,-1] = cust_outer(w2*(x-argu), x-argu)
	P2[1,-1] = cust_outer((x-argu), w2)
	P2[-1,1] = t(P2[1,-1])
	
	q2 = c(-sum(w2*y), -cust_outer(w2*(x-argu), y))

	## Combine them into one objective function, for (intercept1, slope1, intercept2, slope2),
	# written as x'Px + q'x

	P = matrix(0, nrow(P1)+nrow(P2), ncol(P1)+ncol(P2))
	P[1:nrow(P1),1:ncol(P1)] = P1
	P[-(1:nrow(P1)),-(1:ncol(P1))] = P2
	P = (P + t(P))/2

	q = c(q1, q2)

	### Set up gurobi optimization problem

	opt_prob <- list()
	opt_prob$Q = rbind(cbind(rbind(cbind(P,0),0),0),0)
	opt_prob$obj  = 2*c(q, 0, 0)	# Two pseudo variables are used to specify the constraint, see below
	opt_prob$modelsense = "min"
	
	## Set the first pseudo variable as (intercept1 - intercept2)
	## Put constraint on the second pseudo to be at least sqrt(const)

	opt_prob$A = matrix( c(1,rep(0,nrow(P1)-1),-1,rep(0,nrow(P2)-1), -1, 0), nrow=1)
	opt_prob$A = rbind(opt_prob$A, c(rep(0, nrow(P)), 0, 1))
	opt_prob$rhs = c(0, sqrt(max(const,0)))
	opt_prob$sense = c('=','>')
	
	## Set the second pseudo variable to abs of the first pseudo variable

	opt_prob$genconabs = list()
	opt_prob$genconabs[[1]] = list()
	opt_prob$genconabs[[1]]$resvar = nrow(P)+2
	opt_prob$genconabs[[1]]$argvar = nrow(P)+1
	
	## Solve optimization with gurobi solver

	result <- gurobi::gurobi(opt_prob)#$x
	
	result$x[1] - result$x[(1+d)+1]

}


## Local linear regression with bandwidth h calculated at 'arg'

loclin<-function(arg,x,y,h=1,kernel="gauss",type=0)
{
	d<-length(arg)
	n<-length(y)

	if (d>1){

		if (kernel=="bart") 
		   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
		if (kernel=="gauss") 
		   ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
		if (kernel=="uniform") 
		   ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
							  return( ans ) }

		argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
		w<-ker((x-argu)/h)/h^d
		weights<-w/sum(w)

		X<-cbind(matrix(1,n,1),x-argu)
		W<-diag(weights)
		A<-t(X)%*%W%*%X     
		invA<-solve(A,diag(rep(1,d+1))) 
		B<-t(X)%*%W%*%y
		esti<-invA%*%B
		est<-esti[type+1]

	}
	else{  # d==1  #########################################

		if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
		if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }

		x<-matrix(x,length(x),1)
		w<-ker((x-arg)/h)/h^d   
		weights<-w/sum(w)

		X<-cbind(matrix(1,n,1),x-arg)
		W<-diag(c(weights))
		A<-t(X)%*%W%*%X     
		invA<-solve(A,diag(rep(1,d+1))) 
		B<-t(X)%*%W%*%y
		esti<-invA%*%B
		est<-esti[type+1] 

		other<-FALSE
		if (other){
			w<-ker((arg-x)/h); p<-w/sum(w)
			barx<-sum(p*x); bary<-sum(p*y)
			q<-p*(1-((x-barx)*(barx-arg))/sum(p*(x-barx)^2))


			s1<-sum(w*(x-arg))
			s2<-sum(w*(x-arg)^2)
			q<-w*(s2-(x-arg)*s1)/sum(w*(s2-(x-arg)*s1))
		}

	}

	return(est)
}
