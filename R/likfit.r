"lik.fit" <- 	function(data,m,intn,init,control,offset){
	  n <- dim(data)[1]
	  varpos <- 4: (3+m+intn)
	  x <- data[,varpos]
	  varnames <- names(data)[varpos]
	  lbs <- names(x)
	  x <- as.matrix(x)
	  p <- length(varpos)	  
	  
	  d <- data$stat
	  
	  ds <- data$ds
	  h <- data$lambda
	  y <- data$Y-data$start
	  
	  
	  maxiter <- control$maxit
	  if (!missing(init) && !is.null(init)) {
	           if (length(init) != p) 
	               stop("Wrong length for inital values")
	  }
	  else init <- rep(0, p)
	  b <- matrix(init,p,1)
	  b0 <- b
	  	    
	  fit <- mlfit(b,p,x,offset,d,h,ds,y,maxiter,control$epsilon)
	 
	   
	  if (maxiter > 1& fit$nit>=maxiter) 
	                  warning("Ran out of iterations and did not converge")
	              

	  b<-as.vector(fit$b)
	  names(b)<-varnames
	  fit <- list(coefficients=b,var=-solve(fit$sd),iter=fit$nit,loglik=fit$loglik)
	  fit
}


	
	