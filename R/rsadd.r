"rsadd" <-
function (formula = formula(data), data = parent.frame(),ratetable=survexp.us,int, na.action, method="max.lik",init,control,...) 
{	

	
	call <- match.call()
	if (missing(control))
		control <- glm.control(...)
		
	#organising data
	rform <- rformulate(formula,data,ratetable,na.action)
	
	# cutting data at time int
	if (missing(int)) 
		int <- ceiling(max(rform$Y/365.24))
		if(length(int)==1){
		    if(int<=0) stop("The value of 'int' must be positive ")
		    int<-0:int
		}
	else if(int[1]!=0) stop("The first interval in 'int' must start with 0")    
	
	# fit
	if(method=="glm.bin"|method=="glm.poi"){
		fit <- glmxp(rform=rform,interval=int,method=method, control=control)
	}
	else if(method=="max.lik")
		fit <- maxlik(rform=rform,interval=int,method=method,init=init, control=control) 
	else stop("Unknown method")
		
	# output
	data <- data.frame(start = rform$start, Y = rform$Y, stat = rform$status,rform$R)
	if(rform$m!=0) data <- cbind(data,rform$X)
	fit$call <- call
	fit$formula <- formula
	fit$data <- data
	fit$ratetable <- rform$ratetable
	if(length(rform$na.action))fit$na.action <- rform$na.action
	fit$y <- rform$Y.surv 
	fit
}

