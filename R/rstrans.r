"rstrans" <-
function (formula = formula(data), data = parent.frame(), ratetable = survexp.us, 
    int,na.action,init,control,...) 
{
   
   rform <- rformulate(formula,data,ratetable,na.action)
   	   
    if (!missing(int)) {
        rform$status[rform$Y >= int * 365.24] <- 0
        rform$Y <- pmin(rform$Y, int * 365.24)
    }
    if (rform$type == "counting") {
        start <- 1 - srvxp.fit(rform$R, rform$start,rform$ratetable)
    }
    else start <- rep(0, rform$n)
    
    
    stop <- 1 - srvxp.fit(rform$R, rform$Y, rform$ratetable)
        
    if(rform$m==0)data <- data.frame(start = start, stop = stop, status = rform$status)
    else data <- data.frame(start = start, stop = stop, status = rform$status, rform$X)
    
    
    if (!missing(int)) 
        data <- data[start < stop, , drop = FALSE]
    if(rform$m==0) {
    	if(rform$type=="counting")fit<-coxph(Surv(start,stop,status)~1,data=data,init=init,control=control,...)
    	else fit<-coxph(Surv(stop,status)~1,data=data,init=init,control=control,...)
    }
    else{
    	xmat <- as.matrix(data[,4:ncol(data)])
        fit<-coxph(Surv(start,stop,status)~xmat,data=data,init=init,control=control,...)
        names(fit[[1]])<-names(rform$X)
    }
    fit$call <- match.call()
    if(length(rform$na.action))fit$na.action <- rform$na.action
    fit$data <- data
    return(fit)
}
