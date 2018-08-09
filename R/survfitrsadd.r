survfit.rsadd <- 
function (formula, newdata, se.fit = TRUE, conf.int = 0.95, individual = FALSE, 
      conf.type = c("log", "log-log", "plain", "none"),...) 
{
   
    call <- match.call()
    Terms <- terms(formula)  #to rabis, ce je model mal bl smotan - as.factor ali splines ali svasta
    Terms <- delete.response(Terms)
    popdata <- newdata
    newdata <- model.frame(Terms,newdata)
    resp <- list(y=formula$y,x=newdata)
    n <- formula$n
    nvar <- length(formula$coef)
    nx <- nrow(newdata)
    nt <- length(formula$times)
    temp <- list(n=formula$n,time=formula$times,call=call,type="right")
    Lambda0 <- formula$Lambda0
    Lambda0 <- matrix(Lambda0,ncol=nt,nrow=nrow(newdata),byrow=TRUE)
    
    rate <- attr(Terms, "specials")$ratetable
    
    #rat <- attributes(formula$ratetable)$dimid
    rat <- names(attributes(formula$ratetable)$dimnames)
    #mein <- attributes(newdata[,rate])$dimnames[[2]]
    mein <- names(popdata)
    x <- match(rat,mein)
    #R <- as.matrix(newdata[, rate, drop = FALSE])
    R <- as.matrix(popdata)
    R <- R[,x,drop=FALSE]
    R <- data.frame(R)
    names(R) <- rat
   
    
    #newdata <- newdata[,1:(rate-1),drop=FALSE]
    labeli <- attr(attr(newdata,"terms"),"term.labels")
    colnami <- colnames(newdata)
    if(length(rate>0)){
      labeli <- labeli[-rate]
      colnami <- colnami[-rate] 
    }
    newdata <- newdata[,match(colnami,labeli),drop=F]
    if(any(formula$mvalue)>0)newdata <- newdata - matrix(formula$mvalue,nrow=nrow(newdata),byrow=TRUE)
    nx <- ncol(newdata)
    #getl <- function(times,data=R,ratetable=formula$ratetable){
    #	-log(srvxp.fit(data,times,ratetable))
    #}
    #Lambdap <- sapply(formula$times, getl)
   # Lambdap <- NULL
   # for(it in 1:nt){
   # 	Lambdap <- cbind(Lambdap,-log(srvxp.fit(R,formula$times[it],formula$ratetable)))
   # }
   
   Lambdap <- NULL
   for(it in 1:nrow(newdata)){
    Lambdap <- rbind(Lambdap,-log(survexp(~1,data=R[it,,drop=FALSE],times=formula$times,ratetable=formula$ratetable)$surv))
	}
    ebx <- exp(as.matrix(formula$coef %*%as.numeric(newdata)))
    ebx <- matrix(ebx,ncol=nt,nrow=length(ebx))
    Lambda <- Lambdap + Lambda0*ebx
    temp$surv <- t(exp(-Lambda))
    temp$n.event <- rep(1,nt)
    temp$n.risk <- n+1 - cumsum(temp$n.event)
    class(temp) <- c("rs.surv.rsadd", "rs.surv","survfit")
    temp
}