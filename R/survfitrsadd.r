survfit.rsadd <- 
function (object, newdata, se.fit = TRUE, conf.int = 0.95, individual = FALSE, 
      conf.type = c("log", "log-log", "plain", "none")) 
{
   
    call <- match.call()
    Terms <- terms(object)  #to rabis, ce je model mal bl smotan - as.factor ali splines ali svasta
    Terms <- delete.response(Terms)
    newdata <- model.frame(Terms,newdata)
    resp <- list(y=object$y,x=newdata)
    n <- object$n
    nvar <- length(object$coef)
    y <- data$y
    ny <- ncol(y)
    nx <- nrow(newdata)
    nt <- length(object$times)
    temp <- list(n=object$n,time=object$times,call=call,type="right")
    Lambda0 <- object$Lambda0
    Lambda0 <- matrix(Lambda0,ncol=nt,nrow=nrow(newdata),byrow=TRUE)
    rate <- attr(Terms, "specials")$ratetable
    R <- as.matrix(newdata[, rate,drop=FALSE])
    rat <- attributes(object$ratetable)$dimid
    mein <- attributes(newdata[,rate])$dimnames[[2]]
    x <- match(rat,mein)
    R <- R[,x,drop=FALSE]
    newdata <- newdata[,1:(rate-1),drop=FALSE]
    if(any(object$mvalue)>0)newdata <- newdata - matrix(object$mvalue,nrow=nrow(newdata),byrow=TRUE)
    nx <- ncol(newdata)
    #getl <- function(times,data=R,ratetable=object$ratetable){
    #	-log(srvxp.fit(data,times,ratetable))
    #}
    #Lambdap <- sapply(object$times, getl)
   # Lambdap <- NULL
   # for(it in 1:nt){
   # 	Lambdap <- cbind(Lambdap,-log(srvxp.fit(R,object$times[it],object$ratetable)))
   # }
   R <- data.frame(R)
   names(R) <- rat
   Lambdap <- NULL
   for(it in 1:nrow(newdata)){
    Lambdap <- rbind(Lambdap,-log(survexp(~1,data=R[it,,drop=FALSE],times=object$times,ratetable=object$ratetable)$surv))
	}
    ebx <- exp(as.matrix(newdata)%*%object$coef)
    ebx <- matrix(ebx,ncol=nt,nrow=length(ebx))
    Lambda <- Lambdap + Lambda0*ebx
    temp$surv <- t(exp(-Lambda))
    temp$n.event <- rep(1,nt)
    temp$n.risk <- n+1 - cumsum(temp$n.event)
    class(temp) <- c("rs.surv.rsadd", "rs.surv","survfit")
    temp
}