"rsmul" <-
function (formula = formula(data), data = parent.frame(), ratetable = survexp.us, 
    int,na.action,init,method="mul",control,...) 
{

    require(survival)
    rform <- rformulate(formula,data,ratetable,na.action)
    
    if(rform$m==0)U <- data.frame(start = rform$start, Y = rform$Y, status = rform$status, rform$R)
    else  U <- data.frame(start = rform$start, Y = rform$Y, status = rform$status, rform$R, rform$X)
    if (!missing(int)) {
    	    U <- U[which(U$start < int * 365.24), ]
            U$status[U$Y >= int * 365.24] <- 0
            U$Y <- pmin(U$Y, int * 365.24)
    }
    else  int <- ceiling(max(rform$Y/365.24))
    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    if(method=="mul"){
	    U <- survsplit(U, cut = (1:int) * 365.24, end = "Y", event = "status", 
		start = "start", episode = "epi")
	    fk <- (attributes(rform$ratetable)$factor != 1)
	    nfk <- length(fk)
	    U[, 4:(nfk + 3)] <- U[, 4:(nfk + 3)] + 365.24 * (U$epi) %*% 
		t(fk)
	    nsk <- dim(U)[1]
	    xx <- srvxp.fit(U[, 4:(nfk + 3)],rep(365.24, nsk), rform$ratetable)
	    lambda <- -log(xx)/365.24
    }
    else if(method=="mul1"){
    	U$id <- 1:dim(U)[1]
    	my.fun <- function(x, attcut, nfk, fk) {
	        intr <- NULL
	        for (i in 1:nfk) {
	            if (fk[i]) {
	                n1 <- max(findInterval(as.numeric(x[3 + i]) + 
	                  as.numeric(x[1]), attcut[[i]]) + 1, 2)
	                n2 <- findInterval(as.numeric(x[3 + i]) + as.numeric(x[2]), 
	                  attcut[[i]])
	                if (n2 > n1 & length(attcut[[i]] > 1)) {
	                  if (n2 > length(attcut[[i]])) 
	                    n2 <- length(attcut[[i]])
	                  intr <- c(intr, as.numeric(attcut[[i]][n1:n2]) - 
	                    as.numeric(x[3 + i]))
	                }
	            }
	        }
	        intr <- sort(unique(c(intr, as.numeric(x[2]))))
	        intr
	    }
	    attcut <- attributes(rform$ratetable)$cutpoints
	    intr <- apply(U[, 1:(3 + nfk)], 1, my.fun, attcut, nfk, fk)
	    dolg <- unlist(lapply(intr, length))
	    newdata <- lapply(U, rep, dolg)
	    stoptime <- unlist(intr)
	    starttime <- c(-1, stoptime[-length(stoptime)])
	    first <- newdata$id != c(-1, newdata$id[-length(newdata$id)])
	    starttime[first] <- newdata$start[first]
	    last <- newdata$id != c(newdata$id[-1], -1)
	    event <- rep(0, length(newdata$id))
	    event[last] <- newdata$status[last]
	    U <- do.call("data.frame", newdata)
	    U$start <- starttime
	    U$Y <- stoptime
	    U$status <- event
	    U[, 4:(nfk + 3)] <- U[, 4:(nfk + 3)] + (U$start) %*% t(fk)
	    nsk <- dim(U)[1]
	    xx <- srvxp.fit(U[, 4:(nfk + 3)],rep(1, nsk), rform$ratetable)
	    lambda <- -log(xx)/1
   }
   else stop("'method' must be one of 'mul' or 'mul1'")
    U$lambda <- log(lambda)
    if (rform$m == 0) 
	fit<-coxph(Surv(start, Y, status) ~1+ offset(lambda), data = U,init=init,control=control,...)
    else{
    	xmat <- as.matrix(U[,(3+nfk+1):(ncol(U)-2)])
	fit<-coxph(Surv(start, Y, status) ~ xmat+offset(lambda), data = U,init=init,control=control,...)
	names(fit[[1]])<-names(U)[(3+nfk+1):(ncol(U)-2)]
    }
    fit$call <- match.call()
    if(length(rform$na.action))fit$na.action <- rform$na.action
    fit
}
