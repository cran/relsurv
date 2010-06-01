rsfitterem<-function(data,b,maxiter,ratetable,tol,bwin,p,cause,Nie){
 pr.time<-proc.time()[3]
if (maxiter<1) stop("There must be at least one iteration run")
n<-nrow(data)
m <- p
dtimes <- which(data$stat==1)			#the positions of event times in data$Y
td <- data$Y[dtimes]				#event times
ntd <- length(td)				#number of event times
utimes <- which(c(1,diff(td))!=0)		#the positions of unique event times among td
utd <- td[utimes]				#unique event times
nutd <- length(utd)				#number of unique event times
udtimes <- dtimes[utimes]			#the positions of unique event times among data$Y
razteg <- function(x){
	# x is a 0/1 vector, the output is a vector of length sum(x), with the corresponding rep numbers
	n <- length(x)
	repu <- rep(1,n)
	repu[x==1] <- 0
	repu <- rev(cumsum(rev(repu)))
	repu <- repu[x==1]
	repu <- -diff(c(repu,0))+1
	if(sum(repu)!=n)repu <- c(n-sum(repu),repu)  #ce je prvi cas censoring, bo treba se kej narest??
	repu
}
rutd <- rep(0,ntd)
rutd[utimes] <- 1
rutd <- razteg(rutd)				#from  unique event times to event times
rtd <- razteg(data$stat)			#from event times to data$Y

a <- data$a[data$stat==1]

if(bwin[1]!=0){
	#the vector of change points for the smoothing bandwidth
	nt4 <- c(1,ceiling(c(nutd*.25,nutd/2,nutd*.75,nutd)))
	if(missing(bwin))bwin <- rep(1,4)
	else bwin <- rep(bwin,4)
	for(it in 1:4){
		bwin[it] <- bwin[it]*max(diff(utd[nt4[it]:nt4[it+1]]))
	}
	while(utd[nt4[2]]<bwin[1]){		# ce je bwin velik, skrajsamo nt4
	       nt4 <- nt4[-2]
	       if(length(nt4)==1)break
	}
	#the smoothing matrix
	krn <- kernerleftch(utd,bwin,nt4)
}



#forming the new dataset
if(p>0){
whtemp <- data$stat==1&cause==2
dataded <- data[data$stat==1&cause==2,]				#events with unknown cause
datacens <- data[data$stat==0|cause<2,]				#censorings or known cause

datacens$cause <- cause[data$stat==0|cause<2]*data$stat[data$stat==0|cause<2]	

databig <- lapply(dataded, rep, 2)				
databig <- do.call("data.frame", databig)
databig$cause <- rep(2,nrow(databig))
nded <- nrow(databig)
databig$cens <- c(rep(1,nded/2),rep(0,nded/2))

datacens$cens <- rep(0,nrow(datacens))
datacens$cens[datacens$cause<2] <- datacens$cause[datacens$cause<2] 
names(datacens) <- names(databig)

databig <- rbind(databig,datacens)

cause <- cause[data$stat==1]


#model matrix for relative survival
xmat <- as.matrix(data[,7:(6+m)])

#ebx at initial values of b
ebx <- as.vector(exp(xmat%*%b))

#model matrix for coxph
modmat <- as.matrix(databig[,7:(6+m)])
varnames <- names(data)[7:(6+m)]
}
else{
	cause <- cause[data$stat==1]
	ebx <- rep(1,n)
}


#for time-dependent data:
starter <- sort(data$start)
starter1<-c(starter[1],starter[-length(starter)])
          
#the values of interest in the cumsums of the obsolete values (there is at least one value - the 1st)
index <- c(TRUE,(starter!=starter1)[-1])
          
starter <- starter[index]
#the number of repetitions in each cumsum difference - needed for s0 calculation
val1 <- apply(matrix(starter,ncol=1),1,function(x,Y)sum(x>=Y),data$Y)
val1 <- c(val1[1],diff(val1),length(data$Y)-val1[length(val1)])




eb <- ebx[data$stat==1]
s0 <- cumsum((ebx)[n:1])[n:1]
        
ebx.st <- ebx[order(data$start)]
s0.st <- ((cumsum(ebx.st[n:1]))[n:1])[index]
s0.st <- rep(c(s0.st,0),val1)
s0 <- s0 - s0.st


#s0 only at times utd
s0 <- s0[udtimes]


#find the corresponding value of Y for each start!=0 - needed for likelihood calculation
start <- data$start
if(any(start!=0)){
	wstart <- rep(NA,n)
	ustart <- unique(start[start!=0])
	for(its in ustart){
		wstart[start==its] <- min(which(data$Y==its))
	}
}


#tale del je zelo sumljiv - kako se racuna likelihood za ties???
difft <- c(data$Y[data$stat==1][1],diff(td))
difft <- difftu <- difft[difft!=0]
difft <- rep(difft,rutd)
a0 <- a*difft

      

if(sum(Nie==.5)!=0)maxit0 <- maxiter
else maxit0<- maxiter - 3
for(i in 1:maxit0){
	
	#Nie is of length ntd, should be nutd, with the values at times being the sum
	nietemp <- rep(1:nutd,rutd)
	Nies <- as.vector(by(Nie,nietemp,sum))  #shorter Nie - only at times utd
	
	lam0u <- lam0 <- Nies/s0				
	#the smooting of lam0        
        if(bwin[1]!=0)lam0s <- krn%*%lam0
        else lam0s <- lam0/difftu
        
        #extended to all event times 
        lam0s <- rep(lam0s,rutd)
        
                
        #compute Nie, only for those with unknown hazard
    	Nie[cause==2] <- as.vector(lam0s*eb/(a+lam0s*eb))[cause==2]
    	
}

if(maxit0!=maxiter & i==maxit0) i <- maxiter
#likelihood calculation - manjka ti se likelihood za nicelni model!!!
#the cumulative hazard     
Lam0  <- cumsum(lam0)
#extended to all event times
Lam0 <- rep(Lam0,rutd)
if(data$stat[1]==0) Lam0 <- c(0,Lam0)
#extended to all exit times
Lam0 <- rep(Lam0,rtd)
#for time dependent covariates: replace by the difference
if(any(start!=0))Lam0[start!=0] <- Lam0[start!=0] - Lam0[wstart[start!=0]]
       
lam0 <- rep(lam0,rutd)
     	
likely0 <- sum(log(a0 + lam0*eb)) - sum(data$ds + Lam0*ebx)
likely <- likely0
tempind <- Nie<=0|Nie>=1
if(any(tempind)){
	if(any(Nie<=0))Nie[Nie<=0] <- tol
	if(any(Nie>=1))Nie[Nie>=1] <- 1-tol
}
	
if(p>0)databig$wei <- c(Nie[cause==2],1-Nie[cause==2],rep(1,nrow(datacens)))


if(maxiter>=1&p!=0){
for(i in 1:maxiter){
        
        if(p>0){
        b00<-b
        if(i==1)fit <- coxph(Surv(start,Y,cens)~modmat,data=databig,weights=databig$wei,init=b00,x=TRUE,iter.max=maxiter)
        else    fit <- coxph(Surv(start,Y,cens)~modmat,data=databig,weights=databig$wei,x=TRUE,iter.max=maxiter)
             
        
        if(any(is.na(fit$coeff))) stop("X matrix deemed to be singular, variable ",which(is.na(fit$coeff)))
        
        b <- fit$coeff
        
        ebx <- as.vector(exp(xmat%*%b))
        }
        else ebx <- rep(1,n)
       
        eb <- ebx[data$stat==1]

        s0 <- cumsum((ebx)[n:1])[n:1]
        
        ebx.st <- ebx[order(data$start)]
        s0.st <- ((cumsum(ebx.st[n:1]))[n:1])[index]
        s0.st <- rep(c(s0.st,0),val1)
        s0 <- s0 - s0.st
        
        #Nie is of length ntd, should be nutd, with the values at times being the sum
        nietemp <- rep(1:nutd,rutd)
        Nies <- as.vector(by(Nie,nietemp,sum))  #shorter Nie - only at times utd
        #s0 only at times utd
        s0 <- s0[udtimes]
        
	lam0u <- lam0 <- Nies/s0				
	
	
	#the cumulative hazard     
        Lam0  <- cumsum(lam0)
        #extended to all event times
        Lam0 <- rep(Lam0,rutd)
        if(data$stat[1]==0) Lam0 <- c(0,Lam0)
        #extended to all exit times
        Lam0 <- rep(Lam0,rtd)
        #for time dependent covariates: replace by the difference
        if(any(start!=0))Lam0[start!=0] <- Lam0[start!=0] - Lam0[wstart[start!=0]]
       	
        #the smooting of lam0        
        if(bwin[1]!=0)lam0s <- krn%*%lam0
        else lam0s <- lam0/difft
        
        #extended to all event times 
        lam0s <- rep(lam0s,rutd)
        
                
        #compute Nie, only for those with unknown hazard
    	Nie[cause==2] <- as.vector(lam0s*eb/(a+lam0s*eb))[cause==2]
        
                        
        #likelihood calculation - manjka ti se likelihood za nicelni model!!!
       
       	lam0 <- rep(lam0,rutd)
     	
              
        likely <- sum(log(a0 + lam0*eb)) - sum(data$ds + Lam0*ebx)
	
	if(p>0){
        	tempind <- Nie<=0|Nie>=1
	        if(any(tempind)){
	                if(any(Nie<=0))Nie[Nie<=0] <- tol
	                if(any(Nie>=1))Nie[Nie>=1] <- 1-tol
	                #if(which(tempind)!=nev)warning("Weights smaller than 0")                       
	                #if(any(is.na( match(which(tempind),c(1,nev)) )))browser()                      
        	}
        	if(nded==0) break()
        	databig$wei[1:nded] <- c(Nie[cause==2],1-Nie[cause==2])
        	bd <- abs(b-b00)
		if(max(bd)< tol) break()    
        }
        #early stopping time for no covariates???
}
}
iter <- i
#if (maxiter > 1& iter>=maxiter) 
#        warning("Ran out of iterations and did not converge")
if(p>0){
if(nded!=0){
	resi <- resid(fit,type="schoenfeld")
	if(!is.null(dim(resi)))resi <- resi[1:(nded/2),]
	else resi <- resi[1:(nded/2)]
	swei <- fit$weights[1:(nded/2)]
	
	if(is.null(dim(resi))) fishem <- sum((resi^2*swei*(1-swei)))
	else {
		fishem <- apply(resi,1,function(x)outer(x,x))
		fishem <- t(t(fishem)*swei*(1-swei))
		fishem <- matrix(apply(fishem,1,sum),ncol=m)
	}
}
else fishem <- 0
fishcox <- solve(fit$var)
fisher <- fishcox - fishem
fit$var <- solve(fisher)
names(fit$coefficients)<-varnames
fit$lambda0 <- lam0s
}
else fit <- list(lambda0 = lam0s)
fit$lambda0 <- fit$lambda0[utimes]
fit$Lambda0 <- Lam0[udtimes]
fit$times <- utd
fit$Nie <- Nie
fit$bwin <- bwin
fit$iter <- i
class(fit) <- c("rsadd",class(fit))
fit$loglik <- c(likely0,likely)
fit$lam0.ns <- lam0u					
fit
}



em <- function (rform, init, control, bwin) 
{
    data <- rform$data
    n <- nrow(data)
    p <- rform$m
    id <- order(data$Y)
    rform$cause <- rform$cause[id]
    data <- data[id, ]
    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    nev <- length(data$Y[data$stat == 1])
    data$a <- rep(NA, n)
    xx <- srvxp.fit(data[, 4:(nfk + 3)], data$Y - data$start, rform$ratetable)
    data$ds <- -log(xx)
    data1 <- data
    data1[, 4:(nfk + 3)] <- data[, 4:(nfk + 3)] + data$Y %*% t(fk)
    xx <- srvxp.fit(data1[data1$stat == 1, 4:(nfk + 3)], rep(1, 
        nev), rform$ratetable)
    data$a[data$stat == 1] <- -log(xx)
    
    if (p > 0) {
        if (!missing(init) && !is.null(init)) {
            if (length(init) != p) 
                stop("Wrong length for inital values")
        }
        else init <- rep(0, p)
        beta <- matrix(init, p, 1)
    }
    pr.time<-proc.time()[3]
    
    Nie <- rep(.5,sum(data$stat==1))
    Nie[rform$cause[data$stat==1]<2] <-  rform$cause[data$stat==1][rform$cause[data$stat==1]<2]
    
    if(missing(bwin))bwin <- -1
    if(bwin<0){

	if(p>0)data1 <- data[,-c(7:(6+p))]    
	else  data1 <- data
	nfk <- length(attributes(rform$ratetable)$dimid)
	names(data)[4:(3+nfk)] <- attributes(rform$ratetable)$dimid
	expe <- rs.survo(Surv(Y,stat)~1,data,ratetable=rform$ratetable,method="conditional")
	esurv <- -log(expe$surv[expe$n.event!=0])
	if(esurv[length(esurv)]==Inf)esurv[length(esurv)] <-  esurv[length(esurv)-1]
	x <- seq(.1,3,length=5)
	dif <- rep(NA,5)
	options(warn=-1)
	diter <- max(round(max(data$Y)/356.24),3)
	for(it in 1:5){
		fit <- rsfitterem(data1,NULL,diter,rform$ratetable,control$epsilon,x[it],0,rform$cause,Nie)
		dif[it] <- sum((esurv-fit$Lambda0)^2)
	}
	wh <- which.min(dif)
	if(wh==1)x <- seq(x[wh],x[wh+1]-.1,length=5)
	else if(wh==5)x <- c(x, max(data$Y)/ max(diff(data$Y)))
	if(wh!=1)
	x <- seq(x[wh-1]+.1,x[wh+1]-.1,length=5)
    	dif <- rep(NA,5)
    	
   	for(it in 1:5){
		fit <- rsfitterem(data1,NULL,diter,rform$ratetable,control$epsilon,x[it],0,rform$cause,Nie)
		dif[it] <- sum((esurv-fit$Lambda0)^2)
	}
	options(warn=0)
	Nie <- fit$Nie
	bwin <- x[which.min(dif)]
    }
       
        fit <- rsfitterem(data, beta, control$maxit, rform$ratetable, 
                 control$epsilon, bwin, p, rform$cause,Nie)
        
        Nie <- rep(0,nrow(data))
        Nie[data$stat==1] <- fit$Nie 
        fit$Nie <- Nie[order(id)]
         fit$bwin <- list(bwin=fit$bwin,bwinfac=bwin)
         fit
     }


rsadd <- function (formula = formula(data), data = parent.frame(), ratetable = survexp.us, 
    int, na.action, method = "max.lik", init, bwin, centered = FALSE, 
    cause, control, ...) 
{
    call <- match.call()
    if (missing(control)) 
        control <- glm.control(...)
    #if (missing(cause)) 
    #    cause <- rep(FALSE, nrow(data))
    #if (!is.logical(cause)) 
    #    stop("Variable 'cause' must be logical")
    if (missing(cause)) 
           cause <- rep(2, nrow(data))
    if (length(cause) != nrow(data)) 
            stop("Length of cause does not match data dimensions")
        data$cause <- cause
    rform <- rformulate(formula, data, ratetable, na.action, 
        int, centered, cause)
    if (method == "EM") {
        if (!missing(int)) {
            if (length(int) > 1 | any(int <= 0)) 
                stop("Invalid value of 'int'")
        }
    }
    else {
        if (missing(int)) 
            int <- c(0,ceiling(max(rform$Y/365.24)))
        if (length(int) == 1) {
            if (int <= 0) 
                stop("The value of 'int' must be positive ")
            int <- 0:int
        }
        else if (int[1] != 0) 
            stop("The first interval in 'int' must start with 0")
    }
    method <-  match.arg(method,c("glm.bin","glm.poi","max.lik","EM"))

    if (method == "glm.bin" | method == "glm.poi") 
        fit <- glmxp(rform = rform, interval = int, method = method, 
            control = control)
    else if (method == "max.lik") 
        fit <- maxlik(rform = rform, interval = int, init = init, 
            control = control)
    else if (method == "EM") 
        fit <- em(rform, init, control, bwin)
    fit$call <- call
    fit$formula <- formula
    fit$data <- rform$data
    fit$ratetable <- rform$ratetable
    fit$n <- nrow(rform$data)
    if (length(rform$na.action)) 
        fit$na.action <- rform$na.action
    fit$y <- rform$Y.surv
    fit$method <- method
    if (method == "EM") {
        if (!missing(int)) 
            fit$int <- int
        else fit$int <- ceiling(max(rform$Y[rform$status == 1])/365.24)
        fit$terms <- rform$Terms
        if(centered)fit$mvalue <- rform$mvalue
    }
    if (rform$m > 0) 
        fit$linear.predictors <- as.matrix(rform$X) %*% fit$coef[1:ncol(rform$X)]
    fit
}

rformulate <- function (formula, data = parent.frame(), ratetable, na.action, 
    int, centered, cause) 
{
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$ratetable <- m$int <- m$centered <- NULL
    Terms <- if (missing(data)) 
        terms(formula, "ratetable", "cause")
    else terms(formula, "ratetable", "cause", data = data)
    rate <- attr(Terms, "specials")$ratetable
    if (length(rate) > 1) 
        stop("Can have only 1 ratetable() call in a formula")
    if (length(rate) == 0) {
        xx <- function(x) formula(x)
        if (is.ratetable(ratetable)) 
            varlist <- attr(ratetable, "dimid")
        else stop("Invalid rate table")
        ftemp <- deparse(formula)
         ftemp <- paste(ftemp,collapse="")
        formula <- xx(paste(ftemp, "+ ratetable(", paste(varlist, 
            "=", varlist, collapse = ","), ")"))
        Terms <- if (missing(data)) 
            terms(formula, "ratetable")
        else terms(formula, "ratetable", data = data)
        rate <- attr(Terms, "specials")$ratetable
    }
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    n <- nrow(m)
    Y <- model.extract(m, "response")
    if (!is.Surv(Y)) 
        stop("Response must be a survival object")
    Y.surv <- Y
    if (attr(Y, "type") == "right") {
        type <- attr(Y, "type")
        status <- Y[, 2]
        Y <- Y[, 1]
        start <- rep(0, n)
        ncol0 <- 2
    }
    else if (attr(Y, "type") == "counting") {
        type <- attr(Y, "type")
        status <- Y[, 3]
        start <- Y[, 1]
        Y <- Y[, 2]
        ncol0 <- 3
    }
    else stop("Illegal response value")
    if (any(c(Y, start) < 0)) 
        stop("Negative follow up time")
    if(max(Y)<30)
    	warning("The event times must be expressed in days! (Your max time in the data is less than 30 days) \n")
    if (is.ratetable(ratetable)) {
        israte <- TRUE
        rtemp <- match.ratetable(m[, rate], ratetable)
        if(is.null(attributes(ratetable)$factor))attributes(ratetable)$factor <- attributes(ratetable)$type==1
        rtorig <- attributes(ratetable)
        nrt <- length(rtorig$dimid)
        R <- rtemp$R
        if (!is.null(rtemp$call)) {
            ratetable <- eval(parse(text = rtemp$call))
        }
        #checking if the ratetable variables are given in days
        wh.age <- which(attributes(ratetable)$dimid=="age")
        wh.year <- which(attributes(ratetable)$dimid=="year")
        if(length(wh.age)>0){
        	if(max(rtemp$R[,wh.age])<150& median(diff(attributes(ratetable)$cutpoints[[wh.age]]))>12)
        	warning("Age in the ratetable part of the formula must be expressed in days! \n (Your max age is less than 150 days) \n")
        }
        if(length(wh.year)>0){
        	if(min(rtemp$R[,wh.year])>1850 & max(rtemp$R[,wh.year])<2020&class(attributes(ratetable)$cutpoints[[wh.year]])=="date")
        	warning("The calendar year must be expressed in days since 1.1.1960! \n (Your variable seems to be expressed in years) \n")
        }
        #checking if one of the continuous variables is fixed:
        if(nrt!=ncol(R)){
        	nonex <- which(is.na(match(rtorig$dimid,attributes(ratetable)$dimid)))
        	for(it in nonex){
        		if(rtorig$type[it]!=1)warning(paste("Variable ",rtorig$dimid[it]," is held fixed even though it changes in time in the population tables. \n (You may wish to set a value for each individual and not just one value for all)",sep=""))
        	}
        }
        
    }
    else stop("Invalid ratetable argument")
    if (rate == 2) {
        X <- NULL
        mm <- 0
    }
    else {
        #if (length(rate == 1)) {
        #    formula[[3]] <- formula[[3]][[2]]
        #}
        X <- as.data.frame(model.matrix(formula, data = data))[, 
            -1, drop = FALSE]
        mm <- ncol(X)
        X <- X[,1:(mm-nrt),drop=FALSE]
        mm <- ncol(X)
        #mn <- names(X)
        #X <- data.frame(do.call("cbind",m))
	#X <- X[,(ncol0+1):(ncol0+mm),drop=FALSE]
        #names(X) <- mn
    }
    mvalue <- rep(0,mm)
    if (!missing(centered)) {
        if (mm != 0 & centered == TRUE) {
            mvalue <- apply(as.matrix(X),2,mean)
            X <- apply(as.matrix(X), 2, function(x) x - mean(x))
       }
    }
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0) 
        rep(0, n)
    else if (tt == 1) 
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    keep <- Y > start
    cause <- model.extract(m, "cause")
    #status[cause==0] <- 0
    if (!missing(int)) {
        int <- max(int)
        status[Y > int * 365.24] <- 0
        Y <- pmin(Y, int * 365.24)
        keep <- keep & (start < int * 365.24)
    }
    if (any(start > Y) | any(Y < 0)) 
        stop("Negative follow-up times")
    X <- X[keep, , drop = FALSE]
    Y <- Y[keep]
    start <- start[keep]
    status <- status[keep]
    R <- R[keep, ]
    offset <- offset[keep]
    Y.surv <- Y.surv[keep, , drop = FALSE]
    cause <- cause[keep]
    n <- sum(keep)
    data <- data.frame(start = start, Y = Y, stat = status, R)
    if (mm != 0) 
        data <- cbind(data, X)
    out <- list(data = data, R = R, status = status, start = start, 
        Y = Y, X = as.data.frame(X), m = mm, n = n, type = type, Y.surv = Y.surv, 
        Terms = Terms, ratetable = ratetable, offset = offset, formula=formula,
        cause = cause,mvalue=mvalue)
    na.action <- attr(m, "na.action")
    if (length(na.action)) 
        out$na.action <- na.action
    out
}

maxlik <- function (rform, interval, subset, init, control) 
{
    data <- rform$data
    max.time <- max(data$Y)/365.24
    if (max.time < max(interval)) 
        interval <- interval[1:(sum(max.time > interval) + 1)]
    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    data <- cbind(data, offset = rform$offset)
    data <- survsplit(data, cut = interval[-1] * 365.24, end = "Y", 
        event = "stat", start = "start", episode = "epi", interval = interval)
    del <- which(data$start==data$Y)   
    if(length(del))    data <- data[-del,]
    offset <- data$offset
    data$offset <- NULL
    d.int <- diff(interval)
    data[, 4:(nfk + 3)] <- data[, 4:(nfk + 3)] + data$start %*% 
        t(fk)
    data$lambda <- rep(0, nrow(data))
    nsk <- nrow(data[data$stat == 1, ])
    xx <- srvxp.fit(data[data$stat == 1, 4:(nfk + 3)] + (data[data$stat == 
        1, ]$Y - data[data$stat == 1, ]$start) %*% t(fk), rep(1, 
        nsk), rform$ratetable)
    data$lambda[data$stat == 1] <- -log(xx) * 365.24
    xx <- srvxp.fit(data[, 4:(nfk + 3)], data$Y - data$start, 
        rform$ratetable)
    data$epi <- NULL
    data$ds <- -log(xx)
    data$Y <- data$Y/365.24
    data$start <- data$start/365.24
    data <- data[, -(4:(3 + nfk))]
    intn <- length(interval[-1])
    m <- rform$m
    p <-  m + intn
    if (!missing(init) && !is.null(init)) {
        if (length(init) != p) 
            stop("Wrong length for inital values")
    }
    else init <- rep(0, p)
    if(m>0){
    	init0 <- init[-(1:m)]
    	data1 <- data[,-(4:(3+m))]
    }
    else{
    	init0 <- init
    	data1 <- data
    }
    fit0 <- lik.fit(data1, 0, intn, init0,  control, offset)
    if(m>0){
       	init[-(1:m)] <- fit0$coef
    	fit <- lik.fit(data, m, intn, init,  control, offset)
    }
    else fit <- fit0
    fit$int <- interval
    class(fit) <- "rsadd"
    fit
}

lik.fit <- function (data, m, intn, init, control, offset) 
{
    n <- dim(data)[1]
    varpos <- 4:(3 + m + intn)
    x <- data[, varpos]
    varnames <- names(data)[varpos]
    lbs <- names(x)
    x <- as.matrix(x)
    p <- length(varpos)
    d <- data$stat
    ds <- data$ds
    h <- data$lambda
    y <- data$Y - data$start
    maxiter <- control$maxit
    if (!missing(init) && !is.null(init)) {
        if (length(init) != p) 
            stop("Wrong length for inital values")
    }
    else init <- rep(0, p)
    b <- matrix(init, p, 1)
    b0 <- b
    fit <- mlfit(b, p, x, offset, d, h, ds, y, maxiter, control$epsilon)
    if (maxiter > 1 & fit$nit >= maxiter) {
	values <- apply(data[data$stat==1,varpos],2,sum)
	problem <- which.min(values)
	outmes <- "Ran out of iterations and did not converge" 
	if(values[problem]==0)tzero <- ""
	else tzero <- "only "
	if(values[problem]<5)outmes <- paste(outmes, "\n This may be due to the fact that there are ",tzero, values[problem], " events on interval",strsplit(names(values)[problem],"fu")[[1]][2],"\n You can use the 'int' argument to change the follow-up intervals in which the baseline excess hazard is assumed constant",sep="")
	warning(outmes)
    }
    b <- as.vector(fit$b)
    names(b) <- varnames
    fit <- list(coefficients = b, var = -solve(fit$sd), iter = fit$nit, 
        loglik = fit$loglik)
    fit
}



survsplit <- function (data, cut, end, event, start, id = NULL, zero = 0, 
    episode = NULL, interval = NULL) 
{
    ntimes <- length(cut)
    n <- nrow(data)
    p <- ncol(data)
    if (length(interval) > 0) {
        ntimes <- ntimes - 1
        sttime <- c(rep(0, n), rep(cut[-length(cut)], each = n))
        endtime <- rep(cut, each = n)
    }
    else {
        endtime <- rep(c(cut, Inf), each = n)
        sttime <- c(rep(0, n), rep(cut, each = n))
    }
    newdata <- lapply(data, rep, ntimes + 1)
    eventtime <- newdata[[end]]
    if (start %in% names(data)) 
        starttime <- newdata[[start]]
    else starttime <- rep(zero, length = (ntimes + 1) * n)
    starttime <- pmax(sttime, starttime)
    epi <- rep(0:ntimes, each = n)
    if (length(interval) > 0) 
        status <- ifelse(eventtime <= endtime & eventtime >= 
            starttime, newdata[[event]], 0)
    else status <- ifelse(eventtime <= endtime & eventtime > 
        starttime, newdata[[event]], 0)
    endtime <- pmin(endtime, eventtime)
    if (length(interval) > 0) 
        drop <- (starttime > endtime) | (starttime == endtime & 
            status == 0)
    else drop <- starttime >= endtime
    newdata <- do.call("data.frame", newdata)
    newdata <- newdata[!drop, ]
    newdata[, start] <- starttime[!drop]
    newdata[, end] <- endtime[!drop]
    newdata[, event] <- status[!drop]
    if (!is.null(id)) 
        newdata[, id] <- rep(rownames(data), ntimes + 1)[!drop]
    fu <- NULL
    if (length(interval) > 2) {
        for (it in 1:length(interval[-1])) {
            drop1 <- sum(!drop[1:(it * n - n)])
            drop2 <- sum(!drop[(it * n - n + 1):(it * n)])
            drop3 <- sum(!drop[(it * n + 1):(length(interval[-1]) * 
                n)])
            if (it == 1) 
                fu <- cbind(fu, c(rep(1, drop2), rep(0, drop3)))
            else if (it == length(interval[-1])) 
                fu <- cbind(fu, c(rep(0, drop1), rep(1, drop2)))
            else fu <- cbind(fu, c(rep(0, drop1), rep(1, drop2), 
                rep(0, drop3)))
        }
        fu <- as.data.frame(fu)
        names(fu) <- c(paste("fu [", interval[-length(interval)], 
            ",", interval[-1], ")", sep = ""))
        newdata <- cbind(newdata, fu)
    }
    else if (length(interval) == 2) {
        fu <- rep(1, sum(!drop))
        newdata <- cbind(newdata, fu)
        names(newdata)[ncol(newdata)] <- paste("fu [", interval[1], 
            ",", interval[2], "]", sep = "")
    }
    if (!is.null(episode)) 
        newdata[, episode] <- epi[!drop]
    newdata
}


glmxp <- function (rform, data, interval, method, control) 
{
    if (rform$m == 1) 
        g <- as.integer(as.factor(rform$X[[1]]))
    else if (rform$m > 1) {
        gvar <- NULL
        for (i in 1:rform$m) {
            gvar <- append(gvar, rform$X[i])
        }
        tabgr <- as.data.frame(table(gvar))
        tabgr <- tabgr[, 1:rform$m]
        n.groups <- dim(tabgr)[1]
        mat <- do.call("data.frame", gvar)
        names(mat) <- names(tabgr)
        tabgr <- cbind(tabgr, g = as.numeric(row.names(tabgr)))
        mat <- cbind(mat, id = 1:rform$n)
        c <- merge(tabgr, mat)
        g <- c[order(c$id), rform$m + 1]
    }
    else g <- rep(1, rform$n)
    vg <- function(X) {
        n <- dim(X)[1]
        w <- sum((X$event == 0) & (X$fin == 1) & (X$y != 1))
        nd <- sum((X$event == 1) & (X$fin == 1))
        ps <- srvxp.fit(X[, 4:(nfk + 3)], rep(t.int, n), rform$ratetable)
        ld <- n - w/2
        lny <- log(sum(X$y))
        k <- t.int/365.24
        dstar <- sum(-log(ps)/k * X$y)
        ps <- mean(ps)
        if (rform$m == 0) 
            data.rest <- X[1, 7 + nfk + rform$m, drop = FALSE]
        else data.rest <- X[1, c((3 + nfk + 1):(3 + nfk + rform$m), 
            7 + nfk + rform$m)]
        cbind(nd = nd, ld = ld, ps = ps, lny = lny, dstar = dstar, 
            k = k, data.rest)
    }
    nint <- length(interval)
    if (nint < 2) 
        stop("Illegal interval value")
    meje <- interval
    my.fun <- function(x) {
        if (x > 1) {
            x.t <- rep(1, floor(x))
            if (x - floor(x) > 0) 
                x.t <- c(x.t, x - floor(x))
            x.t
        }
        else x
    }
    int <- apply(matrix(diff(interval), ncol = 1), 1, my.fun)
    if (is.list(int)) 
        int <- c(0, cumsum(do.call("c", int)))
    else int <- c(0, cumsum(int))
    int <- int * 365.24
    nint <- length(int)
    X <- cbind(rform$data, grupa = g)
    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    Z <- X[X$start >= int[2], ]
    nz <- dim(Z)[1]
    Z$fin <- rep(0, nz)
    Z$event <- rep(0, nz)
    Z$fu <- rep(0, nz)
    Z$y <- rep(0, nz)
    Z$origstart <- Z$start
    Z$xind <- rep(0, nz)
    if (nrow(Z) > 0) 
        Z[, 4:(nfk + 3)] <- Z[, 4:(nfk + 3)] + matrix(Z$start, 
            ncol = nfk, byrow = FALSE, nrow = nrow(Z)) * matrix(fk, 
            ncol = nfk, byrow = TRUE, nrow = nrow(Z))
    X <- X[X$start < int[2], ]
    X$fin <- (X$Y <= int[2])
    X$event <- X$fin * X$stat
    ford <- eval(substitute(paste("[", a, ",", b, "]", sep = ""), 
        list(a = meje[1], b = meje[2])))
    X$fu <- rep(ford, rform$n - nz)
    t.int <- int[2] - int[1]
    X$y <- (pmin(X$Y, int[2]) - X$start)/365.24
    X$origstart <- X$start
    X$xind <- rep(1, nrow(X))
    gr1 <- by(X, X$grupa, vg)
    grm1 <- do.call("rbind", gr1)
    X <- X[X$fin == 0, ]
    X$start <- rep(int[2], dim(X)[1])
    X <- rbind(X, Z[Z$start < int[3], ])
    Z <- Z[Z$start >= int[3], ]
    temp <- 0
    if (nint > 2) {
        for (i in 3:nint) {
            ni <- dim(X)[1]
            if (ni == 0) {
                temp <- 1
                break
            }
            X$fin <- X$Y <= int[i]
            X$event <- X$fin * X$stat
            l <- sum(int[i - 1] >= meje * 365.24)
            if(l==1)
		ftemp <- eval(substitute(paste("[", a, ",", b, "]", sep = ""),
			 list(a = meje[l], b = meje[l + 1])))

	    else
		ftemp <- eval(substitute(paste("(", a, ",", b, "]", sep = ""),
			 list(a = meje[l], b = meje[l + 1])))
            ford <- c(ford, ftemp)
            X$fu <- rep(ford[i - 1], ni)
            t.int <- int[i] - int[i - 1]
            index <- X$origstart < int[i - 1]
            index1 <- as.logical(X$xind)
            if (sum(index) > 0) 
                X[index, 4:(nfk + 3)] <- X[index, 4:(nfk + 3)] + 
                  matrix(fk * t.int, ncol = nfk, byrow = TRUE, 
                    nrow = sum(index))
            X$xind <- rep(1, nrow(X))
            X$y <- (pmin(X$Y, int[i]) - X$start)/365.24
            gr1 <- by(X, X$grupa, vg)
            grm1 <- rbind(grm1, do.call("rbind", gr1))
            X <- X[X$fin == 0, ]
            X$start <- rep(int[i], dim(X)[1])
            if (i == nint) 
                break
            X <- rbind(X, Z[Z$start < int[i + 1], ])
            X <- X[X$start != X$Y, ]
            Z <- Z[Z$start >= int[i + 1], ]
        }
        l <- sum(int[i - temp] > meje * 365.24)
        interval <- meje[1:(l + 1)]
    }
    else interval <- meje[1:2]
    grm1$fu <- factor(grm1$fu, levels = unique(ford))
    if (method == "glm.bin") {
        ht <- binomial(link = cloglog)
        ht$link <- "Hakulinen-Tenkanen relative survival model"
        ht$linkfun <- function(mu) log(-log((1 - mu)/ps))
        ht$linkinv <- function(eta) 1 - exp(-exp(eta)) * ps
        ht$mu.eta <- function(eta) exp(eta) * exp(-exp(eta)) * 
            ps
        ps <- grm1$ps
        assign(".ps", grm1$ps, env = .GlobalEnv)
        ht$initialize <- expression({
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
            mustart <- (n * y + 0.01)/(n + 0.02)
            mustart[(1 - mustart)/.ps >= 1] <- .ps[(1 - mustart)/.ps >= 
                1] * 0.9
        })
        if (any(grm1$ld - grm1$nd > grm1$ps * grm1$ld)) {
            n <- sum(grm1$ld - grm1$nd > grm1$ps * grm1$ld)
            g <- dim(grm1)[1]
            warnme <- paste("Observed number of deaths is smaller than the expected in ", 
                n, "/", g, " groups of patients", sep = "")
        }
        else warnme <- ""
        
        if (length(interval) == 2 & rform$m == 0) 
            stop("No groups can be formed")
        if (length(interval) == 1 | length(table(grm1$fu)) == 
            1) 
            grm1$fu <- as.integer(grm1$fu)
        if (!length(rform$X)) 
            local.ht <- glm(cbind(nd, ld - nd) ~ -1 + fu + offset(log(k)), 
                data = grm1, family = ht)
        else {
            xmat <- as.matrix(grm1[, 7:(ncol(grm1) - 1)])
            local.ht <- glm(cbind(nd, ld - nd) ~ -1 + xmat + 
                fu + offset(log(k)), data = grm1, family = ht)
        }
        names(local.ht[[1]]) <- c(names(rform$X), paste("fu", 
            levels(grm1$fu)))
    }
    else if (method == "glm.poi") {
        pot <- poisson()
        pot$link <- "glm relative survival model with Poisson error"
        pot$linkfun <- function(mu) log(mu - dstar)
        pot$linkinv <- function(eta) dstar + exp(eta)
        assign(".dstar", grm1$dstar, env = .GlobalEnv)
        if (any(grm1$nd - grm1$dstar < 0)) {
            pot$initialize <- expression({
                if (any(y < 0)) stop(paste("Negative values not allowed for", 
                  "the Poisson family"))
                n <- rep.int(1, nobs)
                mustart <- pmax(y, .dstar) + 0.1
            })
        }
        if (any(grm1$nd - grm1$dstar < 0)) {
            n <- sum(grm1$nd - grm1$dstar < 0)
            g <- dim(grm1)[1]
            warnme <- paste("Observed number of deaths is smaller than the expected in ", 
                n, "/", g, " groups of patients", sep = "")
        }
        else warnme <- ""
        dstar <- grm1$dstar
        if (length(interval) == 2 & rform$m == 0) 
            stop("No groups can be formed")
        if (length(interval) == 1 | length(table(grm1$fu)) == 
            1) 
            grm1$fu <- as.integer(grm1$fu)
        if (!length(rform$X)) 
            local.ht <- glm(nd ~ -1 + fu, data = grm1, family = pot, 
                offset = lny)
        else {
            xmat <- as.matrix(grm1[, 7:(ncol(grm1) - 1)])
            local.ht <- glm(nd ~ -1 + xmat + fu, data = grm1, 
                family = pot, offset = lny)
        }
        names(local.ht[[1]]) <- c(names(rform$X), paste("fu", 
            levels(grm1$fu)))
    }
    else stop(paste("Method '", method, "' not a valid method", 
        sep = ""))
    class(local.ht) <- c("rsadd", class(local.ht))
    local.ht$warnme <- warnme
    local.ht$int <- interval
    local.ht$groups <- local.ht$data
    return(local.ht)
}

residuals.rsadd <- function (object, type = "schoenfeld", ...) 
{
    data <- object$data[order(object$data$Y), ]
    ratetable <- object$ratetable
    beta <- object$coef
    start <- data[, 1]
    stop <- data[, 2]
    event <- data[, 3]
    fk <- (attributes(ratetable)$factor != 1)
    nfk <- length(fk)
    n <- nrow(data)
    scale <- 1
    if (object$method == "EM") 
        scale <- 365.24
    m <- ncol(data)
    rem <- m - nfk - 3
    interval <- object$int
    int <- ceiling(max(interval))
    R <- data[, 4:(nfk + 3)]
    lp <- matrix(-log(srvxp.fit(as.matrix(R), rep(365.24, n), 
        object$ratetable))/scale, ncol = 1)
    fu <- NULL
    if (object$method == "EM") {
        death.time <- stop[event == 1]
        for (it in 1:int) {
            fu <- as.data.frame(cbind(fu, as.numeric(death.time/365.24 < 
                it & (death.time/365.24) >= (it - 1))))
        }
        if(length(death.time)!=length(unique(death.time))){
        	utimes <- which(c(1,diff(death.time))!=0)
        	razteg <- function(x){
		# x is a 0/1 vector, the output is a vector of length sum(x), with the corresponding rep numbers
			n <- length(x)
			repu <- rep(1,n)
			repu[x==1] <- 0
			repu <- rev(cumsum(rev(repu)))
			repu <- repu[x==1]
			repu <- -diff(c(repu,0))+1
			if(sum(repu)!=n)repu <- c(n-sum(repu),repu)  #ce je prvi cas censoring, bo treba se kej narest??
			repu
		}
		rutd <- rep(0,length(death.time))
		rutd[utimes] <- 1
		rutd <- razteg(rutd)				#from  unique event times to event times
	}
	else rutd <- rep(1,length(death.time))
        lambda0 <- rep(object$lambda0,rutd)
    }
    else {
        pon <- NULL
        for (i in 1:(length(interval) - 1)) {
            width <- ceiling(interval[i + 1]) - floor(interval[i])
            lo <- interval[i]
            hi <- min(interval[i + 1], floor(interval[i]) + 1)
            for (j in 1:width) {
                fu <- as.data.frame(cbind(fu, as.numeric(stop/365.24 < 
                  hi & stop/365.24 >= lo)))
                names(fu)[ncol(fu)] <- paste("fu", lo, "-", hi, 
                  sep = "")
                if (j == width) {
                  pon <- c(pon, sum(fu[event == 1, (ncol(fu) - 
                    width + 1):ncol(fu)]))
                  break()
                }
                else {
                  lo <- hi
                  hi <- min(interval[i + 1], floor(interval[i]) + 
                    1 + j)
                }
            }
        }
        m <- ncol(data)
        data <- cbind(data, fu)
        rem <- m - nfk - 3
        lambda0 <- rep(exp(beta[rem + 1:(length(interval) - 1)]), 
            pon)
        fu <- fu[event == 1, , drop = FALSE]
        beta <- beta[1:rem]
    }
    if (int >= 2) {
        for (j in 2:int) {
            R <- R + matrix(fk * 365.24, ncol = ncol(R), byrow = TRUE, 
                nrow = n)
            xx <- srvxp.fit(R, rep(365.24, n), object$ratetable)
            lp <- cbind(lp, -log(xx)/scale)
        }
    }
    z <- as.matrix(data[, (4 + nfk):m])
    out <- resid.com(start, stop, event, z, beta, lp, lambda0, 
        fu, n, rem, int, type)
    out
}

resid.com <- function (start, stop, event, z, beta, lp, lambda0, fup, n, rem, 
    int, type) 
{
    le <- exp(z %*% beta)
    olp <- if (int > 1) 
        apply(lp[n:1, ], 2, cumsum)[n:1, ]
    else matrix(cumsum(lp[n:1])[n:1], ncol = 1)
    ole <- cumsum(le[n:1])[n:1]
    lp.st <- lp[order(start), , drop = FALSE]
    le.st <- le[order(start), , drop = FALSE]
    starter <- sort(start)
    starter1 <- c(starter[1], starter[-length(starter)])
    index <- c(TRUE, (starter != starter1)[-1])
    starter <- starter[index]
    val1 <- apply(matrix(starter, ncol = 1), 1, function(x, Y) sum(x >= 
        Y), stop)
    val1 <- c(val1[1], diff(val1), length(stop) - val1[length(val1)])
    olp.st <- (apply(lp.st[n:1, , drop = FALSE], 2, cumsum)[n:1, 
        , drop = FALSE])[index, , drop = FALSE]
    olp.st <- apply(olp.st, 2, function(x) rep(c(x, 0), val1))
    olp <- olp - olp.st
    olp <- olp[event == 1, ]
    olp <- apply(fup * olp, 1, sum)
    ole.st <- cumsum(le.st[n:1])[n:1][index]
    ole.st <- rep(c(ole.st, 0), val1)
    ole <- ole - ole.st
    ole <- ole[event == 1] * lambda0
    s0 <- ole + olp
    sc <- NULL
    zb <- NULL
    kzb <- NULL
    f1 <- function(x) rep(mean(x), length(x))
    f2 <- function(x) apply(x, 2, f1)
    f3 <- function(x) apply(x, 1:2, f1)
    ties <- length(unique(stop[event == 1])) != length(stop[event == 
        1])
    for (k in 1:rem) {
        zlp <- apply((z[, k] * lp)[n:1, , drop = FALSE], 2, cumsum)[n:1, 
            , drop = FALSE]
        zlp.st <- (apply((z[, k] * lp.st)[n:1, , drop = FALSE], 
            2, cumsum)[n:1, , drop = FALSE])[index, , drop = FALSE]
        zlp.st <- apply(zlp.st, 2, function(x) rep(c(x, 0), val1))
        zlp <- zlp - zlp.st
        zlp <- zlp[event == 1, , drop = FALSE]
        zlp <- apply(fup * zlp, 1, sum)
        zle <- cumsum((z[, k] * le)[n:1])[n:1]
        zle.st <- cumsum((z[, k] * le.st)[n:1])[n:1][index]
        zle.st <- rep(c(zle.st, 0), val1)
        zle <- zle - zle.st
        zle <- zle[event == 1]
        zle <- zle * lambda0
        s1 <- zle + zlp
        zb <- cbind(zb, s1/s0)
        kzb <- cbind(kzb, zle/s0)
    }
    s1ties <- cbind(zb, kzb)
    if (ties) {
        s1ties <- by(s1ties, stop[event == 1], f2)
        s1ties <- do.call("rbind", s1ties)
    }
    zb <- s1ties[, 1:rem, drop = FALSE]
    kzb <- s1ties[, -(1:rem), drop = FALSE]
    sc <- z[event == 1, , drop = FALSE] - zb
    row.names(sc) <- stop[event == 1]
    out.temp <- function(x) outer(x, x, fun = "*")
    krez <- rez <- array(matrix(NA, ncol = rem, nrow = rem), 
        dim = c(rem, rem, sum(event == 1)))
    for (a in 1:rem) {
        for (b in a:rem) {
            zzlp <- apply((z[, a] * z[, b] * lp)[n:1, , drop = FALSE], 
                2, cumsum)[n:1, , drop = FALSE]
            zzlp.st <- (apply((z[, a] * z[, b] * lp.st)[n:1, 
                , drop = FALSE], 2, cumsum)[n:1, , drop = FALSE])[index, 
                , drop = FALSE]
            zzlp.st <- apply(zzlp.st, 2, function(x) rep(c(x, 
                0), val1))
            zzlp <- zzlp - zzlp.st
            zzlp <- zzlp[event == 1, , drop = FALSE]
            zzlp <- apply(fup * zzlp, 1, sum)
            zzle <- cumsum((z[, a] * z[, b] * le)[n:1])[n:1]
            zzle.st <- cumsum((z[, a] * z[, b] * le.st)[n:1])[n:1][index]
            zzle.st <- rep(c(zzle.st, 0), val1)
            zzle <- zzle - zzle.st
            zzle <- zzle[event == 1]
            zzle <- zzle * lambda0
            s2 <- zzlp + zzle
            s20 <- s2/s0
            ks20 <- zzle/s0
            s2ties <- cbind(s20, ks20)
            if (ties) {
                s2ties <- by(s2ties, stop[event == 1], f2)
                s2ties <- do.call("rbind", s2ties)
            }
            rez[a, b, ] <- rez[b, a, ] <- s2ties[, 1]
            krez[a, b, ] <- krez[b, a, ] <- s2ties[, 2]
        }
    }
    juhu <- apply(zb, 1, out.temp)
    if (is.null(dim(juhu))) 
        juhu1 <- array(data = matrix(juhu, ncol = a), dim = c(a, 
            a, length(zb[, 1])))
    else juhu1 <- array(data = apply(juhu, 2, matrix, ncol = a), 
        dim = c(a, a, length(zb[, 1])))
    varr <- rez - juhu1
    kjuhu <- apply(cbind(zb, kzb), 1, function(x) outer(x[1:rem], 
        x[-(1:rem)], fun = "*"))
    if (is.null(dim(kjuhu))) 
        kjuhu1 <- array(data = matrix(kjuhu, ncol = rem), dim = c(rem, 
            rem, length(zb[, 1])))
    else kjuhu1 <- array(data = apply(kjuhu, 2, matrix, ncol = rem), 
        dim = c(rem, rem, length(zb[, 1])))
    kvarr <- krez - kjuhu1
    for (i in 1:dim(varr)[1]) varr[i, i, which(varr[i, i, ] < 
        0)] <- 0
    for (i in 1:dim(kvarr)[1]) kvarr[i, i, which(kvarr[i, i, 
        ] < 0)] <- 0
    varr1 <- apply(varr, 1:2, sum)
    kvarr1 <- apply(kvarr, 1:2, sum)
    if (type == "schoenfeld") 
        out <- list(res = sc, varr1 = varr1, varr = varr, kvarr = kvarr, 
            kvarr1 = kvarr1)
    out
}




rs.br <- function (fit, sc, rho = 0, test = "max", global = TRUE) 
{
    test <- match.arg(test,c("max","cvm"))	
    if (inherits(fit, "rsadd")) {
        if (missing(sc)) 
            sc <- resid(fit, "schoenfeld")
        sresid <- sc$res
        varr <- sc$varr
        sresid <- as.matrix(sresid)
    }
    else {
        coef <- fit$coef
        options(warn = -1)
        sc <- coxph.detail(fit)
        options(warn = 0)
        sresid <- sc$score
        varr <- sc$imat
        if (is.null(dim(varr))) 
            varr <- array(varr, dim = c(1, 1, length(varr)))
        sresid <- as.matrix(sresid)
    }
    if (inherits(fit, "coxph")) {
	if(is.null(fit$data)){
		temp <- fit$y
		class(temp) <- "matrix"
		if(ncol(fit$y)==2)temp <- data.frame(rep(0,nrow(fit$y)),temp)
		if(is.null(fit$x))stop("The coxph model should be called with x=TRUE argument")
		fit$data <- data.frame(temp,fit$x)
		names(fit$data)[1:3] <- c("start","Y","stat")
	}
    }
    data <- fit$data[order(fit$data$Y), ]
    time <- data$Y[data$stat == 1]
    ties <- (length(unique(time)) != length(time))
    keep <- 1:(ncol(sresid))
    options(warn = -1)
    scaled <- NULL
     varnova <- NULL
    if (ncol(sresid) == 1) {
        varr <- varr[1, 1, ]
        scaled <- sresid/sqrt(varr)
    }
    else { for (i in 1:ncol(sresid)) varnova <- cbind(varnova,varr[i,i,])
    	   scaled <- sresid/sqrt(varnova)
    	  }

    options(warn = 0)
    nvar <- ncol(sresid)
    survfit <- getFromNamespace("survfit", "survival")
    temp <- survfit(fit$y~1, type = "kaplan-meier")
    n.risk <- temp$n.risk
    n.time <- temp$time
    if (temp$type == "right") {
        cji <- matrix(fit$y, ncol = 2)
        n.risk <- n.risk[match(cji[cji[, 2] == 1, 1], n.time)]
    }
    else {
        cji <- matrix(fit$y, ncol = 3)
        n.risk <- n.risk[match(cji[cji[, 3] == 1, 2], n.time)]
    }
    n.risk <- sort(n.risk, decreasing = TRUE)
    varnames <- names(fit$coef)[keep]
    u2 <- function(bb) {
        n <- length(bb)
        1/n * (sum(bb^2) - sum(bb)^2/n)
    }
    wc <- function(x, k = 1000) {
        a <- 1
        for (i in 1:k) a <- a + 2 * (-1)^i * exp(-2 * i^2 * pi^2 * 
            x)
        a
    }
    brp <- function(x, n = 1000) {
        a <- 1
        for (i in 1:n) a <- a - 2 * (-1)^(i - 1) * exp(-2 * i^2 * 
            x^2)
        a
    }
    global <- as.numeric(global & ncol(sresid) > 1)
    table <- NULL
    bbt <- as.list(1:(nvar + global))
    for (i in 1:nvar) {
        if (nvar != 1) 
            usable <- which(varr[i, i, ] > 1e-12)
        else usable <- which(varr > 1e-12)
        w <- (n.risk[usable])^rho
        w <- w/sum(w)
        if (nvar != 1) {
            sci <- scaled[usable, i]
        }
        else sci <- scaled[usable]
        if (ties) {
            if (inherits(fit, "rsadd")) {
                sci <- as.vector(by(sci, time[usable], function(x) sum(x)/sqrt(length(x))))
                w <- as.vector(by(w, time[usable], sum))
            }
            else {
                w <- w * as.vector(table(time))[usable]
                w <- w/sum(w)
            }
        }
        sci <- sci * sqrt(w)
        timescale <- cumsum(w)
        bm <- cumsum(sci)
        bb <- bm - timescale * bm[length(bm)]
        if (test == "max") 
            table <- rbind(table, c(max(abs(bb)), 1 - brp(max(abs(bb)))))
        else if (test == "cvm") 
            table <- rbind(table, c(u2(bb), 1 - wc(u2(bb))))
        bbt[[i]] <- cbind(timescale, bb)
    }
    if (inherits(fit, "rsadd")) {
       beta <- fit$coef[1:(length(fit$coef) - length(fit$int) +  1)]
    }
    else beta <- fit$coef
    if (global) {
        qform <- function(matrix, vector) t(vector) %*% matrix %*% 
            vector
        diagonal <- apply(varr, 3, diag)
        sumdiag <- apply(diagonal, 2, sum)
        usable <- which(sumdiag > 1e-12)
        score <- t(beta) %*% t(sresid[usable, ])
        varr <- varr[, , usable]
        qf <- apply(varr, 3, qform, vector = beta)
        w <- (n.risk[usable])^rho
        w <- w/sum(w)
        sci <- score/(qf)^0.5
        if (ties) {
            if (inherits(fit, "rsadd")) {
                sci <- as.vector(by(t(sci), time[usable], function(x) sum(x)/sqrt(length(x))))
                w <- as.vector(by(w, time[usable], sum))
            }
            else {
                w <- w * as.vector(table(time))
                w <- w/sum(w)
            }
        }
        sci <- sci * sqrt(w)
        timescale <- cumsum(w)
        bm <- cumsum(sci)
        bb <- bm - timescale * bm[length(bm)]
        if (test == "max") 
            table <- rbind(table, c(max(abs(bb)), 1 - brp(max(abs(bb)))))
        else if (test == "cvm") 
            table <- rbind(table, c(u2(bb), 1 - wc(u2(bb))))
        bbt[[nvar + 1]] <- cbind(timescale, bb)
        varnames <- c(varnames, "GLOBAL")
    }
    dimnames(table) <- list(varnames, c(test, "p"))
    out <- list(table = table, bbt = bbt, rho = rho)
    class(out) <- "rs.br"
    out
}

rs.zph <- function (fit, sc, transform = "identity", var.type = "sum") 
{
    if (inherits(fit, "rsadd")) {
            if (missing(sc)) 
                sc <- resid(fit, "schoenfeld")
            sresid <- sc$res
            varr <- sc$kvarr
            fvar <- solve(sc$kvarr1)
            sresid <- as.matrix(sresid)
        }
        else {
            coef <- fit$coef
            options(warn = -1)
            sc <- coxph.detail(fit)
            options(warn = 0)
            sresid <- as.matrix(resid(fit, "schoenfeld"))
            varr <- sc$imat
            fvar <- fit$var
    }
   data <- fit$data[order(fit$data$Y), ]
   time <- data$Y
   stat <- data$stat

   if (!inherits(fit, "rsadd")) {
           ties <- as.vector(table(time[stat==1]))
           if(is.null(dim(varr))) varr <- rep(varr/ties,ties)
            else{
    		    varr <- apply(varr,1:2,function(x)rep(x/ties,ties))
    		    varr <- aperm(varr,c(2,3,1))
    	    }
    }
    keep <- 1:(length(fit$coef) - length(fit$int) + 1)
    varnames <- names(fit$coef)[keep]
    nvar <- length(varnames)
    ndead <- length(sresid)/nvar
    if (inherits(fit, "rsadd")) 
        times <- time[stat == 1]
    else times <- sc$time
    if (is.character(transform)) {
        tname <- transform
        ttimes <- switch(transform, identity = times, rank = rank(times), 
            log = log(times), km = {
                fity <- Surv(time, stat)
                temp <- survfit.km(factor(rep(1, nrow(fity))), 
                  fity, se.fit = FALSE)
                t1 <- temp$surv[temp$n.event > 0]
                t2 <- temp$n.event[temp$n.event > 0]
                km <- rep(c(1, t1), c(t2, 0))
                if (is.null(attr(sresid, "strata"))) 
                  1 - km
                else (1 - km[sort.list(sort.list(times))])
            }, stop("Unrecognized transform"))
    }
    else {
        tname <- deparse(substitute(transform))
        ttimes <- transform(times)
    }
    if (var.type == "each") {
            invV <- apply(varr, 3, function(x) try(solve(x), silent = TRUE))
            if (length(invV) == length(varr)){ 
                if(!is.numeric(invV)){
                	usable <- rep(FALSE, dim(varr)[3])
                	options(warn=-1)
                	invV <- as.numeric(invV)
                	usable[1:(min(which(is.na(invV)))-1)] <- TRUE
                	invV <- invV[usable]
                	sresid <- sresid[usable,,drop=FALSE]
                	options(warn=0)
                }
                else usable <- rep(TRUE, dim(varr)[3])
            }
            else {
                usable <- unlist(lapply(invV, is.matrix))
                if (!any(usable)) 
                    stop("All the matrices are singular")
                invV <- invV[usable]
                sresid <- sresid[usable, , drop = FALSE]
        }
        di1 <- dim(varr)[1]
        di3 <- sum(usable)
        u <- array(data = matrix(unlist(invV), ncol = di1), dim = c(di1, 
            di1, di3))
        uv <- cbind(matrix(u, ncol = di1, byrow = TRUE), as.vector(t(sresid)))
        uv <- array(as.vector(t(uv)), dim = c(di1 + 1, di1, di3))
        r2 <- t(apply(uv, 3, function(x) x[1:di1, ] %*% x[di1 + 
            1, ]))
        r2 <- matrix(r2, ncol = di1)
        whr2 <-  apply(r2<100,1,function(x)!any(x==FALSE))
        usable <- as.logical(usable*whr2)
        r2 <- r2[usable,,drop=FALSE]
        u <- u[,,usable]
        dimnames(r2) <- list(times[usable], varnames)
        temp <- list(x = ttimes[usable], y = r2 + outer(rep(1, 
            sum(usable)), fit$coef[keep]), var = u, call = call, 
            transform = tname)
    }
    else if (var.type == "sum") {
       xx <- ttimes - mean(ttimes)
       r2 <- t(fvar %*% t(sresid) * ndead)
       r2 <- as.matrix(r2)
       dimnames(r2) <- list(times, varnames)
       temp <- list(x = ttimes, y = r2 + outer(rep(1, ndead), 
       fit$coef[keep]), var = fvar, transform = tname)
    }
    else stop("Unknown 'var.type'")
    class(temp) <- "rs.zph"
    temp
}

plot.rs.zph <- function (x,resid = TRUE, df = 4, nsmo = 40, var, cex = 1,  add = FALSE, col = 1, 
    lty = 1, xlab, ylab, scale = 1, ...) 
{
    require(splines)
    xx <- x$x
    if(x$transform=="identity")xx <- xx/scale
    yy <- x$y
    d <- nrow(yy)
    df <- max(df)
    nvar <- ncol(yy)
    pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
    temp <- c(pred.x, xx)
    lmat <- ns(temp, df = df, intercept = TRUE)
    pmat <- lmat[1:nsmo, ]
    xmat <- lmat[-(1:nsmo), ]
    qmat <- qr(xmat)
    if (missing(ylab)) 
        ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
    if (missing(xlab)) 
        xlab <- "Time"
    if (missing(var)) 
        var <- 1:nvar
    else {
        if (is.character(var)) 
            var <- match(var, dimnames(yy)[[2]])
        if (any(is.na(var)) || max(var) > nvar || min(var) < 
            1) 
            stop("Invalid variable requested")
    }
    if (x$transform == "log") {
            xx <- exp(xx)
            pred.x <- exp(pred.x)
        }
    else if (x$transform != "identity") {
            xtime <- as.numeric(dimnames(yy)[[1]])/scale
            apr1 <- approx(xx, xtime, seq(min(xx), max(xx), length = 17)[2 * 
                (1:8)])
            temp <- signif(apr1$y, 2)
            apr2 <- approx(xtime, xx, temp)
            xaxisval <- apr2$y
            xaxislab <- rep("", 8)
            for (i in 1:8) xaxislab[i] <- format(temp[i])
    }
    for (i in var) {
        y <- yy[, i]
        yhat <- pmat %*% qr.coef(qmat, y)
        yr <- range(yhat, y)
        if (!add) {
		if (x$transform == "identity") 
		    plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i],...)
		else if (x$transform == "log") 
		    plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i],log = "x", ...)
		else {
		    plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i],axes = FALSE, ...)
		    axis(1, xaxisval, xaxislab)
		    axis(2)
		    box()
		}
        }
        if (resid) 
            points(xx, y, cex = cex, col = col)
        lines(pred.x, yhat, col = col, lty = lty)
    }
}

plot.rs.br <- function (x, var, ylim = c(-2, 2), xlab, ylab, ...) 
{
    bbt <- x$bbt
    par(ask = TRUE)
    if (missing(var)) 
        var <- 1:nrow(x$table)
    ychange <- FALSE
    if (missing(ylab)) 
        ylab <- paste("Brownian bridge for", row.names(x$table))
    else {
        if (length(ylab) == 1 & nrow(x$table) > 1) 
            ylab <- rep(ylab, nrow(x$table))
    }
    if (missing(xlab)) 
        xlab <- "Time"
    for (i in var) {
        timescale <- bbt[[i]][, 1]
        bb <- bbt[[i]][, 2]
        plot(c(0, timescale), c(0, bb), type = "l", ylim = ylim, 
            xlab = xlab, ylab = ylab[i], ...)
        abline(h = 1.36, col = 2)
        abline(h = 1.63, col = 2)
        abline(h = -1.36, col = 2)
        abline(h = -1.63, col = 2)
    }
    par(ask = FALSE)
}


Kernmatch <- function (t, tv, b, tD, nt4) 
{
    kmat <- NULL
    for (it in 1:(length(nt4) - 1)) {
        kmat1 <- (outer(t[(nt4[it] + 1):nt4[it + 1]], tv, "-")/b[it])
        kmat1 <- kmat1^(kmat1 >= 0)
        kmat <- rbind(kmat, pmax(1 - kmat1^2, 0) * (1.5/b[it]))
    }
    kmat
}

kernerleftch <- function (td, b, nt4) 
{
    n <- length(td)
    ttemp <- td[td >= b[1]]
    ntemp <- length(ttemp)
    if (ntemp == n) 
        nt4 <- c(0, nt4[-1])
    else {
        nfirst <- n - ntemp
        nt4 <- c(0, 1:nfirst, nt4[-1])
        b <- c(td[1:nfirst], b)
    }
    krn <- Kernmatch(td, td, b, max(td), nt4)
    krn
}


invtime <- function (y = 0.1, age = 23011, sex = "male", year = 9497, scale = 1, 
    ratetable = survexp.us, lower, upper) 
{
    if (!is.numeric(age)) 
        stop("\"age\" must be numeric", call. = FALSE)
    if (!is.numeric(y)) 
        stop("\"y\" must be numeric", call. = FALSE)
    if (!is.numeric(scale)) 
        stop("\"scale\" must be numeric", call. = FALSE)
    temp <- data.frame(age = age, sex = I(sex), year = year)
    if (missing(lower)) {
        if (!missing(upper)) 
            stop("Argument \"lower\" is missing, with no default", 
                call. = FALSE)
        nyears <- round((110 - age/365.24))
        tab <- data.frame(age = rep(age, nyears), sex = I(rep(sex, 
            nyears)), year = rep(year, nyears))
        vred <- 1 - survexp(c(0, 1:(nyears - 1)) * 365.24 ~ ratetable(age = age, 
            sex = sex, year = year), ratetable = ratetable, data = tab, 
            cohort = FALSE)
        place <- sum(vred <= y)
        if (place == 0) 
            lower <- 0
        else lower <- floor((place - 1) * 365.24 - place)
        upper <- ceiling(place * 365.24 + place)
    }
    else {
        if (missing(upper)) 
            stop("Argument \"upper\" is missing, with no default", 
                call. = FALSE)
        if (!is.integer(lower)) 
            lower <- floor(lower)
        if (!is.integer(upper)) 
            upper <- ceiling(upper)
        if (upper <= lower) 
            stop("'upper' must be higher than 'lower'", call. = FALSE)
    }
    lower <- max(0, lower)
    tab <- data.frame(age = rep(age, upper - lower + 1), sex = I(rep(sex, 
        upper - lower + 1)), year = rep(year, upper - lower + 
        1))
    vred <- 1 - survexp((lower:upper) ~ ratetable(age = age, 
        sex = sex, year = year), ratetable = ratetable, data = tab, 
        cohort = FALSE)
    place <- sum(vred <= y)
    if (place == 0) 
        warning(paste("The event happened on or before day", 
            lower), call. = FALSE)
    if (place == length(vred)) 
        warning(paste("The event happened on or after day", upper), 
            call. = FALSE)
    t <- (place + lower - 1)/scale
    age <- round(age/365.24, 0.01)
    return(list(age, sex, year, Y = y, T = t))
}




rsmul <- function (formula = formula(data), data = parent.frame(), ratetable = survexp.us, 
    int, na.action, init, method = "mul", control, ...) 
{
    require(survival)
    rform <- rformulate(formula, data, ratetable, na.action, 
        int)
    U <- rform$data
    if (missing(int)) 
	    int <- ceiling(max(rform$Y/365.24))
    if(length(int)!=1)int <- max(int)
    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    if (method == "mul") {
        U <- survsplit(U, cut = (1:int) * 365.24, end = "Y", 
            event = "stat", start = "start", episode = "epi")
        fk <- (attributes(rform$ratetable)$factor != 1)
        nfk <- length(fk)
        U[, 4:(nfk + 3)] <- U[, 4:(nfk + 3)] + 365.24 * (U$epi) %*% 
            t(fk)
        nsk <- dim(U)[1]
        xx <- srvxp.fit(U[, 4:(nfk + 3)], rep(365.24, nsk), rform$ratetable)
        lambda <- -log(xx)/365.24
    }
    else if (method == "mul1") {
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
        intr <- apply(U[, 1:(3 + nfk)], 1, my.fun, attcut, nfk, 
            fk)
        dolg <- unlist(lapply(intr, length))
        newdata <- lapply(U, rep, dolg)
        stoptime <- unlist(intr)
        starttime <- c(-1, stoptime[-length(stoptime)])
        first <- newdata$id != c(-1, newdata$id[-length(newdata$id)])
        starttime[first] <- newdata$start[first]
        last <- newdata$id != c(newdata$id[-1], -1)
        event <- rep(0, length(newdata$id))
        event[last] <- newdata$stat[last]
        U <- do.call("data.frame", newdata)
        U$start <- starttime
        U$Y <- stoptime
        U$stat <- event
        U[, 4:(nfk + 3)] <- U[, 4:(nfk + 3)] + (U$start) %*% 
            t(fk)
        nsk <- dim(U)[1]
        xx <- srvxp.fit(U[, 4:(nfk + 3)], rep(1, nsk), rform$ratetable)
        lambda <- -log(xx)/1
    }
    else stop("'method' must be one of 'mul' or 'mul1'")
    U$lambda <- log(lambda)
    if (rform$m == 0) 
        fit <- coxph(Surv(start, Y, stat) ~ 1 + offset(lambda), 
            data = U, init = init, control = control, x = TRUE, 
            ...)
    else {
        xmat <- as.matrix(U[, (3 + nfk + 1):(ncol(U) - 2)])
        fit <- coxph(Surv(start, Y, stat) ~ xmat + offset(lambda), 
            data = U, init = init, control = control, x = TRUE, 
            ...)
        names(fit[[1]]) <- names(U)[(3 + nfk + 1):(ncol(U) - 
            2)]
    }
    class(fit) <- c("rsmul",class(fit))
    fit$lambda <- log(lambda)
    fit$data <- rform$data
    fit$call <- match.call()
    fit$int <- int
    if (length(rform$na.action)) 
        fit$na.action <- rform$na.action
    fit
}

rstrans <- function (formula = formula(data), data = parent.frame(), ratetable = survexp.us, 
    int, na.action, init, control, ...) 
{
    rform <- rformulate(formula, data, ratetable, na.action, 
        int)
    if (missing(int)) 
	    int <- ceiling(max(rform$Y/365.24))
    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    if (rform$type == "counting") {
        start <- 1 - srvxp.fit(rform$R, rform$start, rform$ratetable)
    }
    else start <- rep(0, rform$n)
    stop <- 1 - srvxp.fit(rform$R, rform$Y, rform$ratetable)
     if(any(stop==0&rform$Y!=0))stop[stop==0&rform$Y!=0] <- .Machine$double.eps
     if(length(int)!=1)int <- max(int)
    data <- rform$data
    stat <- rform$status
    if (rform$m == 0) {
        if (rform$type == "counting") 
            fit <- coxph(Surv(start, stop, stat) ~ 1,
                init = init, control = control, x = TRUE, ...)
        else fit <- coxph(Surv(stop, stat) ~ 1, 
            init = init, control = control, x = TRUE, ...)
    }
    else {
        xmat <- as.matrix(data[, (4 + nfk):ncol(data)])
        fit <- coxph(Surv(start, stop, stat) ~ xmat,  
            init = init, control = control, x = TRUE, ...)
        names(fit[[1]]) <- names(rform$X)
    }
    fit$call <- match.call()
    if (length(rform$na.action)) 
        fit$na.action <- rform$na.action
    data$start <- start
    data$Y <- stop
    fit$data <- data
    fit$int <- int
    return(fit)
}
transrate <- function (men, women, yearlim, int.length = 1) 
{
    if (any(dim(men) != dim(women))) 
        stop("The men and women matrices must be of the same size. \n In case of missing values at the end carry the last value forward")
    if ((yearlim[2] - yearlim[1])/int.length + 1 != dim(men)[2]) 
        stop("'yearlim' cannot be divided into intervals of equal length")
    if (!is.matrix(men) | !is.matrix(women)) 
        stop("input tables must be of class matrix")
    dimi <- dim(men)
    temp <- array(c(men, women), dim = c(dimi, 2))
    temp <- -log(temp)/365.24
    temp <- aperm(temp, c(1, 3, 2))
    cp <- as.date(apply(matrix(yearlim[1] + int.length * (0:(dimi[2] - 
        1)), ncol = 1), 1, function(x) {
        paste("1jan", x, sep = "")
    }))
    attributes(temp) <- list(dim = c(dimi[1], 2, dimi[2]), dimnames = list(as.character(0:(dimi[1] - 
        1)), c("male", "female"), as.character(yearlim[1] + int.length * 
        (0:(dimi[2] - 1)))), dimid = c("age", "sex", "year"), 
        factor = c(0, 1, 0),type=c(2,1,3), cutpoints = list((0:(dimi[1] - 1)) * 
            (365.24), NULL, cp), class = "ratetable")
    attributes(temp)$summary <- function (R) 
	{
		x <- c(format(round(min(R[, 1])/365.24, 1)), format(round(max(R[, 
		1])/365.24, 1)), sum(R[, 2] == 1), sum(R[, 2] == 2))
		x2 <- as.character(as.date(c(min(R[, 3]), max(R[, 3]))))
		paste("  age ranges from", x[1], "to", x[2], "years\n", " male:", 
		x[3], " female:", x[4], "\n", " date of entry from", 
		x2[1], "to", x2[2], "\n")
	}
    temp
}

transrate.hld <- function(file, cut.year,race){
	nfiles <- length(file)
	data <- NULL
	for(it in 1:nfiles){
		tdata <- read.table(file[it],sep=",",header=TRUE)
		if(!any(tdata$TypeLT==1)) stop("Currently only TypeLT 1 is implemented")
		names(tdata) <- gsub(".","",names(tdata),fixed=TRUE)
		tdata <- tdata[,c("Country","Year1","Year2","TypeLT","Sex","Age","AgeInt","qx")]
		tdata <- tdata[tdata$TypeLT==1&tdata$AgeInt==1,]
		if(!missing(race))tdata$race <- rep(race[it],nrow(tdata))
		data <- rbind(data,tdata)
	}
	if(length(unique(data$Country))>1)warning("The data belongs to different countries")
	data <- data[order(data$Year1,data$Age),]
	data$qx <- as.character(data$qx)
	options(warn = -1)
	data$qx[data$qx=="."] <- NA
	data$qx <- as.numeric(data$qx)
	options(warn = 0)
	if(missing(cut.year)){
		y1 <-  unique(data$Year1)
		y2 <-  unique(data$Year2)
		if(any(apply(cbind(y1[-1],y2[-length(y2)]),1,diff)!=-1))warning("Data is not given for all the cut.year between the minimum and the maximum, use argument 'cut.year'")
	}
	else
		y1 <- cut.year
	if(length(y1)!=length(unique(data$Year1)))stop("Length 'cut.year' must match the number of unique values of Year1")
	cp <- as.date(apply(matrix(y1,ncol=1),1,function(x){paste("1jan",x,sep="")}))
	dn2 <- as.character(y1)
	amax <- max(data$Age)
	a.fun <- function(data,amax){
		mdata <- data[data$Sex==1,]
		wdata <- data[data$Sex==2,]
		men <-NULL
		women <- NULL
		k <- sum(mdata$Age==0)
		mind <- c(which(mdata$Age[-nrow(mdata)] != mdata$Age[-1]-1),nrow(mdata))
		wind <- c(which(wdata$Age[-nrow(wdata)] != wdata$Age[-1]-1),nrow(wdata))
		mst <- wst <- 1
		for(it in 1:k){
			qx <- mdata[mst:mind[it],]$qx
			lqx <- length(qx)
			if(lqx!=amax+1){
				nmiss <- amax + 1 - lqx
				qx <- c(qx,rep(qx[lqx],nmiss))
			}
			naqx <- max(which(!is.na(qx)))
			if(naqx!=amax+1) qx[(naqx+1):(amax+1)] <- qx[naqx]
			men <- cbind(men,qx)
			mst <- mind[it]+1 
			qx <- wdata[wst:wind[it],]$qx
			lqx <- length(qx)
			if(lqx!=amax+1){
				nmiss <- amax + 1 - lqx
				qx <- c(qx,rep(qx[lqx],nmiss))
			}
			naqx <- max(which(!is.na(qx)))
			if(naqx!=amax+1) qx[(naqx+1):(amax+1)] <- qx[naqx]
			women <- cbind(women,qx)
			wst <- wind[it]+1 
		}
		men<- -log(1-men)/365.241
		women<- -log(1-women)/365.241
		dims <- c(dim(men),2)
		array(c(men,women),dim=dims)
	}
	if(missing(race)){
		out <- a.fun(data,amax)
		dims <- dim(out)
		attributes(out)<-list(
			dim=dims,		
			dimnames=list(as.character(0:amax),as.character(y1),c("male","female")),	
			dimid=c("age","year","sex"),
			factor=c(0,0,1),type=c(2,3,1),
			cutpoints=list((0:amax)*(365.24),cp,NULL),
			class="ratetable"
		)
		
	}
	else{
		race.val <- unique(race)
		if(length(race)!=length(file))stop("Length of 'race' must match the number of files")
		for(it in 1:length(race.val)){
			if(it==1){
				out <- a.fun(data[data$race==race.val[it],],amax)
				dims <- dim(out)
				out <- array(out,dim=c(dims,1))
			}
			else{
				out1 <- array(a.fun(data[data$race==race.val[it],],amax),dim=c(dims,1))
				out <- array(c(out,out1),dim=c(dims,it))
			}
		}
		attributes(out)<-list(
			dim=c(dims,it),		
			dimnames=list(as.character(0:amax),as.character(y1),c("male","female"),race.val),	
			dimid=c("age","year","sex","race"),
			factor=c(0,0,1,1),type=c(2,3,1,1),
			cutpoints=list((0:amax)*(365.24),cp,NULL,NULL),
			class="ratetable"
		)
		attributes(out)$summary <- function (R) 
		{
			x <- c(format(round(min(R[, 1])/365.24, 1)), format(round(max(R[, 
			1])/365.24, 1)), sum(R[, 2] == 1), sum(R[, 2] == 2))
			x2 <- as.character(as.date(c(min(R[, 3]), max(R[, 3]))))
			paste("  age ranges from", x[1], "to", x[2], "years\n", " male:", 
			x[3], " female:", x[4], "\n", " date of entry from", 
			x2[1], "to", x2[2], "\n")
		}
	}
	out
}

transrate.hmd <- function(male,female){
	nfiles <- 2
	men <- read.table(male,sep="",header=TRUE)
	men <- men[,c("Year","Age","qx")]
	y1 <- sort(unique(men$Year))
	ndata <- nrow(men)/111
	if(round(ndata)!=ndata)stop("Each year must contain ages from 0 to 110")
	men <- matrix(men$qx, ncol=ndata)
	men <- matrix(as.numeric(men),ncol=ndata)
	women <- read.table(female,sep="",header=TRUE)
	women <- women[,"qx"]
	if(length(women)!=length(men))stop("Number of rows in the table must be equal for both sexes")
	women <- matrix(women, ncol=ndata)
	women <- matrix(as.numeric(women),ncol=ndata)
		cp <- as.date(apply(matrix(y1,ncol=1),1,function(x){paste("1jan",x,sep="")}))
	dn2 <- as.character(y1)
	tfun <- function(vec){
		ind <- which(vec == 1 | is.na(vec))
		if(length(ind)>0)vec[min(ind):length(vec)] <- 0.999
		vec
	}
	men <- apply(men,2,tfun)
	women <- apply(women,2,tfun)
	men<- -log(1-men)/365.241
	women<- -log(1-women)/365.241
	nr <- nrow(men)-1
	dims <- c(dim(men),2)
	out <- array(c(men,women),dim=dims)
	attributes(out)<-list(
		dim=dims,
		dimnames=list(as.character(0:nr),as.character(y1),c("male","female")),	
		dimid=c("age","year","sex"),
		factor=c(0,0,1),type=c(2,3,1),
		cutpoints=list((0:nr)*(365.24),cp,NULL),
		class="ratetable"
	)
	attributes(out)$summary <- function (R) 
	{
		x <- c(format(round(min(R[, 1])/365.24, 1)), format(round(max(R[, 
		1])/365.24, 1)), sum(R[, 2] == 1), sum(R[, 2] == 2))
		x2 <- as.character(as.date(c(min(R[, 3]), max(R[, 3]))))
		paste("  age ranges from", x[1], "to", x[2], "years\n", " male:", 
		x[3], " female:", x[4], "\n", " date of entry from", 
		x2[1], "to", x2[2], "\n")
	}
	out
}




joinrate <- function(tables,dim.name="country"){
	nfiles <- length(tables)
	if(is.null(names(tables))) names(tables) <- paste("D",1:nfiles,sep="")
	if(any(!unlist(lapply(tables,is.ratetable))))stop("Tables must be in ratetable format")
	if(length(attributes(tables[[1]])$dim)!=3)stop("Currently implemented only for ratetables with 3 dimensions")

	
	for(it in 2:nfiles){
		if(length(attributes(tables[[it]])$dimid)!=3)stop("Each ratetable must have 3 dimensions: age, year and sex")
		mc <- match(attributes(tables[[it]])$dimid,attributes(tables[[1]])$dimid,nomatch=0)
		if(any(mc)==0) stop("Each ratetable must have 3 dimensions: age, year and sex")
		if(any(mc!=1:3)){
			atts <- attributes(tables[[it]])
			tables[[it]] <- aperm(tables[[it]],mc)
			atts$dimid <- atts$dimid[mc]
			atts$dimnames <- atts$dimnames[mc]
			atts$cutpoints <- atts$cutpoints[mc]
			atts$factor <- atts$factor[mc]
			atts$type <- atts$type[mc]
			atts$dim <- atts$dim[mc]
			attributes(tables[[it]]) <- atts
		}
	}
	
	list.eq <- function(l1,l2){
		n <- length(l1)
		rez <- rep(TRUE,n)
		for(it in 1:n){
			if(length(l1[[it]])!=length(l2[[it]]))rez[it] <- FALSE
			else if(any(l1[[it]]!=l2[[it]]))rez[it] <- FALSE
		}
		rez
	}
	
		
	equal <- rep(TRUE,3)
	for(it in 2:nfiles){
		equal <- equal*list.eq(attributes(tables[[1]])$cutpoints,attributes(tables[[it]])$cutpoints)
	}
		
	
	kir <-  which(!equal)
		
	newat <- attributes(tables[[1]])
	imena <- list(d1=NULL,d2=NULL,d3=NULL)
	
	for(jt in kir){
		listy <- NULL
		for(it in 1:nfiles){
			listy <- c(listy,attributes(tables[[it]])$cutpoints[[jt]])
		}
		imena[[jt]] <- names(table(listy)[table(listy) == nfiles])
		if(!length(imena[[jt]]))stop(paste("There are no common cutpoints for dimension", attributes(tables[[1]])$dimid[jt]))
	}
	
	
	for(it in 1:nfiles){
		keep <- lapply(dim(tables[[it]]),function(x)1:x)
		for(jt in kir){
			meci <- which(match(attributes(tables[[it]])$cutpoints[[jt]],imena[[jt]],nomatch=0)!=0)
			
			if(it==1){
				newat$dimnames[[jt]] <- attributes(tables[[it]])$dimnames[[jt]][meci] 
				newat$dim[[jt]] <- length(imena[[jt]])
				newat$cutpoints[[jt]] <- attributes(tables[[it]])$cutpoints[[jt]][meci]
			}
			if(length(meci)>1){if(max(diff(meci)!=1))warning(paste("The cutpoints for ",attributes(tables[[1]])$dimid[jt] ," are not equally spaced",sep=""))}
			keep[[jt]] <- meci		
		}
		tables[[it]] <- tables[[it]][keep[[1]],keep[[2]],keep[[3]]]
	}
	dims <- newat$dim
	out <- array(tables[[1]],dim=c(dims,1))
	for(it in 2:nfiles){
		out1 <- array(tables[[it]],dim=c(dims,1))
		out <- array(c(out,out1),dim=c(dims,it))
	}
	mc <- 1:4
	if(any(newat$factor>1)){
		wh <- which(newat$factor>1)
		mc <- c(mc[-wh],wh)
		out <- aperm(out,mc)
	}
	newat$dim <- c(dims,nfiles)[mc]
	newat$dimid <- c(newat$dimid,dim.name)[mc]
	newat$cutpoints <- list(newat$cutpoints[[1]],newat$cutpoints[[2]],newat$cutpoints[[3]],NULL)[mc]
	newat$factor <- c(newat$factor,1)[mc]
	newat$type <- c(newat$type,1)[mc]
	newat$dimnames <- list(newat$dimnames[[1]],newat$dimnames[[2]],newat$dimnames[[3]],names(tables))[mc]
	attributes(out) <- newat
	out
}
 

srvxp.fit <- function (x, y, ratetable) 
{
    x <- cbind(1:nrow(x), as.matrix(x))
    if (ncol(x) != (1 + length(dim(ratetable)))) 
        stop("x matrix does not match the rate table")
    atts <- attributes(ratetable)
    rfac <- atts$factor
    if (length(rfac) != ncol(x) - 1) 
        stop("Wrong length for rfac")
    ngrp <- nrow(x)
    times <- max(y)
    death <- TRUE
    if (times < 0) 
        stop("Negative time point requested")
    ntime <- 1
    cuts <- atts$cutpoints
    us.special <- (rfac > 1)
    if (any(us.special)) {
        if (sum(us.special) > 1) 
            stop("Two columns marked for special handling as a US rate table")
        cols <- 1 + match(c("age", "year"), attr(ratetable, "dimid"))
        if (any(is.na(cols))) 
            stop("Ratetable does not have expected shape")
        temp <- date.mdy(as.date(x[, cols[2]]) - x[, cols[1]])
        x[, cols[2]] <- x[, cols[2]] - mdy.date(temp$month, temp$day, 
            1960)
        temp <- (1:length(rfac))[us.special]
        nyear <- length(cuts[[temp]])
        nint <- rfac[temp]
        cuts[[temp]] <- round(approx(nint * (1:nyear), cuts[[temp]], 
            nint:(nint * nyear))$y - 1e-04)
    }
    temp <- .C("pyears3", as.integer(death), as.integer(nrow(x)), 
        as.integer(length(atts$dim)), as.integer(rfac), as.integer(atts$dim), 
        as.double(unlist(cuts)), ratetable, as.double(x), as.double(y), 
        as.integer(ntime), as.integer(ngrp), as.double(times), 
        surv = double(ntime * ngrp), n = integer(ntime * ngrp), 
        PACKAGE = "survival")
    temp$surv
}

 mlfit <- function (b, p, x, offset, d, h, ds, y, maxiter, tol) 
{
    for (nit in 1:maxiter) {
        b0 <- b
        fd <- matrix(0, p, 1)
        sd <- matrix(0, p, p)
        if (nit == 1) {
            ebx <- exp(x %*% b) * exp(offset)
            l0 <- sum(d * log(h + ebx) - ds - y * ebx)
        }
        for (it in 1:p) {
            fd[it, 1] <- sum((d/(h + ebx) - y) * x[, it] * ebx)
            for (jt in 1:p) sd[it, jt] = sum((d/(h + ebx) - d * 
                ebx/(h + ebx)^2 - y) * x[, it] * x[, jt] * ebx)
        }
        b <- b - solve(sd) %*% fd
        ebx <- exp(x %*% b) * exp(offset)
        l <- sum(d * log(h + ebx) - ds - y * ebx)
        bd <- abs(b - b0)
        if (max(bd) < tol) 
            break()
    }
    out <- list(b = b, sd = sd, nit = nit, loglik = c(l0, l))
    out
}

print.rs.br <- function (x, digits = max(options()$digits - 4, 3), ...) 
{
    invisible(print(x$table, digits = digits))
    if (x$rho != 0) 
        invisible(cat("Weighted Brownian bridge with rho=", x$rho, 
            "\n"))
}

print.rsadd <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "", "\n")
    if (length(coef(x))) {
        cat("Coefficients")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    if(x$method=="EM")
    		cat("\n", "Expected number of disease specific deaths: ",format(round(sum(x$Nie),2))," = ",format(round(100*sum(x$Nie)/sum(x$data$stat),1)),"% \n" ,sep="")
    if(x$method=="EM"|x$method=="max.lik"){
        	chi <- 2*max((x$loglik[2]-x$loglik[1]),0)
        	if(x$method=="EM")df <- length(x$coef)
        	else df <- length(x$coef)-length(x$int)+1
        	if(df>0){
        		p.val <- 1- pchisq(chi,df)
        		if(x$method=="max.lik")cat("\n")
        		cat("Likelihood ratio test=",format(round(chi,2)),", on ",df," df, p=",format(p.val),"\n",sep="")
        	}
        	else cat("\n")
    }
    cat("n=",nrow(x$data),sep="")
    if(length(x$na.action))cat("  (",length(x$na.action)," observations deleted due to missing)",sep="")
    cat("\n")
     if (length(x$warnme)) 
            cat("\n", x$warnme, "\n\n")
   else cat("\n")
    invisible(x)
}

summary.rsadd <- function (object, correlation = FALSE, symbolic.cor = FALSE, 
    ...) 
{
    if (inherits(object, "glm")) {
        p <- object$rank
        if (p > 0) {
            p1 <- 1:p
            Qr <- object$qr
            aliased <- is.na(coef(object))
            coef.p <- object$coefficients[Qr$pivot[p1]]
            covmat <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
            dimnames(covmat) <- list(names(coef.p), names(coef.p))
            var.cf <- diag(covmat)
            s.err <- sqrt(var.cf)
            tvalue <- coef.p/s.err
            dn <- c("Estimate", "Std. Error")
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "z value", "Pr(>|z|)"))
            df.f <- NCOL(Qr$qr)
        }
        else {
            coef.table <- matrix(, 0, 4)
            dimnames(coef.table) <- list(NULL, c("Estimate", 
                "Std. Error", "t value", "Pr(>|t|)"))
            covmat.unscaled <- covmat <- matrix(, 0, 0)
            aliased <- is.na(coef(object))
            df.f <- length(aliased)
        }
        ans <- c(object[c("call", "terms", "family", "iter", 
            "warnme")], list(coefficients = coef.table, var = covmat, 
            aliased = aliased))
        if (correlation && p > 0) {
            dd <- s.err
            ans$correlation <- covmat/outer(dd, dd)
            ans$symbolic.cor <- symbolic.cor
        }
        class(ans) <- "summary.rsadd"
    }
    else if (inherits(object, "rsadd")) {
        aliased <- is.na(coef(object))
        coef.p <- object$coef
        var.cf <- diag(object$var)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        dn <- c("Estimate", "Std. Error")
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
            "Pr(>|z|)"))
        ans <- c(object[c("call", "terms", "iter", "var")], list(coefficients = coef.table, 
            aliased = aliased))
        if (correlation && sum(aliased) != length(aliased)) {
            dd <- s.err
            ans$correlation <- object$var/outer(dd, dd)
            ans$symbolic.cor <- symbolic.cor
        }
        class(ans) <- "summary.rsadd"
    }
    else ans <- object
    return(ans)
}

print.summary.rsadd <- function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(x$aliased) == 0) {
        cat("\nNo Coefficients\n")
    }
    else {
        cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    if (length(x$warnme)) 
        cat("\n", x$warnme, "\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.col = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}

epa <- function(fit,bwin,times,n.bwin=16,left=FALSE){
	#bwin ... width of the window, relative to the default (1)
	#fit ... EM fit
	#times... times at which the smoothed plot is calculated
	#n.bwin ... number of different windows
	#left ... only predictable smoothing
	utd <- fit$times
	if(missing(times))times <- seq(1,max(utd),length=100)
	if(max(times)>max(utd)){
	warning("Cannot extrapolate beyond max event time")
	times <- pmax(times,max(utd))
	}
	nutd <- length(utd)
	nt4 <- c(1,ceiling(nutd*(1:n.bwin)/n.bwin))
	if(missing(bwin))bwin <- rep(length(fit$times)/100,n.bwin)
	else bwin <- rep(bwin*length(fit$times)/100,n.bwin)
	for(it in 1:n.bwin){
		bwin[it] <- bwin[it]*max(diff(utd[nt4[it]:nt4[it+1]]))
	}
	while(utd[nt4[2]]<bwin[1]){		# ce je bwin velik, skrajsamo nt4
	       nt4 <- nt4[-2]
	       if(length(nt4)==1)break
	}
	#the smoothing matrix
	if(left) krn <- kernerleftch(utd,bwin,nt4)
	else krn <- kern(times,utd,bwin,nt4)
	lams <- pmax(krn%*%fit$lam0.ns,0)
	list(lambda=lams,times=times)				#	, weights=c(fit$times[1],diff(fit$times)))
}

Kern <- function (t, tv, b, tD, nt4) 
{
    Rb <- max(tv)					#Right border	
    kmat <- NULL
    tvs <- tv
    tv <- tv[-1]
    kt <- function(q,t)12*(t+1)/(1+q)^4*( (1-2*q)*t + (3*q^2-2*q+1)/2  )
    totcajti <- NULL
    for (it in 1:(length(nt4) - 1)) {
     	cajti <- t[t>tvs[nt4[it]] & t<=tvs[nt4[it + 1]]]
    	if(length(cajti)){
	    q <- min( cajti/b[it],1,(Rb-cajti)/b[it])
   	    if(q<1 & length(cajti)>1){
   	    	jc <- 1
   	    	while(jc <=length(cajti)){
   	    		qd <- pmin( cajti[jc:length(cajti)]/b[it],1,(Rb-cajti[jc:length(cajti)])/b[it])
			q <- qd[1]
			if(q==1){
			casi <- cajti[jc:length(cajti)][qd==1]
			q <- 1
			jc <- sum(qd==1)+jc
			}
			else{
			casi <- cajti[jc]
			jc <- jc+1
			}
			kmat1 <- outer(casi, tv, "-")/b[it]		#z - to je ok
	                if(q<1){
	                if(casi>b[it]) kmt1 <- -kmat1
	                vr <- kt(q,kmat1)*(kmat1>=-1 & kmat1 <= q)
	                }
   	     		else vr <-  pmax((1 - kmat1^2) * .75,0)
   	 		kmat <- rbind(kmat, vr/b[it])
   	 		totcajti <- c(totcajti,casi)
   	 	}
	   }
 	   else{
    	       kmat1 <- outer(cajti, tv, "-")/b[it]		#z - to je ok
	       q <- min( cajti/b[it],1)
   	       if(q<1)vr <- kt(q,kmat1)*(kmat1>=-1 & kmat1 <= q)
   	       else vr <-  pmax((1 - kmat1^2) * .75,0)
   	       kmat <- rbind(kmat, vr/b[it])
   	       totcajti <- c(totcajti,cajti)
   	   }#else
   	}#if
   	 
    }#for
    kmat
}
 
kern <- function (times,td, b, nt4) 
{
    n <- length(td)
    ttemp <- td[td >= b[1]]
    ntemp <- length(ttemp)
    if (ntemp == n) 
        nt4 <- c(0, nt4[-1])
    td <- c(0,td)
    nt4 <- c(1,nt4+1)
    b <- c(b[1],b)
    krn <- Kern(times, td, b, max(td), nt4)
    krn
}

rs.surv <- function (formula=formula(data), data = parent.frame(), ratetable = survexp.us,na.action,fin.date,method="weighted.ederer",conf.type="log",conf.int=0.95) 
{	
    call <- match.call()
    rform <- rformulate(formula, data, ratetable, na.action)
    
    data <- rform$data
    
    
    #for the Hakulinen method - the potential censoring times - only infinite values implemented for now
    if(method=="hakulinen"){
    R <- rform$R
    coll <- match("year",attributes(ratetable)$dimid)
    year <- R[,coll]	
    if(missing(fin.date)) fin.date <- max(data$Y+year)			#set the final date to the last observed date
    data$Y2 <- data$Y							
    if(length(fin.date==1))data$Y2[data$stat==1] <- fin.date - year[data$stat==1]	 #Y2=potential follow-up time   	  
    else if(length(fin.date==nrow(data))) data$Y2[data$stat==1] <- fin.date[data$stat==1] - year[data$stat==1]	    	  	    	  
    else stop("fin.date must be either one value of a vector of the same length as the data")
    data$stat2 <- rep(0,nrow(data))					#everyone is censored at potential follow-up time
    }
    
    inx.d <- order(data$Y,(1-data$stat))				#the indicator for ordering in time
    data <- data[inx.d,]
    p <- rform$m		#number of covariates defining the strata
    ti <- sort(unique(data$Y))	#unique follow-up times
    if(method=="hakulinen") ti <- sort(unique(pmin(c(ti,data$Y2),max(ti))))
    Ki <- matrix(NA,nrow=nrow(data),ncol=length(ti))
    dNi <- matrix(0,nrow=nrow(data),ncol=length(ti))
    dYi <- matrix(0,nrow=nrow(data),ncol=length(ti))
    dYi.hak <- matrix(0,nrow=nrow(data),ncol=length(ti))
    #cajt <- NULL
    nfk <- length(attributes(rform$ratetable)$dimid)
    #cajt <- proc.time()[3]
    for(jt in 1:length(ti)){
    Ki[,jt] <- srvxp.fit(data[, 4:(nfk + 3)], rep(ti[jt],nrow(data)), rform$ratetable)
    dNi[which(data$Y==ti[jt]),jt] <- data$stat[which(data$Y==ti[jt])]
    dYi[which(data$Y==ti[jt]),jt] <- 1 - data$stat[which(data$Y==ti[jt])]
    if(method=="hakulinen") dYi.hak[which(data$Y2==ti[jt]),jt] <- 1
    }
    #cajt <- c(cajt,proc.time()[3])
    #datalong <- cbind(rep(data[,4],length(ti)),rep(data[,5],length(ti)),rep(data[,6],length(ti)))
    #Ki2 <- srvxp.fit(datalong, rep(ti,each=nrow(data)), rform$ratetable)
    #Ki3 <- matrix(Ki2,nrow=nrow(data))
    #cajt <- c(cajt,proc.time()[3])
    
    #dKi <- -log(Ki) + log(cbind(rep(1,nrow(data)),Ki[,-ncol(Ki)]))		#d\Lambda_{Pi}
    # dKi2 <- cbind(-log(Ki[,1]),t(apply(-log(Ki),1,diff)))
    
    if(p>0)data$Xs <- strata(rform$X[inx.d,])						#the stratification group
    else data$Xs <- rep(1,nrow(data))
  
    se.fac <- sqrt(qchisq(conf.int,1))						#finds the factor needed for the confidence intervals, e.g. 1.96
    out <- NULL
    out$n <- table(data$Xs)

    out$time <- out$n.risk <- out$n.event <- out$n.censor  <- out$surv<- out$std.err <- out$strata <-  NULL
    method <-  match.arg(method,c("weighted.ederer","ederer2","hakulinen"))
   
   
    for(kt in 1:length(out$n)){
       
    	inx <- which(data$Xs==names(out$n)[kt])				#the subjects in this stratum
    	dYis <- dYi[inx,]
    	dNis <- dNi[inx,]
    	Kis <- Ki[inx,]
    	#dKis <- dKi[inx,]
    	dYis.hak <- dYi.hak[inx,]
    	datas <- data[inx,]
    	yinx <- which(apply(dYis+dNis+dYis.hak,2,sum)!=0)			#ti that are relevant for this stratum
    	dYis <- dYis[,yinx]
    	dNis <- dNis[,yinx]
    	Kis <-  Kis[,yinx]
    	tis <- ti[yinx]
    	#dKis <-  dKis[,yinx]
    	dYis.hak <- dYis.hak[,yinx]
        Nis <- t(apply(dNis,1,cumsum))						#jumps to 1 in case of event
        Yis <- t(apply(dYis,1,cumsum))						#jumps to 1 in case of censoring
        Yis.hak <- t(apply(dYis.hak,1,cumsum))					#jumps to 1 in case of potential censoring
        #Yis.hak <-  matrix(1,nrow=nrow(Yis),ncol=ncol(Yis))- Yis 		#jumps to 0 in case of cesnoring
        Yis <-  matrix(1,nrow=nrow(Yis),ncol=ncol(Yis))- Nis - Yis 		#jumps to 0 in case of event or censoring
        Yis <- cbind(rep(1,nrow(Yis)),Yis[,-ncol(Yis)])				#starts with 1, right continuous
    	Yis.hak <- matrix(1,nrow=nrow(Yis),ncol=ncol(Yis))- Yis.hak		#jumps to 0 in case of potential follow-up time
    	Yis.hak <-  cbind(rep(1,nrow(Yis.hak)),Yis.hak[,-ncol(Yis.hak)])	#starts with 1	
    	#delete the columns with 0 at risk (needed for Hakulinen method)
    	yinx <- which(apply(Yis,2,sum)>0)				
    	dYis <- dYis[,yinx]
    	dNis <- dNis[,yinx]
    	Kis <-  Kis[,yinx]
    	tis <- tis[yinx]
    	Yis <- Yis[,yinx]
	Nis <- Nis[,yinx]    	
    	Yis.hak <- Yis.hak[,yinx]
    	   	
    	dKis <- -log(Kis) + log(cbind(rep(1,nrow(Kis)),Kis[,-ncol(Kis)]))	#d\Lambda_{Pi}: -log(S_i(t_i))+ log(S_i(t_{i-1}))
    	out$time <- c(out$time,tis)					
    	out$n.risk <- c(out$n.risk,apply(Yis,2,sum))
    	out$n.event <- c(out$n.event,apply(dNis,2,sum))
    	out$n.censor <- c(out$n.censor,apply(dYis,2,sum))
    	    	
     	if(method=="ederer2"){
    		out$surv <- c(out$surv,exp(-cumsum(apply((dNis-Yis*dKis),2,sum)/apply(Yis,2,sum))))
    		out$std.err <- c(out$std.err,sqrt(cumsum(apply(dNis,2,sum)/(apply(Yis,2,sum))^2)))
    	}
    	else if(method=="weighted.ederer"){
    		out$surv <- c(out$surv,exp(-cumsum(apply(dNis/Kis,2,sum)/apply(Yis/Kis,2,sum) - apply(Yis/Kis*dKis,2,sum)/apply(Yis/Kis,2,sum))))
	    	out$std.err <- c(out$std.err,sqrt(cumsum(apply(dNis/Kis,2,sum)/(apply(Yis/Kis,2,sum))^2)) )
    	}
       	else if(method=="hakulinen"){ #
	       	#d <- -log(apply((Yis.hak*Ki),2,sum)/apply(Yis.hak,2,sum))
	    	#dd <- c(d[1],diff(d))
	    	#e <- apply(Yis.hak * Kis*exp(-dKis), 2, sum)/apply(Yis.hak*Kis, 2, sum)
		#e <- cumprod(e)
		f <- apply(Yis.hak * Kis*dKis, 2, sum)/apply(Yis.hak*Kis, 2, sum)
		#out$surv <- c(out$surv,exp(-cumsum(apply((dNis),2,sum)/apply(Yis,2,sum)))/e)
		out$surv <- c(out$surv,exp(-cumsum( apply(dNis,2,sum)/apply(Yis,2,sum) - f  )))
		#out$surv <- c(out$surv,exp(-cumsum(apply((dNis),2,sum)/apply(Yis,2,sum)-dd)))       		
       		#out$surv <- c(out$surv,exp(-cumsum(apply((dNis),2,sum)/apply(Yis,2,sum)-apply((Yis.hak*dKis),2,sum)/apply(Yis.hak,2,sum))))
	    	out$std.err <- c(out$std.err,sqrt(cumsum(apply(dNis,2,sum)/(apply(Yis,2,sum))^2)))
	    	
	}
    	out$strata <- c(out$strata,length(tis))
    	
    }
    
    if(conf.type=="plain"){
	    out$lower <- as.vector(out$surv - out$std.err*se.fac*out$surv)
	    out$upper <- as.vector(out$surv + out$std.err*se.fac*out$surv)
    }
    else if(conf.type=="log"){
   	    out$lower <- exp(as.vector(log(out$surv) - out$std.err*se.fac))
  	    out$upper <- exp(as.vector(log(out$surv) + out$std.err*se.fac))
    }
    else if(conf.type=="log-log"){
   	    out$lower <- exp(-exp(as.vector(log(-log(out$surv)) - out$std.err*se.fac/log(out$surv))))
  	    out$upper <- exp(-exp(as.vector(log(-log(out$surv)) + out$std.err*se.fac/log(out$surv))))
    }    
    names(out$strata) <- names(out$n)
    if(p==0)out$strata <- NULL					#hide it if no strata (important for the behaviour of plot.survfit)
    out$n <- as.vector(out$n)
    out$conf.type <- conf.type
    out$conf.int <- conf.int
    out$method <-  method
    out$call <- call
    out$type <- "right"
    class(out) <- c("survfit","rs.surv") 
    out
}




rs.survo <- function (formula, data,ratetable=survexp.us,fin.date,method="hakulinen",...) 
{	
    call <- match.call()
     if (!inherits(formula, "formula")) 
            temp <- UseMethod("rs.surv")

    if (mode(call[[2]]) == "call" & call[[2]][[1]] == as.name("Surv")) {
                formula <- eval(parse(text = paste(deparse(call[[2]],width.cutoff = 500L), 
                    1, sep = "~")))
                environment(formula) <- parent.frame()
    }

    else if (mode(call[[2]]) == "call" & any(regexpr("ratetable",call[[2]][[3]])>0)){
                formula <- eval(parse(text = paste(deparse(call[[2]][[2]],width.cutoff = 500L),"~1+",deparse(call[[2]][[3]],width.cutoff = 500L) 
                    , sep = "")))
                environment(formula) <- parent.frame()
    }

    Terms <- if (missing(data)) terms(formula, "ratetable")
    else terms(formula, "ratetable", data = data)
      
    
    rate <- attr(Terms, "specials")$ratetable

    if (length(rate) > 1) 
        stop("Can have only 1 ratetable() call in a formula")
    if (length(rate) == 0) {
        xx <- function(x) formula(x)
        if (is.ratetable(ratetable)) 
            varlist <- attr(ratetable, "dimid")
        else stop("Invalid rate table")
        ftemp <- deparse(formula)
        ftemp <- paste(ftemp,collapse="")
        formula <- xx(paste(ftemp, "+ ratetable(", paste(varlist, 
            "=", varlist, collapse = ","), ")"))
        Terms <- if (missing(data)) 
            terms(formula, "ratetable")
        else terms(formula, "ratetable", data = data)
        rate <- attr(Terms, "specials")$ratetable
    }
     m <- match.call(expand = FALSE)
     m$ratetable <- m$fin.date <- m$method  <- m$... <- NULL
  

    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    n <- nrow(m)
    Y <- model.extract(m, "response")
    status <- Y[,2]
     
    rtemp <- match.ratetable(m[, rate], ratetable)
    R <- rtemp$R
    coll <- match("year",attributes(ratetable)$dimid)
    year <- R[,coll]	
	
    if(missing(fin.date)) fin.date <- max(year+Y[,1])
   
    
    formula.obs <- formula
    
    if (!is.null(rate)) formula.obs[[3]] <- formula.obs[[3]][[2]]
   
    obs.hak <- survfit(formula.obs,data,...)	
    times <- obs.hak$time
    times <- unique(sort(times))
    nc <- 1
    o1 <- summary(obs.hak, time=times)
    if(!is.null(obs.hak$strata))   nc<-length(obs.hak$strata.all)
    
    k <- rep(1,length(o1$time))
    ind <- o1$time < c(0,o1$time[-length(o1$time)])
    kend <- length(k)
    
    for (i in which(ind)){
    	k[i:kend] <- k[i:kend]+1
    }
    method <-  match.arg(method,c("hakulinen","conditional"))

    if(method=="hakulinen"){
 	 data$time1 <- Y[,1]
         wh <- which(status==1)
         data$time1[wh] <- fin.date - year[wh]
         data$status <- status	
         data <- data[data$time1>0,]
         data <- data[!is.na(data$time1),]
    	 xx <- function(x) formula(x)
  	condi <- FALSE
	 survex <- "Surv(time1,status)"
	  formula.exp <- xx(paste(survex,paste(names(m)[-1],collapse="+"),sep="~"))

    }
    else{
   	condi <- TRUE
   	formula.exp <- formula
    }
    deli <- function(x=o1$surv,kji=k,y=exp.hak$surv){
    	xnew <- x
    	for(it in unique(sort(k))){
    		d <- sum(kji==it)
    		nom <- x[kji==it]
    		if(length(unique(k))>1)	denom <- y[1:d,it]
    		else denom <- y[1:d]
    		nom[denom==0] <- 0
    		denom[denom==0] <- 1
    		xnew[kji==it] <- nom/denom
       	}
       	xnew
    }
    
    exp.hak <- survexp(formula.exp,data,conditional=condi,cohort=TRUE,ratetable=ratetable,times=times)
   
    rel.hak <- deli()
    
    se.fac <- sqrt(qchisq(obs.hak$conf.int,1))
   
    #serr <- deli(x=o1$serr,y=exp.hak$surv)

    
    ktime <- NULL
    for(it in unique(sort(k))){
    	d <- sum(k==it)
    	ktime <- c(ktime,times[1:d])
    }
    
    out <- obs.hak
    out$surv <- as.vector(rel.hak)
    out$time <- ktime
    out$n.risk <- o1$n.risk
    out$n.event <- o1$n.event
    out$std.err <- o1$std.err


    if(obs.hak$conf.type=="plain"){
    out$lower <- as.vector(rel.hak - out$std.err*se.fac*rel.hak)
    out$upper <- as.vector(rel.hak + out$std.err*se.fac*rel.hak)
    }
    if(obs.hak$conf.type=="log"){
    out$lower <- exp(as.vector(log(rel.hak) - out$std.err*se.fac))
    out$upper <- exp(as.vector(log(rel.hak) + out$std.err*se.fac))
    }
    if(obs.hak$conf.type=="log-log"){
        out$lower <- exp(-exp(as.vector(log(-log(rel.hak)) - out$std.err*se.fac/log(rel.hak))))
        out$upper <- exp(-exp(as.vector(log(-log(rel.hak)) + out$std.err*se.fac/log(rel.hak))))
    }
    
    out$ntimes.strata <- table(k)
    out$strata <- out$ntimes.strata
    names(out$strata) <- names(obs.hak$strata)
    out$popul <- exp.hak$surv
       
    out
}

