"glmxp" <-
function (rform,interval,method,control) 
{
    #grouping
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
    # calculations withing groups
    vg <- function(X) {
    	n <- dim(X)[1]
        w <- sum((X$event == 0) & (X$fin == 1) & (X$y!=1))
        nd <- sum((X$event == 1) & (X$fin == 1))
        ps <- srvxp.fit(X[, 4:(nfk + 3)], rep(t.int,n), rform$ratetable)
        ld <- n - w/2
        lny <- log(sum(X$y))
        k <- t.int/365.24
        dstar <- sum(-log(ps)/k*X$y)
        ps <- mean(ps)
        if(rform$m==0)data.rest <- X[1, 7 + nfk + rform$m,drop=FALSE]
        else  data.rest <- X[1, c((3 + nfk + 1):(3 + nfk + rform$m), 
        	7 + nfk + rform$m)]
        cbind(nd = nd, ld = ld, ps = ps, lny = lny, dstar = dstar, 
        	k=k,data.rest)
    }
    
    #organisation of intervals
    nint <- length(interval)
    if (nint < 2) 
            stop("Illegal interval value")
    meje <- interval
    my.fun <- function(x){
	if(x>1){
	    x.t<-rep(1,floor(x))
	    if(x-floor(x)>0)x.t<-c(x.t,x-floor(x))
	    x.t
	}
	else x
    }
    int<-apply(matrix(diff(interval),ncol=1),1,my.fun)
    if(is.list(int))int<-c(0,cumsum(do.call("c",int)))	
    else int <- c(0,cumsum(int))
    int <- int * 365.24
    nint <- length(int)
    
    #data set
    if(rform$m==0)X <- data.frame(start = rform$start, Y = rform$Y, stat = rform$status, rform$R, grupa = g)
    else X <- data.frame(start = rform$start, Y =rform$Y, stat = rform$status, rform$R, rform$X, grupa = g)
    
    # ratetable attributes
    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    
    # temp data frame for subjects not yet at risk
    Z <- X[X$start >= int[2], ]
    nz <- dim(Z)[1]
    Z$fin <- rep(0, nz)
    Z$event <- rep(0, nz)
    Z$fu <- rep(0, nz)
    Z$y <- rep(0, nz)
    Z$origstart <- Z$start
    Z$xind <- rep(0,nz)
    if(nrow(Z)>0)Z[, 4:(nfk + 3)] <- Z[, 4:(nfk + 3)] + matrix(Z$start, ncol = nfk, byrow = F, nrow = nrow(Z))*matrix(fk,ncol=nfk,byrow=T,nrow=nrow(Z))
    
    # subjects at risk on first interval
    X <- X[X$start < int[2], ]
    X$fin <- (X$Y < int[2])
    X$event <- X$fin * X$stat
    ford <- eval(substitute(paste("[", a, ",", b, ")", sep = ""), 
        list(a = meje[1], b = meje[2])))
    X$fu <- rep(ford, rform$n - nz)
    t.int <- int[2]-int[1]
    X$y <- (pmin(X$Y, int[2]) - X$start)/365.24
    X$origstart <- X$start
    X$xind <- rep(1,nrow(X))
    gr1 <- by(X, X$grupa, vg)
    grm1 <- do.call("rbind", gr1)
    X <- X[X$fin == 0, ]
    X$start <- rep(int[2], dim(X)[1])
    
    #subjects at risk on second interval
    X <- rbind(X, Z[Z$start < int[3], ])
    
    # temp data frame for subjects not yet at risk
    Z <- Z[Z$start >= int[3], ]
    temp <- 0
    
    # over intervals
    if (nint > 2) {
    	for (i in 3:nint) {
    	    ni <- dim(X)[1]
    	    if (ni == 0)  {temp<-1;break}
            X$fin <- X$Y < int[i]
            X$event <- X$fin * X$stat
            l<-sum(int[i-1]>=meje*365.24)
            ford <- c(ford, eval(substitute(paste("[", a, ",", 
                b, ")", sep = ""), list(a = meje[l], b = meje[l+1]))))
            X$fu <- rep(ford[i-1], ni)
            t.int <- int[i] - int[i - 1]
            index<- X$origstart < int[i-1]
            index1 <- as.logical(X$xind)
            if(sum(index)>0)X[index, 4:(nfk + 3)] <- X[index, 4:(nfk + 3)] + matrix(fk * t.int, ncol = nfk, byrow = TRUE, nrow = sum(index))
            X$xind <- rep(1,nrow(X))
            X$y <- (pmin(X$Y, int[i]) - X$start)/365.24
            gr1 <- by(X, X$grupa, vg)
            grm1 <- rbind(grm1, do.call("rbind", gr1))
            X <- X[X$fin == 0, ]
            X$start <- rep(int[i], dim(X)[1])
            if(i==nint) break
            X <- rbind(X, Z[Z$start < int[i + 1], ])
            X <- X[X$start!=X$Y,]
            Z <- Z[Z$start >= int[i + 1], ]
        }
        l<-sum(int[i-temp]>meje*365.24)
        interval <- meje[1:(l+1)]
    }
    else interval<-meje[1:2]
    grm1$fu <-factor(grm1$fu,levels=unique(ford))
    
    # binomial
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
            warnme <- paste("Observed number of deaths is smaller than the expected in ",n,"/",g, " groups of patients", sep="")
        }
        else warnme <- ""
        if (length(interval) == 2 & rform$m == 0) 
            stop("No groups can be formed")
        
        if (length(interval) == 1|length(table(grm1$fu))==1 ) grm1$fu <- as.integer(grm1$fu)
        if(is.null(rform$X)) local.ht <- glm(cbind(nd,ld-nd)~-1 + fu +offset(log(k)), data=grm1, family=ht )
        else{
        	xmat <- as.matrix(grm1[,7:(ncol(grm1)-1)])
		local.ht <- glm(cbind(nd,ld-nd)~-1+ xmat+ fu +offset(log(k)), data=grm1, family=ht )
	}
	names(local.ht[[1]])<-c(names(rform$X),paste("fu",levels(grm1$fu)))
    }
    #poisson
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
            warnme <- paste("Observed number of deaths is smaller than the expected in ",n,"/",g, " groups of patients",sep="")
           
        }
        else warnme <- ""
        dstar <- grm1$dstar
        if (length(interval)==2 & rform$m== 0) 
            stop("No groups can be formed")
        
        if (length(interval) == 1|length(table(grm1$fu))==1 ) grm1$fu <- as.integer(grm1$fu)
	if(is.null(rform$X)) local.ht <- glm(nd~-1 + fu, data=grm1, family=pot,offset=lny)
	else{
	      	xmat <- as.matrix(grm1[,7:(ncol(grm1)-1)])
		local.ht <- glm(nd~-1+ xmat+ fu, data=grm1, family=pot,offset=lny)
	}
	names(local.ht[[1]])<-c(names(rform$X),paste("fu",levels(grm1$fu)))
    }
    else stop(paste("Method '",method,"' not a valid method",sep=""))
    class(local.ht) <- c("rsadd",class(local.ht))
    local.ht$warnme <- warnme
    local.ht$int <- interval
    local.ht$groups <- local.ht$data
    return(local.ht)
}
