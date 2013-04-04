cmp.rel <-  function (formula = formula(data), data = parent.frame(), ratetable = survexp.us, 
     na.action,tau,conf.int=0.95) 
    
    #formula: for example Surv(time,cens)~1 #not implemented for subgroups - DO IT!
    #data: the observed data set
    #ratetable: the population mortality tables
    #conf.type: confidence interval calculation (plain, log or log-log)
    #conf.int: confidence interval 
    #tau: max. cas do katerega racuna
{

    call <- match.call()
    rform <- rformulate(formula, data, ratetable, na.action)			#get the data ready
    data <- rform$data								#the data set
    se.fac <- sqrt(qchisq(conf.int, 1))						#factor needed for confidence interval

    if(missing(tau)) tau<-max(rform$Y)

    p <- rform$m								#number of covariates
    if (p > 0) #if covariates
    data$Xs <- strata(rform$X[, ,drop=FALSE ])					#make strata according to covariates
    else data$Xs <- rep(1, nrow(data))						#if no covariates, just put 1
    
    tab.strata <- table(data$Xs)						#unique strata values
    ntab.strata <- length(tab.strata)						#number of strata
    
    dtemp <- list(NULL)
    out <- as.list(rep(dtemp,ntab.strata*2))

    for (kt in 1:ntab.strata) {							#for each stratum
        inx <- which(data$Xs == names(tab.strata)[kt])				#individuals within this stratum

	tis <- sort(unique(pmin(tau,rform$Y[inx])))	
	k <- length(tis)
	
	out[[2*kt-1]]$time <- tis
	out[[2*kt]]$time <- tis
	
	temp <- exp.prep(rform$R[inx,,drop=FALSE],rform$Y[inx],ratetable,rform$status[inx],times=tis,fast=TRUE)		#calculate the values for each interval of time

	dLambdap <- temp$yidli/temp$yi
	dLambdao <- temp$dni/temp$yi
	dLambdae <- dLambdao - dLambdap

	sigma <- temp$dni/(temp$yi)^2
	sigmap <- temp$yidli/(temp$yi)^2		
	sigmae <- sigma - sigmap						  #variance
    
  
    	So <- cumprod(1-dLambdao)						  #S(t)
   	Soprej <- c(1,So[-length(So)])						  #S(t-)
    
  	cumince <- cumsum(Soprej*dLambdae)
 	cumincp <- cumsum(Soprej*dLambdap)
 	ve <- vp <- rep(NA,k)
	for(it in 1:k){
		ve[it] <- sum((cumince[it] - cumince[1:it])^2*sigma[1:it]) + sum(So[1:it]*sigmae[1:it]*(So[1:it]-2*(cumince[it]-cumince[1:it])))
		vp[it] <- sum((cumincp[it] - cumincp[1:it])^2*sigma[1:it]) + sum(So[1:it]*sigmap[1:it]*(So[1:it]-2*(cumincp[it]-cumincp[1:it])))
	}		
  
	areae <- sum(diff(c(0,tis))*cumince)/365.241
	areap <- sum(diff(c(0,tis))*cumincp)/365.241
  
  	out[[2*kt-1]]$est <- cumince
  	out[[2*kt-1]]$var <- ve
  	out[[2*kt-1]]$lower <- cumince-1.96*sqrt(ve)
  	out[[2*kt-1]]$upper <- cumince+1.96*sqrt(ve)
  	out[[2*kt-1]]$area <- areae
  	
  	out[[2*kt]]$est <- cumincp
	out[[2*kt]]$var <- vp
	out[[2*kt]]$lower <- cumincp-1.96*sqrt(vp)
	out[[2*kt]]$upper <- cumincp+1.96*sqrt(vp)
	out[[2*kt]]$area <- areap
   }
   if(p>0)names(out) <- paste(rep(c("causeSpec","population"),ntab.strata),rep(names(tab.strata),each=2))
   else names(out) <- c("causeSpec","population")
   class(out) <- "cmp.rel"
   out
}


plot.cmp.rel <- function (x, main = " ", curvlab, ylim = c(0, 1), xlim, wh = 2, 
    xlab = "Time", ylab = "Probability", lty = 1:length(x), xscale=1,
    col = 1, lwd = par("lwd"), ...) 
{
#wh= upper left coordinates of the legend, if of length 1, the legend is placed top left.
    nc <- length(x)			#number of curves
    if (length(lty) < nc)		#if not enough different line types 
        lty <- rep(lty[1], nc)
    else lty <- lty[1:nc]
    if (length(lwd) < nc)		#if not enough different line widths 
        lwd <- rep(lwd[1], nc)
    else lwd <- lwd[1:nc]
    if (length(col) < nc)		#if not enough different colors 
        col <- rep(col[1], nc)
    else col <- col[1:nc]
    if (missing(curvlab)) {		#if no curve labels desired
        if (mode(names(x)) == "NULL") { #and no curve labels prespecified
            curvlab <- as.character(1:nc)
        }
        else curvlab <- names(x)[1:nc]	#use prespecified if they exist
    }
    if (missing(xlim)) {		#if no limits desired
        xmax <- 0
        for (i in 1:nc) {
            xmax <- max(c(xmax, x[[i]][[1]]/xscale))	#take max time over all strata
        }
        xlim <- c(0, xmax)
    }
    plot(x[[1]][[1]]/xscale, x[[1]][[2]], type = "n", ylim = ylim, xlim = xlim, 
        main = main, xlab = xlab, ylab = ylab, bty = "l", ...)		#plot estimates [[1]]=time, [[2]]=est
    if (length(wh) != 2) {
        wh <- c(xlim[1], ylim[2])
    }
    u <- list(...)
    if (length(u) > 0) {
        i <- pmatch(names(u), names(formals(legend)), 0)
        do.call("legend", c(list(x = wh[1], y = wh[2], legend = curvlab, 
            col = col, lty = lty, lwd = lwd, bty = "n", bg = -999999), 
            u[i > 0]))
    }
    else {
        do.call("legend", list(x = wh[1], y = wh[2], legend = curvlab, 
            col = col, lty = lty, lwd = lwd, bty = "n", bg = -999999))
    }
    for (i in 1:nc) {
        lines(x[[i]][[1]]/xscale, x[[i]][[2]], lty = lty[i], col = col[i], 
            lwd = lwd[i], ...)
    }
}


print.cmp.rel <- function (x, ntp = 4, maxtime,xscale=365.241, ...) 
{
    nc <- length(x)

    if (missing(maxtime)) {
        maxtime <- 0
        for (i in 1:nc) maxtime <- max(maxtime, x[[i]]$time)
    }
    tp <- pretty(c(0, maxtime/xscale), ntp + 1)
    cat("Estimates and Variances:\n")
    print(timepoints(x, tp[-c(1, length(tp))],xscale), ...)
    invisible()
}

timepoints <- function (w, times,xscale) 
{
    ng <- length(w)
    times <- sort(unique(times))*xscale
    nt <- length(times)
    storage.mode(times) <- "double"
    storage.mode(nt) <- "integer"
    ind <- matrix(0, ncol = nt, nrow = ng)
    oute <- matrix(NA, ncol = nt, nrow = ng)
    outv <- oute
    storage.mode(ind) <- "integer"
    slct <- rep(TRUE, ng)
    for (i in 1:ng) {
        if (is.null((w[[i]])$est)) {
            slct[i] <- FALSE
        }
        else {
            z <- rep(NA,nt)
            for(kt in 1:nt)z[kt] <- rev(which(w[[i]][[1]]<times[kt]))[1]
            ind[i, ] <- z
            oute[i, ind[i, ] > 0] <- w[[i]][[2]][z]
            if (length(w[[i]]) > 2) 
                outv[i, ind[i, ] > 0] <- w[[i]][[3]][z]
        }
    }
    dimnames(oute) <- list(names(w)[1:ng], as.character(times/xscale))
    dimnames(outv) <- dimnames(oute)
    list(est = oute[slct, , drop = FALSE], var = outv[slct, , 
        drop = FALSE])
}