rformulate <- function(formula,data,ratetable,na.action,subset){
    	
    	call <- match.call()
    	m <- match.call(expand = FALSE)
    	m$ratetable <- NULL
    	Terms <- if (missing(data)) 
    	    terms(formula, "ratetable")
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
	    }
	    else if (attr(Y, "type") == "counting") {
	        type <- attr(Y, "type")
	        status <- Y[, 3]
	        start <- Y[, 1]
	        Y <- Y[, 2]
    	    }
    	    else stop("Illegal response value")
    	    if (any(c(Y, start) < 0)) 
        stop("Negative follow up time")
    if (is.ratetable(ratetable)) {
        israte <- TRUE
        rtemp <- match.ratetable(m[, rate], ratetable)
        R <- rtemp$R
        if (!is.null(rtemp$call)) {
            ratetable <- eval(parse(text = rtemp$call))
        }
    }
    else stop("Invalid ratetable argument")
    if(rate==2) {
    	X<-NULL
    	mm<-0
    }
    else {
    	if(length(rate==1)){
                    	formula[[3]] <- formula[[3]][[2]]
    	}
    	X <-as.data.frame(model.matrix(formula,data=data))[,-1,drop=FALSE]
    	mm <- ncol(X)
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

    out <- list(R=R,status=status,start=start,Y=Y,X=X,m=mm,n=n,type=type,Y.surv=Y.surv,Terms=Terms,ratetable=ratetable,offset=offset)
    na.action <- attr(m, "na.action")
    if (length(na.action))  out$na.action <- na.action
    out
        
    }
    
    