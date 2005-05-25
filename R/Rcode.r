rsadd <- function (formula = formula(data), data = parent.frame(), ratetable = survexp.us, 
    int, na.action, method = "max.lik", init, centered = FALSE, control, ...) 
{
    call <- match.call()
    if (missing(control)) 
        control <- glm.control(...)
   
   
    rform <- rformulate(formula, data, ratetable, na.action, 
        int, centered)
        
	if (missing(int)) 
	    int <- ceiling(max(rform$Y/365.24))
	if (length(int) == 1) {
	    if (int <= 0) 
		stop("The value of 'int' must be positive ")
	    int <- 0:int
	}
	else if (int[1] != 0) 
	    stop("The first interval in 'int' must start with 0")

    if (method == "glm.bin" | method == "glm.poi") 
        fit <- glmxp(rform = rform, interval = int, method = method, 
            control = control)
    else if (method == "max.lik") 
        fit <- maxlik(rform = rform, interval = int, init = init, 
            control = control)
    else stop("Unknown method")
    fit$call <- call
    fit$formula <- formula
    fit$data <- rform$data
    fit$ratetable <- rform$ratetable
    if (length(rform$na.action)) 
        fit$na.action <- rform$na.action
    fit$y <- rform$Y.surv
    fit$method <- method
    if (rform$m > 0) 
        fit$linear.predictors <- as.matrix(rform$X) %*% fit$coef[1:ncol(rform$X)]
    fit
}

rformulate <- function (formula, data = parent.frame(), ratetable, na.action, 
    int, centered) 
{
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$ratetable <- m$int <- m$centered <- NULL
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
    if (rate == 2) {
        X <- NULL
        mm <- 0
    }
    else {
        if (length(rate == 1)) {
            formula[[3]] <- formula[[3]][[2]]
        }
        X <- as.data.frame(model.matrix(formula, data = data))[, 
            -1, drop = FALSE]
        mm <- ncol(X)
    }
    if (!missing(centered)) {
        if (mm != 0 & centered == TRUE) 
            X <- apply(as.matrix(X), 2, function(x) x - mean(x))
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
    
    n <- sum(keep)
    data <- data.frame(start = start, Y = Y, stat = status, R)
    if (mm != 0) 
        data <- cbind(data, X)
    out <- list(data = data, R = R, status = status, start = start, 
        Y = Y, X = X, m = mm, n = n, type = type, Y.surv = Y.surv, 
        Terms = Terms, ratetable = ratetable, offset = offset)
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
    fit <- lik.fit(data, rform$m, length(interval[-1]), init, 
        control, offset)
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
    if (maxiter > 1 & fit$nit >= maxiter) 
        warning("Ran out of iterations and did not converge")
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
            ",", interval[2], ")", sep = "")
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
            k = k, data.rest,g=X$grupa[1])
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
            ncol = nfk, byrow = F, nrow = nrow(Z)) * matrix(fk, 
            ncol = nfk, byrow = T, nrow = nrow(Z))
    X <- X[X$start < int[2], ]
    X$fin <- (X$Y <= int[2])
    X$event <- X$fin * X$stat
    ford <- eval(substitute(paste("[", a, ",", b, ")", sep = ""), 
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
            ford <- c(ford, eval(substitute(paste("(", a, ",", 
                b, "]", sep = ""), list(a = meje[l], b = meje[l + 
                1]))))
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
        if (is.null(rform$X)) 
            local.ht <- glm(cbind(nd, ld - nd) ~ -1 + fu + offset(log(k)), 
                data = grm1, family = ht)
        else {
            xmat <- as.matrix(grm1[, 7:(ncol(grm1) - 2)])
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
        if (is.null(rform$X)) 
            local.ht <- glm(nd ~ -1 + fu, data = grm1, family = pot, 
                offset = lny)
        else {
            xmat <- as.matrix(grm1[, 7:(ncol(grm1) - 2)])
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
    
    m <- ncol(data)
    rem <- m - nfk - 3
    interval <- object$int
    int <- ceiling(max(interval))
    R <- data[, 4:(nfk + 3)]
    lp <- matrix(-log(srvxp.fit(as.matrix(R), rep(365.24, n), 
        object$ratetable))/scale, ncol = 1)
    fu <- NULL
   
   
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
    temp <- survfit(fit$y, type = "kaplan-meier")
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
            usable <- which(varr[i, i, ] > 0)
        else usable <- which(varr != 0)
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
        usable <- which(sumdiag > 0)
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
    times <- time[stat == 1]
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
        if (length(invV) == length(varr)) 
            usable <- rep(TRUE, dim(varr)[3])
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
    xx <- x$x/scale
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
            xtime <- as.numeric(dimnames(yy)[[1]])
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
}




invtime <- function (y = 0.1, age = 23011, sex = "male", year = 9497, scale = 1, 
    ratetable = slop, lower, upper) 
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
    if(length(int)!=1)int <- max(int)
    data <- rform$data
    if (rform$m == 0) {
        if (rform$type == "counting") 
            fit <- coxph(Surv(start, stop, stat) ~ 1, data = data, 
                init = init, control = control, x = TRUE, ...)
        else fit <- coxph(Surv(stop, stat) ~ 1, data = data, 
            init = init, control = control, x = TRUE, ...)
    }
    else {
        xmat <- as.matrix(data[, (4 + nfk):ncol(data)])
        fit <- coxph(Surv(start, stop, stat) ~ xmat, data = data, 
            init = init, control = control, x = TRUE, ...)
        names(fit[[1]]) <- names(rform$X)
    }
    fit$call <- match.call()
    if (length(rform$na.action)) 
        fit$na.action <- rform$na.action
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
        factor = c(0, 1, 0), cutpoints = list((0:(dimi[1] - 1)) * 
            (365.24), NULL, cp), class = "ratetable")
    temp
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
    if (length(x$warnme)) 
        cat("\n", x$warnme, "\n\n")
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

