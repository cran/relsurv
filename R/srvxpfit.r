srvxp.fit<-function(x, y,ratetable)
 {
    x <- cbind(1:nrow(x),as.matrix(x))	
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
