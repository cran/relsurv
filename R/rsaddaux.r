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
    if(length(x$warnme))cat("\n",x$warnme, "\n\n")
    invisible(x)
}

summary.rsadd <- function (object, correlation = FALSE, symbolic.cor = FALSE, 
    ...) 
{   if(inherits(object,"glm")){
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
		dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
		    "t value", "Pr(>|t|)"))
		covmat.unscaled <- covmat <- matrix(, 0, 0)
		aliased <- is.na(coef(object))
		df.f <- length(aliased)
	    }
	    ans <- c(object[c("call", "terms", "family","iter","warnme")], list(coefficients = coef.table, var = covmat,aliased=aliased))
	    if (correlation && p > 0) {
		dd <- s.err
		ans$correlation <- covmat/outer(dd, dd)
		ans$symbolic.cor <- symbolic.cor
	    }
	    class(ans) <- "summary.rsadd"
    }
    else if (inherits(object,"rsadd")){
    	aliased <- is.na(coef(object))
    	coef.p <- object$coef
        var.cf <- diag(object$var)
	s.err <- sqrt(var.cf)
	tvalue <- coef.p/s.err
	dn <- c("Estimate", "Std. Error")
	pvalue <- 2 * pnorm(-abs(tvalue))
	coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
	dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", "Pr(>|z|)"))
	ans <- c(object[c("call", "terms","iter","var")], list(coefficients = coef.table,aliased=aliased))
	if (correlation && sum(aliased)!= length(aliased)) {

			dd <- s.err
			ans$correlation <- object$var/outer(dd, dd)
			ans$symbolic.cor <- symbolic.cor
	    }
	class(ans) <- "summary.rsadd"
    }
    else ans <- object
    return(ans)
}


print.summary.rsadd <- {function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
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
    if(length(x$warnme))cat("\n",x$warnme,"\n")
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
}

