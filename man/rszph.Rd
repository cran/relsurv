\name{rs.zph}
\alias{rs.zph}


\title{Behaviour of Covariates in Time for Relative Survival Regression Models}
\description{
Calculates the scaled partial residuals of a relative survival model (\code{rsadd}, \code{rsmul} or \code{rstrans})}
\usage{
rs.zph(fit,sc,transform="identity",var.type="sum")
}

\arguments{
\item{fit}{
the result of fitting an additive relative survival model, using the
          \code{rsadd}, \code{rsmul} or \code{rstrans}  function. 

	In the case of multiplicative and transformation models the output is 
	identical to \code{cox.zph} function, except no test is performed.
}
\item{sc}{
partial residuals calculated by the \code{resid} function. This is used to save time if several tests
are to be calculated on these residuals and can otherwise be omitted.
}
\item{transform}{
a character string specifying how the survival times should
be transformed. Possible values are \code{"km"}, \code{"rank"}, \code{"identity"} and \code{log}.
The default is \code{"identity"}.

}
\item{var.type}{ 
a character string specifying the variance used to scale the residuals. Possible values are
\code{"each"}, which estimates the variance for each residual separately, and \code{sum}(default), 
which assumes the same variance for all the residuals.
}

}

\value{
 an object of class \code{rs.zph}. This function would usually be followed by a
          plot of the result. The plot gives an estimate of the time-dependent
          coefficient \code{beta(t)}. If the proportional
          hazards assumption is true, \code{beta(t)} will be a horizontal
          line.  
}




\references{
Stare J., Pohar M., Henderson R. "Goodness of fit for relative survival models." \emph{Statistics
in Medicine} 24 3911-3925, 2005.
}


\examples{
data(slopop)
data(rdata)
fit <- rsadd(Surv(time,cens)~sex+ratetable(age=age*365,sex=sex,year=year),
       ratetable=slopop,data=rdata,int=5)
rszph <- rs.zph(fit)
}

\seealso{\code{\link{rsadd}}, \code{rsmul}, \code{rstrans}, \code{\link{resid}}, \code{\link{cox.zph}}.}

\keyword{survival}
