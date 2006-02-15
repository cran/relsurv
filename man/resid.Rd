\name{residuals.rsadd}
\alias{residuals.rsadd}
\alias{resid.com}


\title{Calculate Residuals for a rsadd Fit}
\description{
Calculates partial residuals for an additive relative survival model.
}
\usage{
   \method{residuals}{rsadd}(object,type="schoenfeld",...)
}

\arguments{
\item{object}{
an object inheriting from class \code{rsadd}, representing a
          fitted additive relative survival model. Typically this is the output
          from the \code{rsadd} function. 

}
\item{type}{ character string indicating the type of residual desired. 
Currently only Schoenfeld residuals are implemented. }
\item{...}{other arguments.}
}

\value{
A list of the following values is returned:
\item{res}{a matrix containing the residuals for each variable.}
\item{varr}{the variance for each residual}
\item{varr1}{the sum of \code{varr}.}
\item{kvarr}{the derivative of each residual, to be used in \code{rs.zph} function.}
\item{kvarr1}{the sum of \code{kvarr}.}
}

\examples{
data(slopop)
data(rdata)
fit <- rsadd(Surv(time,cens)~sex+ratetable(age=age*365,sex=sex,year=year),
       ratetable=slopop,data=rdata,int=5)
sresid <- residuals.rsadd(fit)
}

\seealso{\code{\link{rsadd}}.}

\keyword{survival}



