\name{residuals.rsadd}
\alias{residuals.rsadd}
\alias{resid.com}


\title{Calculate Residuals for a "rsadd" Fit}
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

\references{
Package. Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.

Goodness of fit: Stare J.,Pohar Perme M., Henderson R. (2005) "Goodness of fit of relative survival models." Statistics in Medicine, \bold{24}: 3911--3925.
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



