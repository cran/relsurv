\name{rs.br}
\alias{rs.br}
\alias{plot.rs.br}
\alias{print.rs.br}


\title{Test the Proportional Hazards Assumption for Relative Survival
 Regression Models}
\description{
Test the proportional hazards assumption for relative survival  models  (\code{rsadd}, \code{rsmul} or \code{rstrans})
by forming a Brownian Bridge. }
\usage{
rs.br(fit,sc,rho=0,test="max",global=TRUE)
}

\arguments{
\item{fit}{
the result of fitting a relative survival model, using the
          \code{rsadd}, \code{rsmul} or \code{rstrans} function. 
}
\item{sc}{
partial residuals calculated by the \code{resid} function. This is used to save time if several tests
are to be calculated on these residuals and can otherwise be omitted.
}
\item{rho}{
a number controlling the weigths of residuals. The weights are the number of individuals at
risk at each event time to the power \code{rho}. The default is \code{rho=0}, which sets all
weigths to 1.
}
\item{test}{ 
a character string specifying the test to be performed on Brownian bridge.
Possible values are \code{"max"} (default), which tests the maximum absolute value of the 
bridge, and \code{cvm}, which calculates the Cramer Von Mises statistic.
}
\item{global}{should a global Brownian bridge test be performed, in addition to the
          per-variable tests}

}

\value{
 an object of class \code{rs.br}. This function would usually be followed by both a
print and a plot of the result. The plot gives a Brownian bridge for each of the variables.
The horizontal lines are the 95% and 99% confidence intervals for the maximum absolute value
of the Brownian bridge}




\references{
\item{Goodness of fit}{
 Stare J.,Pohar Perme M., Henderson R. "Goodness of fit of relative survival models." \emph{Statistics in Medicine}, 24: 3911-3925, 2005. 
}

\item{package}{Pohar M., Stare J. "Relative survival analysis in R." \emph{Computer Methods and Programs in Biomedicine,} 81: 272-278, 2006.}

\item{relative survival}{Pohar, M., Stare, J. "Making relative survival analysis relatively easy."
\emph{Computers in biology and medicine}, 37: 1741-1749, 2007.}
}



\examples{
data(slopop)
data(rdata)
fit <- rsadd(Surv(time,cens)~sex+ratetable(age=age*365,sex=sex,
	     year=year),ratetable=slopop,data=rdata,int=5)
rsbr <- rs.br(fit)
print(rsbr)
plot(rsbr)
}

\seealso{\code{\link{rsadd}}, \code{rsmul}, \code{rstrans}, \code{\link{resid}}.}

\keyword{survival}
