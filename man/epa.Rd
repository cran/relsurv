\name{epa}
\alias{epa}

\title{Excess hazard function smoothing}
\description{
An Epanechnikov kernel function based smoother for smoothing the baseline excess hazard calculated by the \code{rsadd} function with the \code{EM} method.
}
\usage{
epa(fit,bwin,times,n.bwin=16,left=FALSE)
}
\arguments{
\item{fit}{
Fit from the additive relative survival model using the \code{EM} method.}
\item{bwin}{
The relative width of the smoothing window (default is 1).
}
\item{times}{
The times at which the smoother is to be evaluated. If missing, it is evaluated at all event times.
}
\item{n.bwin}{
Number of times that the window width may change.
}
\item{left}{If \code{FALSE} (default) smoothing is performed symmetrically, if \code{TRUE} only leftside neighbours are considered. }

}

\details{
The function performs Epanechnikov kernel smoothing. The follow up time is divided (according to percentiles of event times) into several intervals (number of intervals defined by \code{n.bwin}) in which the width is calculated as a factor of the maximum span between event times. 
Boundary effects are also taken into account on both sides.
}

\value{
A list with two components:
\item{lambda}{the smoothed excess baseline hazard function}
\item{times}{the times at which the smoothed excess baseline hazard is evaluated.}
}

\examples{
data(slopop)
data(rdata)
#fit an additive model with the EM method
fit <- rsadd(Surv(time,cens)~sex+age,rmap=list(age=age*365.241),
		ratetable=slopop,data=rdata,int=5,method="EM")
sm <- epa(fit)
plot(sm$times,sm$lambda)
}


\references{
Package. Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.


EM algorithm: Pohar Perme M., Henderson R., Stare, J. (2009) "An approach to estimation in relative survival regression." Biostatistics, \bold{10}: 136--146.

}



\seealso{
\code{\link{rsadd}},
}

\keyword{survival}
