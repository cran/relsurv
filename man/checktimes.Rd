\name{checktimes}
\alias{checktimes}


\title{Summary of a cumulative incidence curve}

\description{
     Returns a list containing the estimated values at required times. 
}
\usage{
	checktimes(w, times, xscale = 1,area=FALSE) 

}
\arguments{
\item{w}{output of the function \code{cmp.rel}.}
\item{times}{the times at which the output is required.}
\item{xscale}{The time scale in which the times are specified. The default value is \code{1}, i.e. days.}
\item{area}{Should area under the curves at time \code{tau} be printed out? Default is \code{FALSE}.}
}

\details{
The variance is calculated using numerical integration. If the required time is not a time at which the value was estimated, the value at the last time before it is reported. The density of the time points is set by the \code{precision} argument in the \code{cmp.rel} function.
}
\value{
A list of values is returned.
}


\examples{
data(slopop)
data(rdata)
#calculate the cumulative incidence curve
#note that the variable year is given in days since 01.01.1960 and that 
#age must be multiplied by 365 in order to be expressed in days.
fit <- cmp.rel(Surv(time,cens)~sex+ratetable(age=age*365,sex=sex,
      year=year),ratetable=slopop,data=rdata,tau=3652.41)
checktimes(fit,c(1,3),xscale=365.241)
}


\seealso{
\code{cmp.rel}
}

\keyword{survival}
