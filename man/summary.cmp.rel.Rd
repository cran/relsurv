\name{summary.cmp.rel}
\alias{summary.cmp.rel}


\title{Summary of the crude probability of death}

\description{
     Returns a list containing the estimated values at required times. 
}
\usage{
	\method{summary}{cmp.rel}(object, times, scale = 365.241,area=FALSE,...) 

}
\arguments{
\item{object}{output of the function \code{cmp.rel}.}
\item{times}{the times at which the output is required.}
\item{scale}{The time scale in which the times are specified. The default value is \code{1}, i.e. days.}
\item{area}{Should area under the curves at time \code{tau} be printed out? Default is \code{FALSE}.}
\item{...}{Additional arguments, currently not implemented}
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
#calculate the crude probability of death and summarize it
fit <- cmp.rel(Surv(time,cens)~sex,rmap=list(age=age*365),
      ratetable=slopop,data=rdata,tau=3652.41)
summary(fit,c(1,3),scale=365.241)
}


\seealso{
\code{cmp.rel}
}

\keyword{survival}
