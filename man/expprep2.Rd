\name{expprep2}
\alias{expprep2}



\title{expprep2 function}

\description{
     Helper calculation function using C code. Saved also as exp.prep (unexported function).
}
\usage{
     expprep2(x, y,ratetable,status,times,fast=FALSE,ys,prec,cmp=F,netweiDM=FALSE)
}
\arguments{
\item{x}{ matrix of demographic covariates - each individual has one line}
\item{y}{ follow-up time for each individual (same length as nrow(x))}
\item{ratetable}{ rate table used for calculation}
\item{status}{ status for each individual (same length as nrow(x)!), not needed if we only need Spi, status needed for rs.surv}
\item{times}{ times at which we wish to evaluate the quantities, not needed if we only need Spi, times needed for rs.surv}
\item{fast}{ for mpp method only}
\item{ys}{ entry times (if empty, individuals are followed from time 0)}
\item{prec}{ deprecated}
\item{cmp}{ should cmpfast.C be used}
\item{netweiDM}{ should new netwei script be used}
}

\details{
Helper function used in rs.surv and other relsurv functions.
}

\value{
List containing the calculated hazards and probabilities using the population mortality tables.
}

\keyword{survival}

\seealso{rs.surv}
