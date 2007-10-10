\name{rs.surv}
\alias{rs.surv}
\alias{rs.surv.rsadd}


\title{Compute a Relative Survival Curve}

\description{
     Computes an estimate of the relative survival curve using the
     Kaplan-Meier method for the observed and the Hakulinen method for
     the expected survival.
}
\usage{
     rs.surv(formula, data,ratetable=survexp.us,fin.date,method="hakulinen",...) }
\arguments{
\item{formula}{
a formula object, with the response as a \code{Surv} object on the
          left of a \code{~} operator, and, if desired, terms separated  by
          the \code{+} operator on the right. If the variables are not
          organized and named in the same way as in the population
          tables, a \code{ratetable} term must be added to match each
          subject  to his/her expected cohort.  The time must be in
          days.  For a single survival curve the \code{~ 1} part of the
          formula is not required. 
}
\item{data}{
a data.frame in which to interpret the variables named in the
          \code{formula}. 

}
\item{ratetable}{
a table of event rates, organized as a \code{ratetable} object,
          such as \code{survexp.us}.  
}
\item{fin.date}{ The date of the study ending, used for calculating the
          potential follow-up times in the Hakulinen method. If
          missing, it is calculated as \code{max(year+time)}. 
}
\item{method}{The method for calculating the expected survival. The options are \code{hakulinen}(default) and \code{conditional}, see \code{survexp} for details.}
\item{...}{
other arguments will be passed to the \code{survfit} function that
          calculates the observed survival.
}
}

\details{
The potential censoring times needed for the calculation of the expected survival by the Hakulinen method
are calculated automatically. The times of censoring are left as they are, the times of events are replaced with 
\code{fin.date - year}.
}

\value{
 a \code{survfit} object; see the help on \code{survfit.object} for details.
 The \code{survfit} methods are used for \code{print},
     \code{plot}, \code{lines}, and \code{points}.

}

\examples{
data(slopop)
data(rdata)
rs.surv(Surv(time,cens)~sex+ratetable(age=age*365,sex=sex,
      year=year),ratetable=slopop,data=rdata)
}


\seealso{
\code{survfit},
\code{survexp}
}

\keyword{survival}
