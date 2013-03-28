\name{cmp.rel}
\alias{cmp.rel}
\alias{plot.cmp.rel}
\alias{print.cmp.rel}



\title{Compute crude mortality curves}

\description{
     Estimates the probability of dying due to disease and due to population reasons
}
\usage{
     cmp.rel(formula, data, ratetable = survexp.us, na.action,tau,conf.int=0.95) 

}
\arguments{
\item{formula}{
a formula object, with the response as a \code{Surv} object on the
          left of a \code{~} operator, and, if desired, terms separated  by
          the \code{+} operator on the right. If the variables are not
          organized and named in the same way as in the population
          tables, a \code{ratetable} term must be added to match each
          subject  to his/her expected cohort.   For a single survival curve the \code{~ 1} part of the
          formula is not required. 

 NOTE: The time must be in days, and the same
is required for the ratetable variables (the variables used in the population tables), for example age and year (year must be 
given in the \code{date} format, i.e. in number of days since 01.01.1960). 
}
\item{data}{
a data.frame in which to interpret the variables named in the
          \code{formula}. 

}
\item{ratetable}{
a table of event rates, organized as a \code{ratetable} object,
          such as \code{survexp.us}.  
}
\item{na.action}{a missing-data filter function, applied to the model.frame,
          after any subset argument has been used.  Default is
          \code{options()$na.action}.} 


\item{tau}{the maximum follow-up time of interest, all times larger than \code{tau} shall be censored. Equals maximum observed time by default 
}
\item{conf.int}{the level for a two-sided confidence interval on the survival curve(s). Default is 0.95.}

}

\details{
NOTE: All times used in the formula argument must be specified in days. This is true for the follow-up time as well as for
any variables needed \code{ratetable} object, like \code{age} and \code{year}. On the contrary, the \code{int} argument requires
interval specification in years.

The area under the curve is calculated on the interval [0,\code{tau}]. 

}

\value{
An object of class \code{cmp.rel}.
Objects of this class have methods for the functions \code{print} and \code{plot}.
An object of class \code{cmp.rel} is composed of several lists, each pertaining the the cumulative hazar function for one risk and one strata. Each of the lists contains the following objects:
\item{time}{the time-points at which the curves are estimated}
\item{est}{the estimate}
\item{var}{the variance of the estimate}
\item{lower}{the lower limit of the confidence interval}
\item{upper}{the upper limit of the confidence interval}
\item{area}{the area under the curve calculated on the interval [0,\code{tau}]}
}



\references{
Package. Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.
}




\examples{
data(slopop)
data(rdata)
#calculate the relative survival curve
#note that the variable year is given in days since 01.01.1960 and that 
#age must be multiplied by 365 in order to be expressed in days.
fit <- cmp.rel(Surv(time,cens)~sex+ratetable(age=age*365,sex=sex,
      year=year),ratetable=slopop,data=rdata,tau=3652.4)
fit
plot(fit,col=c(1,1,2,2),xscale=365.24)
}


\seealso{
\code{rs.surv}
}

\keyword{survival}
