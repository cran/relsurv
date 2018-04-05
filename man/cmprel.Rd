\name{cmp.rel}
\alias{cmp.rel}
\alias{print.cmp.rel}



\title{Compute crude probability of death}

\description{
     Estimates the crude probability of death due to disease and due to population reasons
}
\usage{
     cmp.rel(formula, data, ratetable = relsurv::slopop, na.action,tau,
	     conf.int=0.95,precision=1,add.times) 
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
          such as \code{slopop}.  
}
\item{na.action}{a missing-data filter function, applied to the model.frame,
          after any subset argument has been used.  Default is
          \code{options()$na.action}.} 


\item{tau}{the maximum follow-up time of interest, all times larger than \code{tau} shall be censored. Equals maximum observed time by default 
}
\item{conf.int}{the level for a two-sided confidence interval on the survival curve(s). Default is 0.95.}
\item{precision}{the level of precision used in the numerical integration of variance. Default is 1, which means that daily intervals are taken, the value may be decreased to get a higher precision or increased to achieve a faster calculation. The calculation intervals always include at least all times of event and censoring as border points.}
\item{add.times}{specific times at which the value of estimator and its variance should be evaluated. Default is all the event and censoring times.}
}

\details{
NOTE: All times used in the formula argument must be specified in days. This is true for the follow-up time as well as for
any variables needed \code{ratetable} object, like \code{age} and \code{year}. 

Note that numerical integration is required to calculate the variance estimator. The integration precision is set with argument \code{precision}, which defaults to daily intervals, a default that should give enough precision for any practical purpose.

The area under the curve is calculated on the interval [0,\code{tau}]. 

Function \code{summary} may be used to get the output at specific points in time.

}

\value{
An object of class \code{cmp.rel}.
Objects of this class have methods for the functions \code{print} and \code{plot}. The \code{summary} function can be used for printing output at required time points.
An object of class \code{cmp.rel} is composed of several lists, each pertaining the cumulative hazard function for one risk and one strata. Each of the lists contains the following objects:
\item{time}{the time-points at which the curves are estimated}
\item{est}{the estimate}
\item{var}{the variance of the estimate}
\item{lower}{the lower limit of the confidence interval}
\item{upper}{the upper limit of the confidence interval}
\item{area}{the area under the curve calculated on the interval [0,\code{tau}]}
\item{index}{indicator of event and censoring times among all the times in the output. The times added via paramater \code{add.times} are also included}
\item{add.times}{the times added via parameter \code{add.times}}
}




\examples{
data(slopop)
data(rdata)
#calculate the crude probability of death
#note that the variable year is in the date format - the number of days since 01.01.1960 and that 
#age must be multiplied by 365 in order to be expressed in days.
fit <- cmp.rel(Surv(time,cens)~sex+ratetable(age=age*365,sex=sex,
      year=year),ratetable=slopop,data=rdata,tau=3652.41)
fit
plot(fit,col=c(1,1,2,2),xscale=365.241,xlab="Time (years)")
}


\seealso{
\code{rs.surv}, \code{summary.cmp.rel}
}

\keyword{survival}
