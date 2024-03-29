% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmprel.r
\name{cmp.rel}
\alias{cmp.rel}
\alias{print.cmp.rel}
\title{Compute crude probability of death}
\usage{
cmp.rel(
  formula = formula(data),
  data = parent.frame(),
  ratetable = relsurv::slopop,
  na.action,
  tau,
  conf.int = 0.95,
  precision = 1,
  add.times,
  rmap
)
}
\arguments{
\item{formula}{a formula object, with the response as a \code{Surv} object
on the left of a \code{~} operator, and, if desired, terms separated by the
\code{+} operator on the right. If no strata are used, \code{~1} should be
specified.

NOTE: The follow-up time must be in days.}

\item{data}{a data.frame in which to interpret the variables named in the
\code{formula}.}

\item{ratetable}{a table of event rates, organized as a \code{ratetable}
object, such as \code{slopop}.}

\item{na.action}{a missing-data filter function, applied to the model.frame,
after any subset argument has been used.  Default is
\code{options()$na.action}.}

\item{tau}{the maximum follow-up time of interest, all times larger than
\code{tau} shall be censored. Equals maximum observed time by default}

\item{conf.int}{the level for a two-sided confidence interval on the
survival curve(s). Default is 0.95.}

\item{precision}{the level of precision used in the numerical integration of
variance. Default is 1, which means that daily intervals are taken, the
value may be decreased to get a higher precision or increased to achieve a
faster calculation. The calculation intervals always include at least all
times of event and censoring as border points.}

\item{add.times}{specific times at which the value of estimator and its
variance should be evaluated. Default is all the event and censoring times.}

\item{rmap}{an optional list to be used if the variables are not organized
and named in the same way as in the \code{ratetable} object. See details
below.}
}
\value{
An object of class \code{cmp.rel}. Objects of this class have
methods for the functions \code{print} and \code{plot}. The \code{summary}
function can be used for printing output at required time points. An object
of class \code{cmp.rel} is composed of several lists, each pertaining the
cumulative hazard function for one risk and one strata. Each of the lists
contains the following objects: \item{time}{the time-points at which the
curves are estimated} \item{est}{the estimate} \item{var}{the variance of
the estimate} \item{lower}{the lower limit of the confidence interval}
\item{upper}{the upper limit of the confidence interval} \item{area}{the
area under the curve calculated on the interval [0,\code{tau}]}
\item{index}{indicator of event and censoring times among all the times in
the output. The times added via paramater \code{add.times} are also
included} \item{add.times}{the times added via parameter \code{add.times}}
}
\description{
Estimates the crude probability of death due to disease and due to
population reasons
}
\details{
NOTE: The follow-up time must be specified in days. The \code{ratetable}
being used may have different variable names and formats than the user's
data set, this is dealt with by the \code{rmap} argument. For example, if
age is in years in the data set but in days in the \code{ratetable} object,
age=age*365.241 should be used. The calendar year can be in any date format
(date, Date and POSIXt are allowed), the date formats in the
\code{ratetable} and in the data may differ.

Note that numerical integration is required to calculate the variance
estimator. The integration precision is set with argument \code{precision},
which defaults to daily intervals, a default that should give enough
precision for any practical purpose.

The area under the curve is calculated on the interval [0,\code{tau}].

Function \code{summary} may be used to get the output at specific points in
time.
}
\examples{


data(slopop)
data(rdata)
#calculate the crude probability of death
#note that the variable year must be given in a date format and that 
#age must be multiplied by 365.241 in order to be expressed in days.
fit <- cmp.rel(Surv(time,cens)~sex,rmap=list(age=age*365.241),
		ratetable=slopop,data=rdata,tau=3652.41)
fit
plot(fit,col=c(1,1,2,2),xscale=365.241,xlab="Time (years)")
#if no strata are desired:
fit <- cmp.rel(Surv(time,cens)~1,rmap=list(age=age*365.241),
		ratetable=slopop,data=rdata,tau=3652.41)



}
\references{
Package: Pohar Perme, M., Pavlic, K. (2018) "Nonparametric
Relative Survival Analysis with the R Package relsurv". Journal of
Statistical Software. 87(8), 1-27, doi: "10.18637/jss.v087.i08"
}
\seealso{
\code{rs.surv}, \code{summary.cmp.rel}
}
\keyword{survival}
