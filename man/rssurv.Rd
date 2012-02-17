\name{rs.surv}
\alias{rs.surv}
\alias{rs.survo}


\title{Compute a Relative Survival Curve}

\description{
     Computes an estimate of the relative survival curve using the
     Ederer II method, Pohar-Perme method or the Hakulinen method
}
\usage{
     rs.surv(formula, data,ratetable=survexp.us,na.action,fin.date,
	     method="pohar-perme",conf.type="log",conf.int=0.95,type="kaplan-meier") 
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


\item{fin.date}{the date of the study ending, used for calculating the
          potential follow-up times in the Hakulinen method. If
          missing, it is calculated as \code{max(year+time)}. 
}
\item{method}{the method for calculating the relative survival. The options are \code{pohar-perme}(default), \code{ederer} and \code{hakulinen}.}
\item{conf.type}{one of \code{plain}, \code{log} (the default), or \code{log-log}. The first option causes the standard intervals curve +- k *se(curve), where k is determined from conf.int. The
log option calculates intervals based on the cumulative hazard or log(survival). The last option bases intervals on the log hazard or log(-log(survival)). }
\item{conf.int}{the level for a two-sided confidence interval on the survival curve(s). Default is 0.95.}
\item{type}{defines how survival estimates are to be calculated given the hazards.  The default (\code{kaplan-meier})
calculates the product integral, whereas the option \code{fleming-harrington} exponentiates the negative cumulative hazard. Analogous to the usage in \code{survfit}.
, \code{ederer} and \code{hakulinen}.
}

}

\details{
NOTE: All times used in the formula argument must be specified in days. This is true for the follow-up time as well as for
any variables needed \code{ratetable} object, like \code{age} and \code{year}. On the contrary, the \code{int} argument requires
interval specification in years.

The potential censoring times needed for the calculation of the expected survival by the Hakulinen method
are calculated automatically. The times of censoring are left as they are, the times of events are replaced with 
\code{fin.date - year}.

The calculation of the Pohar-Perme estimate is more time consuming since more data are needed from the population tables. 
The old version of the function, now named \code{rs.survo} can be used as a faster version for the Hakulinen and Ederer II estimate.
}

\value{
 a \code{survfit} object; see the help on \code{survfit.object} for details.
 The \code{survfit} methods are used for \code{print},
     \code{plot}, \code{lines}, and \code{points}.

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
rs.surv(Surv(time,cens)~sex+ratetable(age=age*365,sex=sex,
      year=year),ratetable=slopop,data=rdata)
}


\seealso{
\code{survfit},
\code{survexp}
}

\keyword{survival}
