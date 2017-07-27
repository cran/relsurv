\name{rs.diff}
\alias{rs.diff}
\alias{print.rsdiff}




\title{Test Net Survival Curve Differences}

\description{
     Tests if there is a difference between two or more net survival curves using a log-rank type test.
}
\usage{
     rs.diff(formula, data, ratetable = relsurv::slopop, na.action,precision=1) 
}
\arguments{
\item{formula}{
a formula object, with the response as a \code{Surv} object on the
          left of a \code{~} operator, and, if desired, terms separated  by
          the \code{+} operator on the right. If the variables are not
          organized and named in the same way as in the population
          tables, a \code{ratetable} term must be added to match each
          subject  to his/her expected cohort.   
A formula expression as for other survival models, of the form \code{Surv(time, status) ~ predictors}.
Each combination of predictor values defines a subgroup. A \code{strata} term may be used to produce a stratified test. 

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


\item{precision}{Precision for numerical integration. Default is 1, which means that daily intervals are taken, the value may be decreased to get a higher precision or increased to achieve a faster calculation. The calculation intervals always include at least all times of event and censoring as border points.
}

}

\details{
NOTE: All times used in the formula argument must be specified in days. This is true for the follow-up time as well as for
any variables needed \code{ratetable} object, like \code{age} and \code{year}.
}

\value{
 a \code{rsdiff} object; can be printed with \code{print}.

}

\references{
Graffeo, N., Castell, F., Belot, A. and Giorgi, R. (2016) "A log-rank-type test to compare net survival distributions. Biometrics. doi: 10.1111/biom.12477"
Pavlic, K., Pohar Perme, M. (2017) "On comparison of net survival curves. BMC Med Res Meth. doi: 10.1186/s12874-017-0351-3"
}




\examples{
data(slopop)
data(rdata)
#calculate the relative survival curve
#note that the variable year is given in days since 01.01.1960 and that 
#age must be multiplied by 365.241 in order to be expressed in days.
rs.diff(Surv(time,cens)~sex+ratetable(age=age*365.241,sex=sex,
      year=year),ratetable=slopop,data=rdata)
}


\seealso{
\code{rs.surv},
\code{survdiff}
}

\keyword{survival}
