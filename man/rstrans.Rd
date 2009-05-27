\name{rstrans}
\alias{rstrans}

\title{Fit Cox Proportional Hazards Model in Transformed Time}
\description{
The function transforms each person's time to his/her probability
of dying at that time according to the ratetable. It then fits the
Cox proportional hazards model with the transformed times as a 
response. It can also be used for calculatin the transformed times (no covariates are needed in the formula for that purpose).  
}
\usage{
rstrans(formula, data, ratetable, int,na.action,init,control,...)
}

\arguments{
\item{formula}{
a formula object, with the response on the left of a \code{~} operator, and
the terms on the right. The terms consist of predictor variables separated by
the \code{+} operator, along with a \code{ratetable} term. The \code{ratetable} term
matches each subject to his/her expected cohort.If the variables are organized and named
in the same way as in the population tables, the \code{ratetable} term can be omitted.
 The response must be a survival object as returned by the \code{Surv} function.
 
NOTE: The time must be in days, and the same is required for the ratetable 
variables (the variables used in the population tables), for example age and year (year must be 
given in the \code{date} format, i.e. in number of days since 01.01.1960). 
}
\item{data}{
a data.frame in which to interpret the variables named
 in the \code{formula}.  
}
\item{ratetable}{
a table of event rates, such as \code{survexp.us}. 
}
\item{int}{ 
the number of follow-up years used for calculating survival(the rest is censored). If missing, it is set the 
the maximum observed follow-up time.
}
\item{na.action}{a missing-data filter function, applied to the model.frame,
          after any subset argument has been used.  Default is
          \code{options()$na.action}. }

\item{init}{vector of initial values of the iteration.  Default initial
          value is zero for all variables. }
\item{control}{a list of parameters for controlling the fitting process. 
          See the documentation for \code{coxph.control} for details.
}
\item{...}{other arguments will be passed to \code{coxph.control}.}

}

\value{
 an object of class \code{coxph}. See \code{coxph.object} and \code{coxph.detail} for details.  
\item{y}{
 an object of class \code{Surv} containing the transformed times (these times do not depend on covariates). 
}
}


\details{
NOTE: All times used in the formula argument must be specified in days. This is true for the follow-up time as well as for
any variables needed \code{ratetable} object, like \code{age} and \code{year}. On the contrary, the \code{int} argument requires
interval specification in years.
}


\references{
Method: Stare J., Henderson R., Pohar M. (2005) "An individual measure for relative survival." Journal
of the Royal Statistical Society: Series C, \bold{54} 115--126.

Package. Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.
}


\examples{
data(slopop)
data(rdata)

#fit a Cox model using the transformed times
#note that the variable year is given in days since 01.01.1960 and that 
#age must be multiplied by 365 in order to be expressed in days.
fit <- rstrans(Surv(time,cens)~sex+as.factor(agegr)+ratetable(age=age*365,
        sex=sex,year=year),ratetable=slopop,data=rdata)


#check the goodness of fit
rs.br(fit)
}






\seealso{\code{\link{rsmul}}, \code{\link{invtime}}, \code{\link{rsadd}}, \code{\link{survexp}}.}

\keyword{survival}
