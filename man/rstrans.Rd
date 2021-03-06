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
rstrans(formula, data, ratetable, int,na.action,init,control,rmap,...)
}

\arguments{
\item{formula}{
a formula object, with the response as a \code{Surv} object on the
          left of a \code{~} operator, and, if desired, terms separated  by
          the \code{+} operator on the right. 

 NOTE: The follow-up time must be in days.
}
\item{data}{
a data.frame in which to interpret the variables named
 in the \code{formula}.  
}
\item{ratetable}{
a table of event rates, such as \code{slopop}. 
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
\item{rmap}{an optional list to be used if the variables are not
          organized and named in the same way as in the \code{ratetable} object. See details below.}

\item{...}{other arguments will be passed to \code{coxph.control}.}

}

\value{
 an object of class \code{coxph}. See \code{coxph.object} and \code{coxph.detail} for details.  
\item{y}{
 an object of class \code{Surv} containing the transformed times (these times do not depend on covariates). 
}
}


\details{
NOTE: The follow-up time must be specified in days. The \code{ratetable} being used may have different variable names and formats than the user's data set, this is dealt with by the \code{rmap} argument. For example, if age is in years in the data set but in days in the \code{ratetable} object, age=age*365.241 should be used. The calendar year can be in any date format (date, Date and POSIXt are allowed), the date formats in the \code{ratetable} and in the data may differ. 
A side product of this function are the transformed times - stored in teh \code{y} object of the output. To get these times, covariates are of course irrelevant.
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
#age must be multiplied by 365.241 in order to be expressed in days.
fit <- rstrans(Surv(time,cens)~sex+as.factor(agegr),rmap=list(age=age*365.241,
        sex=sex,year=year),ratetable=slopop,data=rdata)


#check the goodness of fit
rs.br(fit)
}






\seealso{\code{\link{rsmul}}, \code{\link{invtime}}, \code{\link{rsadd}}, \code{\link{survexp}}.}

\keyword{survival}
