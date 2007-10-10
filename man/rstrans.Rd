\name{rstrans}
\alias{rstrans}

\title{Fit Cox Proportional Hazards Model in Transformed Time}
\description{
The function transforms each person's time to his/her probability
of dying at that time according to the ratetable. It then fits the
Cox proportional hazards model with the transformed times as a 
response. A possible side effect are the transformed times.  
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
  The response must be a survival object as
returned by the \code{Surv} function. The time must be in days.  
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
 an object of class \code{Surv} containing the transformed times.
}
}



\references{
Stare J., Henderson R., Pohar M. "An individual measure for relative survival." \emph{Journal
of the Royal Statistical Society: Series C} , 54 115-126, 2005.

Pohar M., Stare J. "Relative survival analysis in R." \emph{Computer Methods and Programs in Biomedicine,} 81: 272-278, 2006.
}


\examples{
data(slopop)
data(rdata)
fit <- rstrans(Surv(time,cens)~sex+as.factor(agegr)+ratetable(age=age*365,
        sex=sex,year=year),ratetable=slopop,data=rdata)


#check the goodness of fit
rs.br(fit)
}






\seealso{\code{\link{rsmul}}, \code{\link{invtime}}, \code{\link{rsadd}}, \code{\link{survexp}}.}

\keyword{survival}
