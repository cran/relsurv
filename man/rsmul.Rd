\name{rsmul}
\alias{rsmul}
\title{Fit Andersen et al Multiplicative Regression Model for Relative Survival}
\description{
  Fits the Andersen et al multiplicative regression model in relative survival. An extension of the 
coxph function using relative survival.
}
\usage{
rsmul(formula, data, ratetable = survexp.us, int,na.action,init,method,control,...)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a \code{~} operator, and
the terms on the right. The terms consist of predictor variables separated by
the \code{+} operator, along with a \code{ratetable} term. The \code{ratetable} term
matches each subject to his/her expected cohort. If the variables are organized and named
in the same way as in the population tables, the \code{ratetable} term can be omitted.
 The response must be a survival object as
returned by the \code{Surv} function. The time must be in days.
}
\item{data}{
a data.frame in which to interpret the variables named in
the \code{formula}.
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
\item{method}{the default method \code{mul} assumes hazard to be constant on yearly intervals. Method \code{mul1} uses the ratetable
to determine the time points when hazard changes. The \code{mul1} method is therefore more accurate, but at the same time can
be more computationally intensive.}
\item{control}{a list of parameters for controlling the fitting process. 
          See the documentation for \code{coxph.control} for details.
}
\item{...}{Other arguments will be passed to \code{coxph.control}.}
}



\value{
an object of class \code{coxph}.
}
\references{
Andersen, P.K., Borch-Johnsen, K., Deckert, T., Green, A.,
Hougaard, P., Keiding, N. and Kreiner, S. "A Cox regression
model for relative mortality and its application to diabetes
mellitus survival data.", \emph{Biometrics,} 921-932, 1985.
}

\examples{
data(slopop)
data(rdata)
rsmul(Surv(time,cens)~sex+as.factor(agegr)+ratetable(age=age*365.24,sex=sex,year=year),ratetable=slopop,data=rdata)
}




\seealso{\code{\link{rsadd}}, \code{\link{rstrans}}. }

\keyword{survival}


