\name{rsadd}
\alias{rsadd}
\alias{print.rsadd}
\alias{summary.rsadd}
\alias{print.summary.rsadd}
\alias{maxlik}
\alias{glmxp}
\alias{lik.fit}


\title{Fit an Additive model for Relative Survival}
\description{
The function fits an additive model to the data. The methods implemented are the maximum likelihood method,
a glm model with a \code{binomial} error and a glm model with a \code{poisson} error. 
}
\usage{
rsadd(formula, data=parent.frame(), ratetable = survexp.us,
      int, na.action, method, init,centered,control,...)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a \code{~} operator, and
the terms on the right. The terms consist of predictor variables separated by
the \code{+} operator, along with a \code{ratetable} term. The \code{ratetable} term
matches each subject to his/her expected cohort. If the variables are organized and named
in the same way as in the population tables, the \code{ratetable} term can be omitted.
The response must be a survival object as returned by the \code{Surv} function. The time must be in days.
}
\item{data}{
a data.frame in which to interpret the variables named in
the \code{formula}.
}
\item{ratetable}{
a table of event rates, organized as a \code{ratetable} object, such as \code{survexp.us}. 
}
\item{int}{
either a single value denoting the number of follow-up years or a vector
specifying the intervals (in years) in which the hazard is constant (the times that are
bigger than \code{max(int)} are censored. If missing, the intervals are set to be one year 
long and include the maximum observed follow-up time.
}
\item{na.action}{a missing-data filter function, applied to the model.frame,
          after any subset argument has been used.  Default is
          \code{options()$na.action}. }
\item{method}{\code{glm.bin} or \code{glm.poi} for a glm model and \code{max.lik} for the maximum likelihood model (default). }
\item{init}{vector of initial values of the iteration.  Default initial
          value is zero for all variables. }

\item{centered}{if \code{TRUE}, all the variables are centered before fitting and the baseline excess hazard
	is calculated accordingly. Default is \code{FALSE}. }

\item{control}{a list of parameters for controlling the fitting process. 
          See the documentation for \code{glm.control} for details.
}
\item{...}{other arguments will be passed to \code{glm.control}.}
}

\Details{
The methods using glm are methods for grouped data. The groups are formed according to the covariate values.
This should be taken into account when fitting a model.
}

\value{
An object of class \code{rsadd}. In the case of \code{method="glm.bin"} and \code{method="glm.poi"} the class also
inherits from \code{glm} which  inherits from the class \code{lm}.
Objects of this class have methods for the functions \code{print} and \code{summary}.
An object of class \code{rsadd} is a list containing at least the following components:
\item{data}{the original data, along with the variables defined in the rate table}
\item{warnme}{the number of groups with the observed number of deaths is smaller than the expected.}
\item{ratetable}{the ratetable used.}
\item{int}{the maximum time (in years) used. All the events at and after this value are censored.}
\item{method}{the fitting method that was used.}
\item{linear.predictors}{the vector of linear predictors, one per subject.}
Additionally the glm models include:
\item{groups}{life tables for groups specified by the covariates.}
}

\examples{
data(slopop)
data(rdata)
rsadd(Surv(time,cens)~sex+as.factor(agegr)+ratetable(age=age*365,sex=sex,year=year),
      ratetable=slopop,data=rdata,int=5)
}


\references{
Dickman P.W., Sloggett A., Hills M., Hakulinen
T. "Regression models for relative survival." \emph{Statistics
in Medicine,} 2004.
}



\seealso{
\code{\link{rstrans}},
\code{\link{rsmul}}
}

\keyword{survival}
