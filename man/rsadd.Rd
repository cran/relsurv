\name{rsadd}
\alias{rsadd}

\title{Fit an Additive model for Relative Survival}
\description{
The function fits an additive model to the data. The methods implemented are the maximum likelihood method, the semiparametric
method, a glm model with a \code{binomial} error and a glm model with a \code{poisson} error. 
}
\usage{
rsadd(formula, data=parent.frame(), ratetable = relsurv::slopop,
      int, na.action, method, init,bwin,centered,cause,control,rmap,...)
}
\arguments{
\item{formula}{
a formula object, with the response as a \code{Surv} object on the
          left of a \code{~} operator, and, if desired, terms separated  by
          the \code{+} operator on the right. \code{Surv(start,stop,event)} 
	    outcomes are also possible for time-dependent covariates and 
	    left-truncation for \code{method='EM'}.

 NOTE: The follow-up time must be in days.
}
\item{data}{
a data.frame in which to interpret the variables named in
the \code{formula}.
}
\item{ratetable}{
a table of event rates, organized as a \code{ratetable} object, such as \code{slopop}. 
}
\item{int}{
either a single value denoting the number of follow-up years or a vector
specifying the intervals (in years) in which the hazard is constant (the times that are
bigger than \code{max(int)} are censored. If missing, only one interval (from time 0 to maximum observation time) is assumed.
 The EM method does not need the intervals,
only the maximum time can be specified (all times are censored after this time point).
}
\item{na.action}{a missing-data filter function, applied to the model.frame,
          after any subset argument has been used.  Default is
          \code{options()$na.action}. }
\item{method}{\code{glm.bin} or \code{glm.poi} for a glm model, \code{EM} for the EM algorithm and \code{max.lik} for the maximum likelihood model (default). }
\item{init}{vector of initial values of the iteration.  Default initial
          value is zero for all variables. }
\item{bwin}{controls the bandwidth used for smoothing in the EM algorithm. The follow-up time is divided into quartiles and
 \code{bwin} specifies a factor by which the maximum between events time length on each interval is multiplied. The default
 \code{bwin=-1} lets the function find an appropriate value. If \code{bwin=0}, no smoothing is applied.
}
\item{centered}{if \code{TRUE}, all the variables are centered before fitting and the baseline excess hazard
	is calculated accordingly. Default is \code{FALSE}. }
\item{cause}{A vector of the same length as the number of cases. \code{0} for population deaths, \code{1} for disease specific
deaths, \code{2} (default) for unknown. Can only be used with the \code{EM} method. }
\item{control}{a list of parameters for controlling the fitting process. 
          See the documentation for \code{glm.control} for details.
}
\item{rmap}{an optional list to be used if the variables are not
          organized and named in the same way as in the \code{ratetable} object. See details below.}

\item{...}{other arguments will be passed to \code{glm.control}.}
}

\details{
NOTE: The follow-up time must be specified in days. The \code{ratetable} being used may have different variable names and formats than the user's data set, this is dealt with by the \code{rmap} argument. For example, if age is in years in the data set but in days in the \code{ratetable} object, age=age*365.241 should be used. The calendar year can be in any date format (date, Date and POSIXt are allowed), the date formats in the \code{ratetable} and in the data may differ. 

The maximum likelihood method and both glm methods assume a fully parametric model with a piecewise constant baseline
excess hazard function. The intervals on which the baseline is assumed constant should be passed via argument \code{int}. The
EM method is semiparametric, i.e. no assumptions are made for the baseline hazard and therefore no intervals need to be specified.

The methods using glm are methods for grouped data. The groups are formed according to the covariate values.
This should be taken into account when fitting a model. The  glm method returns life tables for groups specified by the covariates in \code{groups}.

The EM method output includes the smoothed baseline excess hazard \code{lambda0}, the cumulative baseline excess hazard 
\code{Lambda0} and  \code{times} at which they are estimated. The individual probabilites of dying due to the excess risk
are returned as \code{Nie}. 
The EM method fitting procedure requires some local smoothing of the baseline excess hazard. The default \code{bwin=-1}
 value lets the function find an appropriate value for the 
smoothing band width. While this ensures an unbiased estimate, the procedure time is much longer. As the value found by 
the function is independent of the covariates in the model, the value can be read from the output (\code{bwinfac}) and 
used for refitting different models to the same data to save time.
}

\value{
An object of class \code{rsadd}. In the case of \code{method="glm.bin"} and \code{method="glm.poi"} the class also
inherits from \code{glm} which  inherits from the class \code{lm}.
Objects of this class have methods for the functions \code{print} and \code{summary}.
An object of class \code{rsadd} is a list containing at least the following components:
\item{data}{the data as used in the model, along with the variables defined in the rate table}
\item{ratetable}{the ratetable used.}
\item{int}{the maximum time (in years) used. All the events at and after this value are censored.}
\item{method}{the fitting method that was used.}
\item{linear.predictors}{the vector of linear predictors, one per subject.}
}

\examples{
data(slopop)
data(rdata)
#fit an additive model
#note that the variable year is given in days since 01.01.1960 and that
#age must be multiplied by 365.241 in order to be expressed in days.
fit <- rsadd(Surv(time,cens)~sex+as.factor(agegr)+ratetable(age=age*365.241),
	    ratetable=slopop,data=rdata,int=5)

#check the goodness of fit
rs.br(fit)

#use the EM method and plot the smoothed baseline excess hazard
fit <- rsadd(Surv(time,cens)~sex+age,rmap=list(age=age*365.241),
	 ratetable=slopop,data=rdata,int=5,method="EM")
sm <- epa(fit)
plot(sm$times,sm$lambda,type="l")
}


\references{
Package. Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.


EM algorithm: Pohar Perme M., Henderson R., Stare, J. (2009) "An approach to estimation in relative survival regression." Biostatistics, \bold{10}: 136--146.

}



\seealso{
\code{\link{rstrans}},
\code{\link{rsmul}}
}

\keyword{survival}
