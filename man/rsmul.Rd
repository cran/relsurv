\name{rsmul}
\alias{rsmul}
\title{Fit Andersen et al Multiplicative Regression Model for Relative Survival}
\description{
  Fits the Andersen et al multiplicative regression model in relative survival. An extension of the 
coxph function using relative survival.
}
\usage{
rsmul(formula, data, ratetable = relsurv::slopop, int,na.action,init,
      method,control,rmap,...)
}
\arguments{
\item{formula}{
a formula object, with the response as a \code{Surv} object on the
          left of a \code{~} operator, and, if desired, terms separated  by
          the \code{+} operator on the right. 

 NOTE: The follow-up time must be in days.
}
\item{data}{
a data.frame in which to interpret the variables named in
the \code{formula}.
}
\item{ratetable}{
a table of event rates, such as \code{slopop}. 
}
\item{int}{
the number of follow-up years used for calculating survival(the data are censored after this time-point). If missing, it is set the 
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
\item{rmap}{an optional list to be used if the variables are not
          organized and named in the same way as in the \code{ratetable} object. See details below.}

\item{...}{Other arguments will be passed to \code{coxph.control}.}
}



\value{
an object of class \code{coxph} with an additional item:
\item{basehaz}{Cumulative baseline hazard (population values are seen as offset) at centered values of covariates.}
}

\details{
NOTE: The follow-up time must be specified in days. The \code{ratetable} being used may have different variable names and formats than the user's data set, this is dealt with by the \code{rmap} argument. For example, if age is in years in the data set but in days in the \code{ratetable} object, age=age*365.241 should be used. The calendar year can be in any date format (date, Date and POSIXt are allowed), the date formats in the \code{ratetable} and in the data may differ. 
}

\references{
Method: Andersen, P.K., Borch-Johnsen, K., Deckert, T., Green, A.,
Hougaard, P., Keiding, N. and Kreiner, S. (1985) "A Cox regression
model for relative mortality and its application to diabetes
mellitus survival data.", Biometrics, \bold{41}: 921--932.



Package. Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.
}

\examples{
data(slopop)
data(rdata)
#fit a multiplicative model
#note that the variable year is given in days since 01.01.1960 and that 
#age must be multiplied by 365.241 in order to be expressed in days.
fit <- rsmul(Surv(time,cens)~sex+as.factor(agegr),rmap=list(age=age*365.241),
            ratetable=slopop,data=rdata)


#check the goodness of fit
rs.br(fit)

}




\seealso{\code{\link{rsadd}}, \code{\link{rstrans}}. }

\keyword{survival}


