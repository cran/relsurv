\name{survfit.rsadd}
\alias{survfit.rsadd}



\title{Compute a Predicited Survival Curve}

\description{
Computes a predicted survival curve based on the additive model estimated by rsadd function.
}
\usage{
\method{survfit}{rsadd}(formula, newdata, se.fit = TRUE, conf.int = 0.95, individual = FALSE, 
    conf.type = c("log", "log-log", "plain", "none"),...)
}
\arguments{
\item{formula}{a rsadd object}
\item{newdata}{a data frame with the same variable names as those that appear in the rsadd formula. The curve(s) produced will be representative of a cohort who's covariates correspond to the values in newdata.}
\item{se.fit}{a logical value indicating whether standard errors should be computed. Default is \code{TRUE}.}
\item{conf.int}{the level for a two-sided confidence interval on the survival curve(s). Default is 0.95.}
\item{individual}{a logical value indicating whether the data frame represents different time epochs for only one individual (T), or whether multiple rows indicate multiple individuals (F, the default). If the former only one curve will be produced; if the latter there will be one curve per row in newdata.}
\item{conf.type}{One of \code{none}, \code{plain}, \code{log} (the default), or \code{log-log}. The first option causes confidence intervals not to be generated. The second causes the standard intervals curve +- k *se(curve), where k is determined from conf.int. The log option calculates intervals based on the cumulative hazard or log(survival). The last option bases intervals on the log hazard or log(-log(survival)).  }
\item{...}{Currently not implemented}
}

\details{
When predicting the survival curve, the ratetable values for future years will be equal to those of the last given year. The same ratetables will be used for fitting and predicting. To predict a relative survival curve, use \code{rs.surv.rsadd}.
}

\value{
 a \code{survfit} object; see the help on \code{survfit.object} for details.
 The \code{survfit} methods are used for \code{print},
     \code{plot}, \code{lines}, and \code{points}.

}

\references{
Package: Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine,\bold{81}: 272--278.  

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.  
}




\examples{
data(slopop)
data(rdata)
#BTW: work on a smaller dataset here to run the example faster
fit <- rsadd(Surv(time,cens)~sex,rmap=list(age=age*365.241),
	ratetable=slopop,data=rdata[1:500,],method="EM")
survfit.rsadd(fit,newdata=data.frame(sex=1,age=60,year=17000))

}

\seealso{
\code{survfit},
\code{survexp}, \code{\link{rs.surv}}
}

\keyword{survival}
