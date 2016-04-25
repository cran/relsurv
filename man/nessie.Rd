\name{nessie}
\alias{nessie}



\title{Net Expected Sample Size Is Estimated}

\description{
     Calculates how the sample size decreases in time due to population mortality
}
\usage{
     nessie(formula, data, ratetable = relsurv::slopop,times)  
}
\arguments{
\item{formula}{
a formula object, same as in \code{rs.surv}. The right-hand side of the formula object includes the variable that defines the subgroups (a variable of type \code{factor}) by which the expected sample size is to be calculated and, if required, the \code{ratetable} argument as in all the other functions of the package.
}
\item{data}{
a data.frame in which to interpret the variables named in the
          \code{formula}. 

}
\item{ratetable}{
a table of event rates, organized as a \code{ratetable} object,
          such as \code{slopop}.  
}
\item{times}{Times at which the calculation should be evaluated - in years!} 
}

\details{
The function calculates the sample size we can expect at a certain time point if the patients die only due to population causes (population survival * initial sample size in a certain category), i.e. the number of individuals that remains at risk at given timepoints after the individuals who die due to population causes are removed.
The result should be used as a guideline for the sensible length of follow-up interval when calculating the net survival.

The first column of the output reports the number of individuals at time 0.
The last column of the output reports the conditional expected (population) survival time for each subgroup.
}

\value{
A list of values.
}




\examples{
data(slopop)
data(rdata)
rdata$agegr <-cut(rdata$age,seq(40,95,by=5))
nessie(Surv(time,cens)~agegr+ratetable(age=age*365.241,sex=sex,year=year),
	ratetable=slopop,data=rdata,times=c(1,3,5,10,15))
}


\seealso{
\code{rs.surv}
}

\keyword{survival}
