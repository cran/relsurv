\name{transrate}
\alias{transrate}
\title{Reorganize Data into a Ratetable Object}
\description{
 The function assists in reorganizing certain types of data
 into a ratetable object.
}
\usage{
transrate(men,women,yearlim,int.length=1)
}

\arguments{
  \item{men}{
a matrix containing the yearly (conditional) probabilities
of one year survival for men. Rows represent age (increasing
1 year per line,starting with 0), the columns represent cohort
years (the limits are in \code{yearlim}, the increase is in
\code{int.length}.
  }
  \item{women}{
a matrix containing the yearly (conditional) probabilities
of one year survival for women.
  }
  \item{yearlim}{the first and last cohort year given in the tables.}
  \item{int.length}{the length of intervals in which cohort years are given.}
}
\details{ 
This function only applies for ratetables that are
organized by age, sex and year.
}
\value{An object of class \code{ratetable}.}
\examples{
men <- cbind(exp(-365.241*exp(-14.5+.08*(0:100))),exp(-365*exp(-14.7+.085*(0:100))))
women <- cbind(exp(-365.241*exp(-15.5+.085*(0:100))),exp(-365*exp(-15.7+.09*(0:100))))
table <- transrate(men,women,yearlim=c(1980,1990),int.length=10)
}

\references{
Package. Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.

}


\seealso{\code{\link{ratetable}}.}
\keyword{survival}
