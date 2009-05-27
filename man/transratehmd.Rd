\name{transrate.hmd}
\alias{transrate.hmd}
\title{Reorganize Data obtained from Human Mortality Database into a Ratetable Object}
\description{
 The function assists in reorganizing the .txt files obtained from Human Mortality Database 
(http://www.mortality.org)  into a ratetable object.
}
\usage{
transrate.hmd(male,female)
}

\arguments{
  \item{male}{
	a .txt file, containing the data on males. 
  }
  \item{female}{
	a .txt file, containing the data on females.   }
}
\details{ 
This function works with tables organised in the format provided by the Human Mortality 
Database. The period lifetables, organized by date of death, must be downloaded separately for men and
women for each country and the top line deleted. The function only works for tables that have yearly
cutpoints for age as well as calendar year (1x1). The only two columns needed are calendar year (Year) 
and probability of death (qx). Death probabilities must be calculated up to age 110 (in yearly intervals).
}
\value{An object of class \code{ratetable}.}

\references{
Package. Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272--278

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741--1749.
}



\examples{
\dontrun{
auspop <- transrate.hmd("mltper_1x1.txt","fltper_1x1.txt")
}
}
\seealso{\code{\link{ratetable}}, \code{\link{transrate.hld}}, \code{\link{joinrate}}, \code{\link{transrate}}.}
\keyword{survival}
