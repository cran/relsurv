\name{joinrate}
\alias{joinrate}
\title{Join ratetables}
\description{
 The function joins two or more objects organized as \code{ratetable} by adding a new dimension. 
}
\usage{
joinrate(tables,dim.name="country")
}

\arguments{
  \item{tables}{
	a list of ratetables. If names are given, they are included as \code{dimnames}. 
  }
  \item{dim.name}{
	the name of the added dimension.   }
}
\details{ 
This function joins two or more \code{ratetable} objects by adding a new dimension. The cutpoints of all the
rate tables are compared and only the common intervals kept. If the intervals defined by the cutpoints are not of
the same length, a warning message is displayed.  Each rate table must have 3 dimensions, i.e. age, sex and year 
(the order is not important).

}
\value{An object of class \code{ratetable}.}

\references{
Package: Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272-278.

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy." Computers in biology and medicine, \bold{37}: 1741-1749.  
}



\examples{
#newpop <- joinrate(list(Arizona=survexp.az,Florida=survexp.fl,
#                   Minnesota=survexp.mn),dim.name="state")
}
\seealso{\code{\link{ratetable}}, \code{\link{transrate.hld}}, \code{\link{transrate.hmd}}, \code{\link{transrate}}.}
\keyword{survival}
					   