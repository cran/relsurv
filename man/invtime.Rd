\name{invtime}
\alias{invtime}

\title{Inverse transforming of time in Relative Survival }
\description{
This function can be used when predicting in Relative Survival using the transformed time regression
model (using \code{rstrans} function). It inverses the time from Y to T in relative survival using the
given ratetable. The times Y can be produced with the \code{rstrans} function, in which case, this is the reverse
function. This function does the transformation for one person at a time. 
}
\details{
Works only with ratetables that are split by age, sex and year.
Transforming can be computationally intensive, use lower and/or upper to guess the interval of the result
and thus speed up the function.
}
\usage{
invtime(y, age, sex, year, scale, ratetable, lower, upper)
}
\arguments{
\item{y}{time in Y.}
\item{age}{age of the individual.  Must be in days.}
\item{sex}{sex of the individual. Must be coded in the same way as in the \code{ratetable}.}
\item{year}{date of diagnosis. Must be in a date format}
\item{scale}{
numeric value to scale the results. If \code{ratetable} is in units/day, \code{scale = 365.241} causes the output
to be reported in years.
}
\item{ratetable}{a table of event rates, such as \code{survexp.us}. }
\item{lower}{the lower bound of interval where the result is expected. This argument is optional,
but, if given, can shorten the time the function needs to calculate the result. }
\item{upper}{
the upper bound of interval where the result is expected. See \code{lower} 
}
}
\value{
  A list of values
\item{T}{the original time}
\item{Y}{the transformed time}  
}

\references{
Package: Pohar M., Stare J. (2006) "Relative survival analysis in R."  Computer Methods and Programs in Biomedicine, \bold{81}: 272-278.  

Relative survival: Pohar, M., Stare, J. (2007) "Making relative survival analysis relatively easy."
Computers in biology and medicine, \bold{37}: 1741-1749.  
}



\examples{
data(slopop)
invtime(y = 0.1, age = 23011, sex = 1, year = 9497, ratetable = slopop)
}




\seealso{ \code{\link{rstrans}}  }
\keyword{ survival }
