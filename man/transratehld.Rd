\name{transrate.hld}
\alias{transrate.hld}
\title{Reorganize Data obtained from Human Life-Table Database into a Ratetable Object}
\description{
 The function assists in reorganizing the .txt files obtained from Human Life-Table Database 
(http://www.lifetable.de -> Data by Country)  into a ratetable object.
}
\usage{
transrate.hld(file,cut.year,race)
}

\arguments{
  \item{file}{
	a vector of file names which the data are to be read from. Must be in
	.tex format and in the same format as the files in Human Life-Table Database.
  }
  \item{cut.year}{
	a vector of cutpoints for years. Must be specified when	the year spans in the
	files are not consecutive.
  }
  \item{race}{a vector of race names for the input files.}
}
\details{ 
This function works with any table organised in the format provided by the Human Life-Table
Database, but currently only works with TypeLT 1 (i.e. age intervals of length 1). The age must always 
start with value 0, but can end at different values (when that happens, the last value is carried forward).
 The rates between the cutpoints are taken to be constant.
}
\value{An object of class \code{ratetable}.}
\examples{
\dontrun{
finpop <- transrate.hld(c("FIN_1981-85.txt","FIN_1986-90.txt","FIN_1991-95.txt"))
}
\dontrun{
nzpop <- transrate.hld(c("NZL_1980-82_Non-maori.txt","NZL_1985-87_Non-maori.txt",
	 "NZL_1980-82_Maori.txt","NZL_1985-87_Maori.txt"),
	 cut.year=c(1980,1985),race=rep(c("nonmaori","maori"),each=2))
}
}

\references{
Pohar M., Stare J. "Relative survival analysis in R." \emph{Computer Methods and Programs in Biomedicine,} 81: 272-278, 2006.
}

\seealso{\code{\link{ratetable}}, \code{\link{transrate.hmd}}, \code{\link{joinrate}}, \code{\link{transrate}}.}
\keyword{survival}
