\name{plot.rs.zph}
\alias{plot.rs.zph}


\title{Graphical Inspection of Proportional Hazards Assumption in Relative Survival Models}
\description{
Displays a graph of the scaled partial residuals, along with a smooth curve.
}
\usage{
\method{plot}{rs.zph}(x, resid=TRUE, df = 4, nsmo = 40,var,cex=1,add=FALSE,col=1,
    lty=1,xlab,ylab,scale=1,...)
}

\arguments{
\item{x}{
result of the \code{rs.zph} function.
}
\item{resid}{a logical value, if \code{TRUE} the residuals are included on the
          plot, as well as the smooth fit.}
\item{df}{
the degrees of freedom for the fitted natural spline, \code{df=2}
          leads to a linear fit. 
}
\item{nsmo}{
number of points used to plot the fitted spline. 
}
\item{var}{ 
the set of variables for which plots are desired.  By
          default, plots are produced in turn for each variable of a
          model.  Selection of a single variable allows other features
          to be added to the plot, e.g., a horizontal line at zero or a
          main title.
}
\item{cex}{a numerical value giving the amount by which plotting text
          and symbols should be scaled relative to the default.
}
\item{add}{logical, if \code{TRUE} the plot is added to an existing plot}
\item{col}{a specification for the default plotting color.}
\item{lty}{the line type.}
\item{xlab}{x axis label.}
\item{ylab}{y axis label.}
\item{scale}{units for x axis, default is 1, i.e. days.}
\item{...}{Additional arguments passed to the \code{plot} function. }
}


\references{
\item{Goodness of fit}{
 Stare J.,Pohar Perme M., Henderson R. "Goodness of fit of relative survival models." \emph{Statistics in Medicine}, 24: 3911-3925, 2005. 
}

\item{package}{Pohar M., Stare J. "Relative survival analysis in R." \emph{Computer Methods and Programs in Biomedicine,} 81: 272-278, 2006.}

\item{relative survival}{Pohar, M., Stare, J. "Making relative survival analysis relatively easy."
\emph{Computers in biology and medicine}, 37: 1741-1749, 2007.}
}




\examples{
data(slopop)
data(rdata)
fit <- rsadd(Surv(time,cens)~sex+as.factor(agegr)+ratetable(age=age*365,
             sex=sex,year=year),ratetable=slopop,data=rdata,int=5)
rszph <- rs.zph(fit)
plot(rszph)
}

\seealso{\code{\link{rs.zph}}, \code{\link{plot.cox.zph}}.}

\keyword{survival}



