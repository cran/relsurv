\name{srvxp.fit}
\alias{srvxp.fit}

\title{Work Function to Compute Expected Survival}
\description{
  Compute expected survival for a group of people. This function is not to be called by
      the user.
}
\usage{
srvxp.fit(x, y, ratetable)
}

\arguments{
  \item{x}{a  data.frame or a matrix. Each row represents one person. The   columns must match the
  	dimensions of the \code{ratetable}, in the correct order. 
 }
  \item{y}{a vector of survival times of the same length as the number of rows of x.}
 \item{ratetable}{
 a table of event rates, organized as a \code{ratetable} object, such as \code{survexp.us}. 
}
}
\value{A vector of survival times.}

\keyword{survival}
