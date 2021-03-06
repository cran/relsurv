\name{survsplit}
\alias{survsplit}

\title{Split a Survival Data Set at Specified Times}
\description{
  Given a survival data set and a set of specified cut times, the function splits
 each record into multiple records at each cut time. The new data 
 set is be in \code{counting process} format, with a start time, stop
time, and event status for each record. More general than \code{survSplit}
 as it also works with the data already in the \code{counting process} format.
}
\usage{
survsplit(data, cut, end, event, start, id = NULL, zero = 0,
          episode = NULL,interval=NULL)
}

\arguments{
  \item{data}{data frame. }
  \item{cut}{vector of timepoints to cut at.}
  \item{end}{character string with name of event time variable. }
  \item{event}{character string with name of censoring indicator. }
  \item{start}{character string with name of start variable (will be
	created if it does not exist). }
  \item{id}{character string with name of new id variable to create
(optional). }
  \item{zero}{If \code{start} doesn't already exist, this is the time
that the original records start. May be a vector or single value. }
  \item{episode}{character string with name of new episode variable (optional).}
  \item{interval}{this argument is used by \code{max.lik} function}
}
\value{New, longer, data frame.}

\seealso{\code{\link{survSplit}}.}
\keyword{survival}
