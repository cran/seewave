\name{duration}
\alias{duration}
\title{Duration of a time wave}

\description{Returns the duration (in second) of a time wave}
\usage{duration(wave, f)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
}

\value{A numeric vector of length 1 returning the duration in second.}

\author{Jerome Sueur}

\examples{
data(tico)
duration(tico)
}

\keyword{ts}
