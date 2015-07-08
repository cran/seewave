\name{squarefilter}

\alias{squarefilter}

\title{Frequency square filter}

\description{This function prepares the amplitude profile of a square frequency filter.}

\usage{squarefilter(f, from = NULL, to = NULL, bandpass = TRUE, wl = 1024)}

\arguments{
  \item{f}{a numeric vector of length 1 for the sampling frequency of the object to be filtered (in Hz).}
  \item{from}{a numeric vector for the start frequencies (in Hz) where to apply the filter.}
  \item{to}{a numeric vector of the end frequencies (in Hz) where to apply the filter.}
  \item{bandpass}{if \code{TRUE} a band-pass filter is prepared between
  start and end frequencies (arguments \code{from} and \code{to}), if \code{FALSE} a bandstop filter is prepared.}
  \item{wl}{window length of the impulse filter (even number of points).}
}

\value{The function returns a two-column matrix, the first column is the frequency in kHz and the second column is the amplitude of the filter (frequency response of the filter).}

\author{Laurent Lellouch}

\note{This function can be used to prepare bandpass or bandstop filters to be used with \code{\link{fir}} and \code{\link{ffilter}}. See examples.}

\seealso{\code{\link{fir}}, \code{\link{drawfilter}}, \code{\link{ffilter}}, \code{\link{combfilter}}}

\examples{
f <- 44100
a <- noisew(f = f, d = 1)
p <- squarefilter(f, from = c(100, 1000, 4000), to = c(500, 3000, 8000))
plot(p, type="l")
h <- fir(a, f = f, custom = p, wl = 1024, output = 'Wave')
spectro(h)
}

\keyword{ts}
\keyword{filter}
