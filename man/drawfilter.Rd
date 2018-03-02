\name{drawfilter}

\alias{drawfilter}

\title{Draw the amplitude profile of a frequency filter}

\description{This function lets the user drawing the amplitude profile of a frequency filter.}

\usage{drawfilter(f, n = 256, continuous = TRUE, discrete = TRUE)}

\arguments{
  \item{f}{a numeric vector of length 1 for the sampling frequency of
  the object to be filtered (in Hz).}
  \item{n}{a numeric vector of length 1 for the length (i.e. number of points) of the filter. By default = 256 to fit with a FIR with wl = 512.}
  \item{continuous}{a logical (\code{TRUE} by default) to draw a continuous filter.}
  \item{discrete}{a logical (\code{TRUE} by default) to draw a discrete filter.}
}

\details{If the same frequency of a discrete filter is selected twice then the sum of the amplitudes of the two selections is used.
  If both arguments \code{continuous} and \code{discrete} are set to \code{TRUE} and if frequencies selected overlap between the two filters then only the frequencies of the discrete filter are considered.}

\value{The function returns a two-column matrix, the first column is the frequency in kHz and the second column is the amplitude of the filter.}

\author{Laurent Lellouch}

\note{This function can be used to prepare bandpass or bandstop custom filters to be used with \code{\link{fir}} and \code{\link{ffilter}}. See examples.}

\seealso{\code{\link{fir}}, \code{\link{squarefilter}}, \code{\link{combfilter}}, \code{\link{ffilter}}, \code{\link{drawenv}}}

\examples{
\dontrun{
f <- 8000
a <- noisew(f=f, d=1)
## bandpass continuous and discrete
cont.disc <- drawfilter(f=f/2)
a.cont.disc <- fir(a, f=f, custom=cont.disc)
spectro(a.cont.disc, f=f)
## bandpass continuous only
cont <- drawfilter(f=f/2, discrete=FALSE)
a.cont <- fir(a, f=f, custom=cont)
spectro(a.cont, f=f)
## bandstop continuous only
cont.stop <- drawfilter(f=f/2, discrete=FALSE)
a.cont.stop <- fir(a, f=f, custom=cont.stop, bandpass=FALSE)
spectro(a.cont.stop, f=f)
## bandpass discrete only
disc <- drawfilter(f=f/2, continuous=FALSE)
a.disc <- fir(a, f=f, custom=disc, bandpass=FALSE)
spectro(a.disc, f=f)
}
}

\keyword{ts}
\keyword{filter}
