\name{smoothw}

\alias{smoothw}

\title{A function to tentavily smooth a time wave}

\description{This function tries to smooth with a sum sliding window a time wave, and then to remove residual noise.}

\usage{smoothw(wave, f, wl, padding=TRUE, output = "matrix")}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{window length in number of points (samples).}
  \item{padding}{a logical, if \code{TRUE} add 0 values at the start and
  end of the file to match \code{wave} length (duration).}
  \item{output}{character string, the class of the object to return,
    either \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
}

\details{A window slides along the signal and sums up the sample
  amplitude values. Zero values are added at the end of the wave to keep
  wave length (duration).}

\value{A new wave is returned. The class of the returned object is set
  with the argument \code{output}. If \code{padding} is \code{TRUE}, the
  new wave starts and ends up with \code{0} values to match the size of \code{wave}.}

\author{Jerome Sueur}

\section{Warning }{This function should be used with care as this kind of filter may change the frequency content of the sound. See the examples section for an illustration.}

\seealso{\code{\link{fir}}, \code{\link{filter}}}

\examples{
# An example to show that smoothw() may change
# the frequency content of your sound
data(orni)
orni2 <- smoothw(orni, wl=2, out="Wave")
orni10 <- smoothw(orni, wl=10, out="Wave")
orni50 <- smoothw(orni, wl=50, out="Wave")
orni100 <- smoothw(orni, wl=100, out="Wave")
meanspec(orni)
lines(meanspec(orni2, plot=FALSE), col=2)
lines(meanspec(orni10, plot=FALSE), col=3)
lines(meanspec(orni50, plot=FALSE), col=4)
lines(meanspec(orni100, plot=FALSE), col=5)
legend("topright", col=1:5, lty=1, legend=c("original","wl=2","wl=10","wl=50","wl=100"))
}

\keyword{ts}
