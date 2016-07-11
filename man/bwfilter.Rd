\name{bwfilter}

\alias{bwfilter}

\title{Butterworth frequency filter}

\description{This functin is a Butterworth frequency filter that filters
  out a selected frequency section of of a time wave (low-pass, high-pass, low-stop, high-stop, bandpass or bandstop frequency filter).}

\usage{
bwfilter(wave, f, n = 1, from = NULL, to = NULL,
bandpass = TRUE, listen = FALSE, output = "matrix")
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{n}{Order of the filter. See details.}
  \item{from}{start frequency (in Hz) where to apply the filter.}
  \item{to}{end frequency (in Hz) where to apply the filter.}
  \item{bandpass}{if \code{TRUE} a band-pass filter is applied between
  \code{from} and  \code{to}, if not \code{NULL} a band-stop filter is applied
    between \code{from} and  \code{to} (by default \code{NULL}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
}

\details{
The order of the filter determines the value of the roll-off value, that
is the dB decrease per octave of the transfer function. A filter of order \emph{n} will have a
transfer function with a roll-off value of - \emph{n*6} dB.
}

\value{A new wave is returned. The class
  of the returned object is set with the argument \code{output}.
}

\references{Stoddard, P. K. (1998). Application of filters in bioacoustics.
  \emph{In}: Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds), \emph{Animal acoustic 
    communication}. Springer, Berlin, Heidelberg,pp. 105-127.
}

\author{Jerome Sueur, functions \code{filter()} and
  \code{filtfilt()} from the package \code{signal}.
}

\note{This function mainly uses the functions \code{filter()} and
\code{filtfilt()} from the package \code{signal}
}

\seealso{\code{\link{ffilter}}, \code{\link{bwfilter}}, \code{\link{preemphasis}}, \code{\link{lfs}}, \code{\link{afilter}}}

\examples{
require(signal)
f <- 8000
a <- noisew(f=f, d=1)
## low-pass
# 1st order filter
res <- bwfilter(a, f=f, n=1, to=1500)
spectro(res, f=f)
res <- bwfilter(a, f=f, n=2, to=1500)
# 2nd order filter
spectro(res,f=f)
# 8th order filter
res <- bwfilter(a, f=f, n=8, to=1500)
spectro(res,f=f)
## high-pass
res <- bwfilter(a, f=f, from=2500)
spectro(res,f=f)
## band-pass
res <- bwfilter(a, f=f, from=1000, to=2000)
spectro(res,f=f)
## band-stop
res <- bwfilter(a, f=f, from=1000, to=2000,bandpass=FALSE)
spectro(res,f=f)
}

\keyword{ts}
\keyword{filter}
