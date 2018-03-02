\name{rmam}

\alias{rmam}

\title{Remove the amplitude modulations of a time wave}

\description{This functions removes the amplitude modulation of a time wave
  through the Hilbert amplitude envelope.}

\usage{rmam(wave, f, plot = FALSE, listen = FALSE, output = "matrix", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot
    of the nwe time wave (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{The new time wave is obtained by dividing the original time wave by
  its Hilbert amplitude envelope.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\references{Mbu Nyamsi, R. G., Aubin, T. & Bremond, J. C. 1994
  On the extraction of some time dependent parameters of an acoustic signal
  by means of the analytic signal concept.
  Its application to animal sound study. \emph{Bioacoustics}, 5: 187-203.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{hilbert}}.}

\examples{
# generate a new sound with amplitude modulation
a<-synth(f=8000, d=1, cf=1500, am=c(50,10))
# remove the amplitude modulation and plot the result
rmam(a,f=8000,plot=TRUE)
}

\keyword{ts}