\name{rmam}

\alias{rmam}

\title{Remove the amplitude modulations of a time wave}

\description{This functions removes the amplitude modulation of a time wave
through the Hilbert amplitude envelope.}

\usage{rmam(wave, f, plot = FALSE, listen = FALSE, Sample = FALSE, ...)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot
  of the nwe time wave (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.  
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{The new time wave is obtained by dividing the original time wave by
its Hilbert amplitude envelope.}

\value{If \code{plot} is \code{FALSE}, a new time wave is returned as a
one-column matrix or as a \code{\link[sound]{Sample}} object
if \code{Sample} is \code{TRUE}.}

\references{Mbu Nyamsi, R. G., Aubin, T. & Br�mond, J. C. 1994
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