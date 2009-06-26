\name{hilbert}

\alias{hilbert}

\title{Hilbert transform and analytic signal}

\description{This function returns the analytic signal of a time wave
through Hilbert transform.}

\usage{hilbert(wave, f)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
}

\details{The analytic signal is useful to get the amplitude envelope 
(see argument \code{henv} of \code{\link{oscillo}} and the instantaneous
phase or frequency (see \code{\link{ifreq}}) of a time wave.}

\value{\code{hilbert} returns the analytic signal as a complex matrix.
The imaginary part of this matrix is the Hilbert transform.}

\references{Mbu Nyamsi, R. G., Aubin, T. & Bremond, J. C. 1994
On the extraction of some time dependent parameters of an acoustic signal
by means of the analytic signal concept.
Its application to animal sound study. \emph{Bioacoustics}, 5: 187-203.}

\author{Jonathan Lees \email{jonathan.lees@unc.edu}}

\note{To get the Hilbert component only, use \code{Im(Hilbert(wave))}.}

\seealso{\code{\link{ifreq}}, argument \code{henv} of \code{\link{oscillo}}}

\examples{
a<-synth(f=8000, d=1, cf=1000)
aa<-hilbert(a, f=8000)
}

\keyword{ts}
