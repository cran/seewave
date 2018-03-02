\name{hilbert}

\alias{hilbert}

\title{Hilbert transform and analytic signal}

\description{This function returns the analytic signal of a time wave
  through Hilbert transform.}

\usage{hilbert(wave, f, fftw = FALSE)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
    library \code{fftw} for faster computation. See Notes of the function \code{spectro}.}
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

\author{Jonathan Lees \email{jonathan.lees@unc.edu}. Implementation of
  'fftw' argument by Jean Marchal and Francois Fabianek.}

\note{To get the Hilbert component only, use \code{Im(Hilbert(wave))}.}

\seealso{\code{\link{ifreq}}}

\examples{
a<-synth(f=8000, d=1, cf=1000)
aa<-hilbert(a, f=8000)
}

\keyword{ts}
