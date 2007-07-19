\name{lfs}

\alias{lfs}

\title{Linear Frequency Shift}

\description{
This function linearly shifts all the frequency content of a time wave. 
}

\usage{lfs(wave, f, shift, wl = 512, wn = "hanning")}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{shift}{positive or negative frequency shift to apply (in Hz).)}
  \item{wl}{window length for the analysis (even number of points, by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
}

\details{
A short-term Fourier transform is first applied to the signal
(see \code{\link{spectro}}), then the frequency shift is applied and the new
signal is eventually generated using the reverse of the Fourier Transform
(\code{\link{fft}}).\cr
There is therefore neither temporal modifications nor 
amplitude modifications.
}

\value{
\code{fls} returns a one-column matrix describing the new wave.
}

\references{
Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr} and 
Thierry Aubin \email{thierry.aubin@u-psud.fr}}

\seealso{\code{\link{ffilter}}, \code{\link{spectro}}}

\examples{
data(orni)
a<-lfs(orni,f=22050,shift=1000)
spectro(a,f=22050)
# to be compared with the original signal
spectro(orni,f=22050)
}

\keyword{ts}
