\name{ceps}

\alias{ceps}

\title{Cepstrum or real cepstrum}

\description{This functions returns the cepstrum of a time wave allowing fundamental
frequency detection.}

\usage{
ceps(wave, f, wl = 512, at = FALSE, from = FALSE, to = FALSE, tpeaks = FALSE,
fpeaks = FALSE, identify = FALSE, col = "black", cex = 1, colpeaks = "red", cexpeaks = 0.75,
fontpeaks = 1, plot = TRUE, qlab = "Quefrency (s)", alab = "Amplitude",
type = "l", ...)
}

\arguments{
  \item{wave}{data describing a time wave or a \code{\link[Sound]{Sample}}
  object generated loading a wav file with \code{\link[Sound]{loadSample}}
  (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{at}{position where to compute the cepstrum (in s).}
  \item{from}{start mark where  to compute the cepstrum (in s).}
  \item{to}{end mark where to compute the cepstrum (in s).}
  \item{tpeaks}{returns peaks value for a given span
      according to time scale (s)(see \code{details}).}
  \item{fpeaks}{returns peaks value for a given span
      according to frequency scale (Hz)(see \code{details}).}
  \item{identify}{to identify frequency values on the plot with the help of a cursor.}
  \item{col}{colour of the cepstrum.}
  \item{cex}{pitch size of the cepstrum.}
  \item{colpeaks}{colour of peaks value plotted on the cepstrum.}
  \item{cexpeaks}{character size of peaks value plotted on the cepstrum.}
  \item{fontpeaks}{font of peaks value plotted on the cepstrum.}
  \item{plot}{logical, if \code{TRUE} plots the cepstrum.}
  \item{qlab}{title of the quefrency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
      See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
The independent variable of a cepstral graph is called the quefrency.
The quefrency is a measure of time, though not in the sense of a signal
in the time domain. A correspondence with the frequency domain is obtained
by simply computing the reverse of the temporal x corrdinate. For instance if
a peak appears at 0.005 s, this reveals a frequency peak at 200 Hz (=1/0.005).\cr
If \code{at}, \code{from} or \code{to} are \code{FALSE} then \code{ceps} 
computes the cepstrum of the whole signal.\cr
\code{tpeaks} and \code{fpeaks} setting corresponds to \code{dimension}
of \code{\link{embed}}.\cr  When using \code{identify}, press \sQuote{stop}
tools bar button to return values in the console. \code{tpeaks} and \code{fpeaks} 
just differ in the unit of the results. 
}

\value{
The cepstrum of a time wave is the Fourier transform of the logarithm
of the Fourier transform. The cepstrum of a wave \emph{s} is then calcuated as follows:\cr
\deqn{C(s) = Re[FFT^{-1}(\log{|FFT(s)|)]}}{% 
      C(s) = Re(invfft(log(abs(fft(s)))))}
}

\references{Oppenheim, A.V. and Schafer, R.W. 2004. From frequency
to quefrency: a history of the cepstrum. \emph{Signal Processing Magazine IEEE},
21: 95-106.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr} and 
Martin Maechler \email{maechler@stat.math.ethz.ch} for \code{peaks}.}

\note{Cepstral analysis is mainly used in speech processing.
This analysis allows to extract the fundamental frequency, see the examples.\cr
This function is based on \code{\link{fft}}.} 

\seealso{\code{\link{autoc}}}

\examples{
data(sheep)
ceps(sheep,f=22050,at=0.4,wl=1024)
# peaks detection in Hertz, the fundamental is at 416 Hz.
ceps(sheep,f=22050,at=0.4,wl=1024,fpeaks=63)
}


\keyword{dplot}
\keyword{ts}
