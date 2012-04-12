\name{ceps}

\alias{ceps}

\title{Cepstrum or real cepstrum}

\description{This function returns the cepstrum of a time wave allowing fundamental
  frequency detection.}

\usage{
ceps(wave, f, wl = 512, at = NULL, from = NULL, to = NULL,
tidentify = FALSE, fidentify = FALSE, col = "black", cex = 1, plot = TRUE,
qlab = "Quefrency (bottom: s, up: Hz)", alab = "Amplitude",
qlim = NULL, alim = NULL, type = "l", ...)
}

\arguments{
  \item{wave}{an R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{at}{position where to compute the cepstrum (in s).}
  \item{from}{start position where to compute the cepstrum (in s).}
  \item{to}{end position to compute the cepstrum (in s).}
  \item{tidentify}{to identify time values on the plot with the help of a cursor.}
  \item{fidentify}{to identify frequency values on the plot with the help of a cursor.}
  \item{col}{colour of the cepstrum.}
  \item{cex}{pitch size of the cepstrum.}
  \item{plot}{logical, if \code{TRUE} plots the cepstrum.}
  \item{qlab}{title of the quefrency axis (in s).}
  \item{alab}{title of the amplitude axis.}
  \item{qlim}{range of quefrency axis.}
  \item{alim}{range of amplitude axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}  
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  The cepstrum of a time wave is the inverse Fourier transform of the logarithm
  of the Fourier transform. The cepstrum of a wave \emph{s} is then calculated as follows:\cr
  \deqn{C(s) = Re[FFT^{-1}(\log{(|FFT(s)|)]}}{% 
    C(s) = Re(invfft(log(abs(fft(s)))))}\cr

  The independent variable of a cepstral graph is called the quefrency.
  The quefrency is a measure of time, though not in the sense of a signal
  in the time domain. A correspondence with the frequency domain is obtained
  by simply computing the reverse of the temporal x coordinate. For instance if
  a peak appears at 0.005 s, this reveals a frequency peak at 200 Hz (=1/0.005).
  This explain the two scales plotted when \code{plot} is \code{TRUE}. \cr
  If \code{at}, \code{from} or \code{to} are \code{FALSE} then \code{ceps} 
  computes the cepstrum of the whole signal.\cr
  When using \code{tidentify} or \code{tidentify}, press \sQuote{stop}
  tools bar button to return values in the console.
}

\value{
  When \code{plot} is \code{FALSE}, \code{ceps} returns the cesptral profile as a two-column matrix, the first column corresponding to quefrency (\emph{x}-axis) and the second
  corresponding to amplitude (\emph{y}-axis).
}

\references{Oppenheim, A.V. and Schafer, R.W. 2004. From frequency
  to quefrency: a history of the cepstrum. \emph{Signal Processing Magazine IEEE},
  21: 95-106.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{Cepstral analysis is mainly used in speech processing.
  This analysis allows to extract the fundamental frequency, see the examples.\cr
  This function is based on \code{\link{fft}}.}

\section{Warning}{The argument \code{peaks} is no more available
  (version > 1.5.6). See the function \code{\link{fpeaks}}
  for peak(s) detection.} 

\seealso{\code{\link{cepstro}}, \code{\link{fund}}, \code{\link{autoc}}}

\examples{
data(sheep)
ceps(sheep,f=8000,at=0.4,wl=1024)
}

\keyword{dplot}
\keyword{ts}
