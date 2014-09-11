\name{istft}

\alias{istft}

\title{Inverse of the short-term Fourier transform}

\description{This function returns a wave object from a complex STFT matrix by computing the inverse of the short-term
  Fourier transform (STFT)}

\usage{istft(stft, f, wl, ovlp=75, wn="hanning", output = "matrix")}

\arguments{
  \item{stft}{a complex matrix resulting of a short-term Fourier transform.}
  \item{f}{sampling frequency of the original \code{wave} object (in Hz)}
  \item{wl}{FFT window length for the analysis (even number of points).}
  \item{ovlp}{overlap between successive FFT windows (in \%, by default
  75\%, see the Details section).}
  \item{wn}{character string specifying the FFT window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"audioSample"} or \code{"ts"}.}
}

\details{The function is based on the inverse of the FFT (see \code{\link{fft}}) and on
  the overlap add (OLA) method.\cr
  The overlap percentage must satisfy the Perfect Reconstruction OLA-constraint. For
  the most windows, this constraint is: \cr
  \deqn{ovlp = 100 \times (1 - \frac{1}{4 \times n}),}{%
  ovlp = 100 * (1 - 1/(4 * n)), }\cr
with \emph{n} being a positive integer.\cr
A default value is set to 75\%. We suggest not to change it.
}

\value{A new wave is returned. The class of the returned object is set with the argument \code{output}.}

\author{Original Matlab code by Hristo Zhivomirov (Technical University of Varna, Bulgaria),
  translated and adapted to R by Jerome Sueur}

\note{The \code{stft} input data must be complex. \cr
  This function is used by \code{\link{ffilter}}, \code{\link{lfs}} to
  respectively filter in frequency and shift in frequency a sound. \cr
  The function can be used to reconstruct or modify a sound. See examples.}

\seealso{\code{\link{spectro}}, \code{\link{ffilter}}, \code{\link{lfs}}}

\examples{
\dontrun{
# STFT and iSTFT parameters
wl <- 1024
ovlp <- 75
# reconstruction of the tico sound from the stft complex data matrix
data(tico)
data <- spectro(tico, wl=wl, ovlp=ovlp, plot=FALSE, norm=FALSE, dB=NULL, complex=TRUE)$amp
res <- istft(data, ovlp=ovlp, wn="hanning", wl=wl, f=22050, out="Wave")
spectro(res)
# a strange frequency filter
n <- noisew(d=1, f=44100)
data <- spectro(n, f=44100, wl=wl, ovlp=ovlp, plot=FALSE, norm=FALSE, dB=NULL, complex=TRUE)$amp
data[64:192, 6:24] <- 0 
nfilt <- istft(data, f=8000, wl=wl, ovlp=ovlp, output="Wave")
spectro(nfilt, wl=wl, ovlp=ovlp)
}
}

\keyword{ts}
