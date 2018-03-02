\name{stft.ext}

\alias{stft.ext}

\title{Short-term Fourier transform using fftw and libsndfile C libraries}

\description{
This function performs the short-term Fourier transform externally,
i. e. outside R, on an audio .wav file}

\usage{
stft.ext(file, wl = 512, ovlp = 0, mean = FALSE,
norm = FALSE, dB = FALSE, verbose = FALSE) 
}

\arguments{
  \item{file}{a character string vector of length 1, corresponding to
  the audio file name. Only \code{.wav} files are supported.}
  \item{wl}{a numeric vector of length 1, window length for the analysis (even number of points) (by default = 512).}
  \item{ovlp}{a numeric vector of length 1, overlap between two successive windows (in \%).}
  \item{mean}{logical, if \code{TRUE} calculates the mean spectrum.}
  \item{norm}{logical, if \code{TRUE} the spectrogram matrix and the
  mean spectrum if \code{mean} is \code{TRUE} is/are  normalized to a maximum of 1.}
  \item{dB}{logical, if \code{TRUE} the values are in deciBel with a maximum of 0.}
  \item{verbose}{logical, if \code{TRUE} displays audio file meta-information.}
}

\value{
  If \code{mean} is \code{FALSE}, the function returns a two-column
  matrix corresponding to the amplitude values of the spectrogram. Each
  column is a discrete Fourier transform of length \code{wl}. \cr

  If \code{mean} is \code{TRUE}, the function returns a list with
  two elements:\cr 
  \item{mean}{a vector of length \code{wl/2} corresponding to the mean spectrum}\cr 
  \item{amp}{the spectrogram matrix.}
}

\author{Alicia Stotz}

\note{This function was developped to speed up the process of computing
  a spectrogram or a mean spectrum. The function should be faster than
  \code{spectro} and \code{meanspec} respectively as it does not import
  the \code{.wav} file into \code{R} and uses \code{C}
  fftw3 (\url{http://www.fftw.org/}) and libsndfile
  (\url{http://www.mega-nerd.com/libsndfile/}) libraries. The function
  is typically made for batch processing or long files. Gain in process
  time is significant for files > 10 minutes and for more thant 1000
  files lasting one second.\cr
  A Hanning window shape is used by default.
} 

\seealso{\code{\link{spectro}}, \code{\link{meanspec}}, \code{\link{fft}}.}

\examples{
# tico data
data(tico)
# write a local .wav file
savewav(tico, file="tico.wav")
# spectrogram - not normalised - linear scale - meta-information returned
res <- stft.ext(file="tico.wav", verbose=TRUE)
# spectrogram - normalised - linear scale - no meta-information
res <- stft.ext(file="tico.wav", norm=TRUE)
# spectrogram - dB scale - no meta-information
res <- stft.ext(file="tico.wav", dB=TRUE)
# see how it looks like (no scale)
filled.contour(t(res))
# spectrogram and mean spectrum - normalised - linear scale
res <- stft.ext(file="tico.wav", norm = TRUE, mean = TRUE)
# remove .wav file
unlink("tico.wav")
}

\keyword{ts}
