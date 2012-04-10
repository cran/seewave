\name{csh}

\alias{csh}

\title{Continuous spectral entropy}

\description{This function computes the continuous spectral entropy (H) of a time wave.}

\usage{
csh(wave, f, wl = 512, wn = "hanning", ovlp = 0, fftw = FALSE, threshold = NULL,
plot = TRUE, xlab = "Times (s)", ylab = "Spectral Entropy",
ylim = c(0, 1.1), type = "l", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes of the \code{spectro}.}
  \item{threshold}{amplitude threshold for signal detection (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the spectral entropy against time
    (by default \code{TRUE}).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{ylim}{the range of y values.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{See \code{\link{sh}} for computing method.}

\value{When \code{plot} is \code{FALSE}, \code{csh} returns a two-column matrix, the
  first column being time in seconds (\emph{x}-axis) and the second column being the spectral entropy (\emph{y}-axis) computed along time.\cr
  NA corresponds to pause sections in \code{wave} (see \code{threshold}). 
}

\references{Toh, A. M., Togneri, R. & Nordholm, S. 2005 
  Spectral entropy as speech features for speech recognition.
  \emph{Proceedings of PEECS}, pp. 60-65.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{The spectral entropy of a noisy signal will tend towards 1 whereas
  the spectral entropy of a pure tone signal will tend towards 0.}

\seealso{\code{\link{sh}}, \code{\link{th}}}

\examples{
data(orni)
csh(orni,f=22050,wl=512,ovlp=50)
# using the threshold argument can lead to some edge effets
# here sh=1 at the end of echemes
csh(orni,f=22050,wl=512,ovlp=50,threshold=5)
}

\keyword{ts}