\name{csh}

\alias{csh}

\title{Continuous spectral entropy}

\description{This function computes the continuous spectral entropy (H) of a time wave.}

\usage{
csh(wave, f, wl = 512, wn = "hanning", ovlp = 0, threshold = FALSE,
plot = TRUE, xlab = "Times (s)", ylab = "Spectral Entropy",
ylim = c(0, 1.1), ...)
}

\arguments{
  \item{wave}{data describing a time wave or a \code{\link[Sound]{Sample}}
    object created loading a wav file with \code{\link[Sound]{loadSample}}
    (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz)}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{threshold}{amplitude threshold for signal detection (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the spectral entropy against time
    (by default \code{TRUE}).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{ylim}{the range of y values.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{See \code{\link{sh}} for computing method.}

\value{If \code{plot} is \code{FALSE}, \code{csh} function returns a vector of 
numeric data corresponding to the spectral entropy along the time wave.}

\references{Toh, A. M., Togneri, R. & Nordholm, S. 2005 
Spectral entropy as speech features for speech recognition.
\emph{Proceedings of PEECS}, pp. 60-65.}

\author{Jérôme Sueur \email{jerome.sueur@univ-tours.fr}}

\note{The spectral entropy of a noisy signal will tend towards 1 whereas
the spectral entropy of a pure tone signal will tend towards 0.}

\seealso{\code{\link{sh}}}

\examples{
data(orni)
csh(orni,f=22050,wl=512,ovlp=50,type="l")
# using the threshold argument can lead to some edge effets
# here sh=1 at the end of echemes
csh(orni,f=22050,wl=512,ovlp=50,threshold=5,type="l")
}

\keyword{ts}