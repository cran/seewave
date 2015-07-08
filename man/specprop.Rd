\name{specprop}

\alias{specprop}

\title{Spectral properties}

\description{
  This function returns a list of statistical properties of a frequency spectrum.}

\usage{specprop(spec, f=NULL,
str = FALSE, flim=NULL, mel=FALSE,
plot = FALSE, type = "l", xlab=NULL, ylab = NULL,
col.mode = 2, col.quartiles = 4, ...)}

\arguments{
  \item{spec}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}
  \item{f}{sampling frequency of \code{spec} (in Hz).}
  \item{str}{logical, if \code{TRUE} returns the results in a structured table.}
  \item{flim}{a vector of length 2 to specifgy the frequency limits of the analysis (in kHz)}
  \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}
  \item{plot}{if \code{1} returns the spectrum , if \code{2} returns 
    the cumulative spectrum, both of them with the first quartile, the third quartile,
    the median and the mode plotted (by default \code{FALSE}).}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{xlab}{label of the x axis.}
  \item{ylab}{label of the y axis.}
  \item{col.mode}{colour of the mode segments (by default blue).}
  \item{col.quartiles}{colour of the quartiles segments (by default red).}
  \item{\dots}{other arguments to be passed to \code{plot}}
}

\details{
  The spectrum is converted in a probability mass function (PMF).\cr
  If a selected value has to be selected with $,
  the argument \code{str} has to be set to \code{FALSE}.
}

\value{
  A list of 15 values is returned
  \item{mean}{mean frequency (see \code{\link{mean}})}
  \item{sd}{standard deviation of the mean (see \code{\link{sd}})}
  \item{sem}{standard error of the mean}
  \item{median}{median frequency (see \code{\link{median}})}
  \item{mode}{mode frequency, \emph{i.e.} the dominant frequency}
  \item{Q25}{first quartile (see \code{\link{quantile}})}
  \item{Q75}{third quartile (see \code{\link{quantile}})}
  \item{IQR}{interquartile range (see \code{\link{IQR}})}
  \item{cent}{centroid, see \code{note}}
  \item{skewness}{skewness, a measure of asymmetry, see \code{note}}
  \item{kurtosis}{kurtosis, a measure of peakedness, see \code{note}}
  \item{sfm}{spectral flatness measure (see \code{\link{sfm}})}
  \item{sh}{spectral entropy (see \code{\link{sh}})}
  \item{prec}{frequency precision of the spectrum}
}

\note{
  Centroid is computed according to:\cr
  \deqn{C = \sum_{i = 1}^{N} x_i\times{y_i}}{%
    C = sum(x*y)}\cr
  with:\cr
  \emph{x} = frequencies,
  \emph{y} = relative amplitude of the \emph{i} frequency,\cr
  \emph{N} = number of frequencies.\cr\cr
  Skewness is computed according to:\cr
  \deqn{S = \frac{\sum_{i = 1}^{N}(x_i-\bar{x})^{3}}{N-1}\times{\frac{1}{\sigma^{3}}}}{%
    S = sum((x-mean(x))^3)/(N-1)/sd^3}.\cr
  S < 0 when the spectrum is skewed to left,\cr
  S = 0 when the spectrum is symetric,\cr
  S > 0 when the spectrum is skewed to right.\cr
  Spectrum asymmetry increases with |S|.\cr\cr 
  Kurtosis is computed according to:\cr
  \deqn{K = \frac{\sum_{i = 1}^{N}(x_i-\bar{x})^{4}}{N-1}\times{\frac{1}{\sigma^{4}}}}{%
    K = sum((x-mean(x))^4)/(N-1)/sd^4}.\cr
  K < 3 when the spectrum is platikurtic, \emph{i.e.} it has fewer items at the
  center and at the tails than the normal curve but has more items in the shoulders,\cr
  K = 3 when the spectrum shows a normal shape,\cr
  K > 3 when the spectrum is leptokurtic, \emph{i.e.} it has more items near
  the center and at the tails, with fewer items in the shoulders
  relative to normal distribution with the same mean and variance.
}

\author{Jerome Sueur and Caroline Simonis, and a patch by Jesse Ross
  (Dec. 2012)}

\examples{
data(orni)
a<-meanspec(orni,f=22050,plot=FALSE)
specprop(a,f=22050)
# to get a single measure of the list
specprop(a,f=22050)$mode
# to get the results structured
specprop(a,f=22050,str=TRUE)
# to limit the analysis between 4 and 6 kHz
specprop(a,f=22050,flim=c(4,6),str=TRUE)
# plots
specprop(a,f=22050,plot=1)
specprop(a,f=22050,plot=2)
# (htk-)mel scale
require(tuneR)
mel <- melfcc(orni, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
melspec.mean <- apply(mel$aspectrum, MARGIN=2, FUN=mean)
specprop(melspec.mean, f=22050, mel=TRUE)
# be aware that flim is always given in kHz even if mel=TRUE
specprop(melspec.mean, f=22050, flim=c(4,6), mel=TRUE, plot=TRUE)
}

\keyword{dplot}
\keyword{ts}
