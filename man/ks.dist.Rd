\name{ks.dist}

\alias{ks.dist}

\title{Kolmogorov-Smirnov distance}

\description{Compare two distributions (e.g. two frequency spectra) by
  computing the Kolmogorov-Smirnov distance}

\usage{
ks.dist(spec1, spec2, f = NULL, plot = FALSE, type = "l",
col = c("blue", "red"),
lty = c(2, 4), flab = "Frequency (kHz)",
alab = "Cumulated amplitude", flim = c(0, f/2000),
alim = c(0, 1), title = TRUE, legend = TRUE, ...)
}

\arguments{
  \item{spec1}{any distribution, especially a spectrum obtained with \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This can be either a two-column matrix (col1 = frequency, col2 = amplitude) or a vector (amplitude).}
  \item{spec2}{any distribution, especially a spectrum obtained with
    \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This can be
    either a two-column matrix (col1 = frequency, col2 = amplitude) or a
    vector (amplitude).} 
  \item{f}{sampling frequency of waves used to obtain \code{spec1}
    and \code{spec2} (in Hz). Not necessary if \code{spec1} and/or \code{spec2}
    is a two columns matrix obtained with \code{\link{spec}}
    or \code{\link{meanspec}}.}
  \item{plot}{logical, if \code{TRUE} plots both cumulated spectra and
    their maximal distance (\emph{i.e.} the K-S distance.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{col}{a vector of length 2 for the colour of \code{spec1} and \code{spec2}.}
  \item{lty}{a vector of length 2 for the line type of \code{spec1} and
    \code{spec2} if \code{type}="l".}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{the range of frequency values.}
  \item{alim}{range of amplitude axis.}
  \item{title}{logical, if \code{TRUE}, adds a title with D and F values.}
  \item{legend}{logical, if \code{TRUE} adds a legend to the plot.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{The Kolmogorov distance is the maximal distance between the
  cumulated spectra. The function returns this distance and the
  corresponding frequency. This is an adaptation of the statistic
  computed by the  non-parametric  Kolmogorov-Smirnov test (see \code{\link{ks.test}}).}

\value{The function returns a list of two items
  \item{D}{the Kolomogorv-Smirnov distance}
  \item{F}{the frequency (in KHz) where the Kolmogorov-Smirnov distance was found}
}

\author{Jerome Sueur}

\note{There is no p-value associated to the K-S distance.}

\seealso{\code{\link{kl.dist}}, \code{\link{simspec}}, \code{\link{diffspec}}}

\examples{
# Comparison of two spectra and plot of the cumulated spectra with the K-S distance
data(tico)
tico1 <- spec(tico, at=0.65, plot=FALSE)
tico2 <- spec(tico, at=1.1, plot=FALSE)
ks.dist(tico1, tico2, plot=TRUE)
}

\keyword{distrubtion}
\keyword{ts}
