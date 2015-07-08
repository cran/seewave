\name{diffcumspec}

\alias{diffcumspec}

\title{Difference between two cumulative frequency spectra}

\description{This function compares two distributions (e.g. two frequency spectra) by
  computing the difference between two cumulative frequency spectra}

\usage{
diffcumspec(spec1, spec2, f = NULL, mel = FALSE, 
plot = FALSE, type = "l", lty = c(1, 2), col = c(2, 4, 8),
flab = NULL, alab = "Cumulated amplitude",
flim = NULL, alim = NULL,
title = TRUE, legend = TRUE, ...)
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
 \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}
 \item{plot}{logical, if \code{TRUE} plots both cumulative spectra and
    their distance.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{col}{a vector of length 3 for the colour of \code{spec1},
  \code{spec2}, and the difference between each of them.}
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

\details{
  Both spectra are transformed into cumulative distribution functions
  (CDF). \cr
  Spectral difference is then computed according to:\cr
  \deqn{D_{cf}(x, y) = \frac{\sum_{i=1}^{n}|X_{i} - Y_{i}|}{n}, with
  with X and Y the spectrum CDFs, and D \in [0,1].}{%
    D = 1/n * sum(abs(X-Y)), with X and Y the spectrum CDFs, and 0 < D
  < 1.}
}

\value{A numeric vector of length 1 returning the difference between the
  two spectra. No unit.}


\note{This metric is sensitive not only to the spectral overlap between
  but also to the mean frequential distance between the
  different frequency peaks.}

\references{Lellouch L, Pavoine S, Jiguet F, Glotin H, Sueur J (2014) Monitoring temporal change of bird communities with dissimilarity acoustic indices. \emph{Methods in Ecology and Evolution}, in press.
}

\author{Laurent Lellouch, Jerome Sueur}

\seealso{\code{\link{kl.dist}}, \code{\link{ks.dist}}, \code{\link{simspec}},
  \code{\link{diffspec}}, \code{\link{logspec.dist}}, \code{\link{itakura.dist}}}

\examples{
## Hz scale
data(tico)
data(orni)
orni.hz <- meanspec(orni, plot=FALSE)
tico.hz <- meanspec(tico, plot=FALSE)
diffcumspec(orni.hz, tico.hz, plot=TRUE)
## mel scale
require(tuneR)
orni.mel <- melfcc(orni, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
orni.mel.mean <- apply(orni.mel$aspectrum, MARGIN=2, FUN=mean)
tico.mel <- melfcc(tico, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
tico.mel.mean <- apply(tico.mel$aspectrum, MARGIN=2, FUN=mean)
diffcumspec(orni.mel.mean, tico.mel.mean, f=22050, mel=TRUE, plot=TRUE)
}

\keyword{distribution}
\keyword{ts}
