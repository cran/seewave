\name{diffspec}

\alias{diffspec}

\title{Difference between two frequency spectra}

\description{This function estimates the surface difference between
  two frequency spectra.}

\usage{
diffspec(spec1, spec2, f = NULL, mel = FALSE,
plot = FALSE, type="l", 
lty=c(1, 2), col =c(2, 4, 8),
flab = NULL, alab = "Amplitude",
flim = NULL, alim = NULL, title = TRUE, legend = TRUE, ...)
}

\arguments{
  \item{spec1}{a first data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This
    can be either a two-column matrix (col1 = frequency, col2 = amplitude)
    or a vector (amplitude).}
  \item{spec2}{a first data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This
    can be either a two-column matrix (col1 = frequency, col2 = amplitude)
    or a vector (amplitude).}
  \item{f}{sampling frequency of waves used to obtain \code{spec1}
    and \code{spec2} (in Hz). Not necessary if \code{spec1} and/or \code{spec2}
    is a two-column matrix obtained with \code{\link{spec}} or \code{\link{meanspec}}.}
  \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}
  \item{plot}{logical, if \code{TRUE} plots both spectra and their surface
    difference (by default \code{FALSE}).}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{lty}{a vector of length 2 for the line type of \code{spec1} and
    \code{spec2} if \code{type}="l".}
  \item{col}{a vector of length 3 for the colour of \code{spec1},
    \code{spec2}, and the surface difference between each of them.}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{the range of frequency values.}
  \item{alim}{range of amplitude axis.}
  \item{title}{logical, if \code{TRUE}, adds a title with D value.}
  \item{legend}{logical, if \code{TRUE} adds a legend to the plot.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  Both spectra are first transformed as probability mass functions (PMF).\cr
  Spectral difference is then computed according to:\cr
  \deqn{D = \frac{\sum{|spec1-spec2|}}{2}, with D \in [0,1].}{% 
    D = sum(abs(spec1-spec2))/2}, with 0 < D < 1.
}

\value{
  The difference is returned. This value is without unit.
  When \code{plot} is \code{TRUE}, both spectra and their difference surface are
  plotted on the same graph.
}

\references{
  Sueur, J., Pavoine, S., Hamerlynck, O. and Duvail, S. (2008).
  Rapid acoustic survey for biodiversity appraisal. \emph{PLoS One},
  3(12): e4065.}

\author{Jerome Sueur, Sandrine Pavoine and Laurent Lellouch}

\note{This method can be used as a relative distance estimation
  between different spectra.\cr
  The dB value obtained can be very different from the one visually estimated
  when looking at the graph (\code{plot=TRUE}).}

\seealso{\code{\link{spec}}, \code{\link{meanspec}}, \code{\link{corspec}},
  \code{\link{simspec}}, \code{\link{diffcumspec}}, \code{\link{diffenv}},  \code{\link{kl.dist}},
  \code{\link{ks.dist}}, \code{\link{logspec.dist}}, \code{\link{itakura.dist}}}

\examples{
a <- noisew(f=8000,d=1)
b <- synth(f=8000,d=1,cf=2000)
c <- synth(f=8000,d=1,cf=1000)
d <- noisew(f=8000,d=1)
speca <- spec(a,f=8000,wl=512,at=0.5,plot=FALSE)
specb <- spec(b,f=8000,wl=512,at=0.5,plot=FALSE)
specc <- spec(c,f=8000,wl=512,at=0.5,plot=FALSE)
specd <- spec(d,f=8000,wl=512,at=0.5,plot=FALSE)
diffspec(speca,speca,f=8000)
#[1] 0 => similar spectra of course !
diffspec(speca,specb)
diffspec(speca,specc,plot=TRUE)
diffspec(specb,specc,plot=TRUE)
diffspec(speca,specd,plot=TRUE)
## mel scale
require(tuneR)
data(orni)
data(tico)
orni.mel <- melfcc(orni, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
orni.mel.mean <- apply(orni.mel$aspectrum, MARGIN=2, FUN=mean)
tico.mel <- melfcc(tico, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
tico.mel.mean <- apply(tico.mel$aspectrum, MARGIN=2, FUN=mean)
diffspec(orni.mel.mean, tico.mel.mean, f=22050, mel=TRUE, plot=TRUE)
}

\keyword{dplot}
\keyword{ts}

