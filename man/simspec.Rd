\name{simspec}

\alias{simspec}

\title{Similarity between two frequency spectra}

\description{This function estimates the similarity between
  two frequency spectra.}

\usage{
simspec(spec1, spec2, f = NULL, mel = FALSE,
norm = FALSE, PMF = FALSE,
plot = FALSE, type = "l",
lty =c(1, 2, 3), col = c(2, 4, 1),
flab = NULL, alab = "Amplitude (percentage)",
flim = NULL, alim = NULL,
title = TRUE, legend = TRUE, ...)
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
    is a two columns matrix obtained with \code{\link{spec}}
    or \code{\link{meanspec}}.}
  \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}  
  \item{norm}{a logical, if \code{TRUE} \code{spec1} and \code{spec2}
  are normalised (scaled) between 0 and 1.}
  \item{PMF}{a logical, if \code{TRUE} \code{spec1} and \code{spec2} are
  transformed into probability mass functions.}
  \item{plot}{logical, if \code{TRUE} plots both spectra and similarity function 
    (by default \code{FALSE}).}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{lty}{a vector of length 3 for the line type of \code{spec1},
    \code{spec2} and  of the similarity function if \code{type}="l".}
  \item{col}{a vector of length 3 for the colour of \code{spec1},
   \code{spec2}, and the similarity function.}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{the range of frequency values.}
  \item{alim}{range of amplitude axis.}
  \item{title}{logical, if \code{TRUE}, adds a title with S value.}
  \item{legend}{logical, if \code{TRUE} adds a legend to the plot.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  Spectra similarity is assessed according to:\cr
  \deqn{S = \frac{100/N} \times{\sum_{i=1}^N{\frac{\min{spec1(i),spec2(i)}}
        {\max{spec1(i),spec2(i)}}}}}{% 
    S = 100*sum(min(spec1,spec2)/max(spec1,spec2))/N, with \emph{S} in \%.}
  with \emph{S} in \%.
}

\value{
  The similarity index is returned. This value is in \%.\cr
  When \code{plot} is \code{TRUE}, both spectra and the similarity function are
  plotted on the same graph. The similarity index is the mean of this function.
}

\references{Deecke, V. B. and Janik, V. M. 2006. Automated categorization of
  bioacoustic signals: avoiding perceptual pitfalls. \emph{Journal of the Acoustical
    Society of America}, 119: 645-653.}

\author{Jerome Sueur, improved by Laurent Lellouch}

\seealso{\code{\link{spec}}, \code{\link{meanspec}}, \code{\link{corspec}},
  \code{\link{diffspec}}, \code{\link{diffenv}},  \code{\link{kl.dist}},
   \code{\link{ks.dist}}, \code{\link{logspec.dist}}, \code{\link{itakura.dist}}}

\examples{
a<-noisew(f=8000,d=1)
b<-synth(f=8000,d=1,cf=2000)
c<-synth(f=8000,d=1,cf=1000)
d<-noisew(f=8000,d=1)
speca<-spec(a,f=8000,at=0.5,plot=FALSE)
specb<-spec(b,f=8000,at=0.5,plot=FALSE)
specc<-spec(c,f=8000,at=0.5,plot=FALSE)
specd<-spec(d,f=8000,at=0.5,plot=FALSE)
simspec(speca,speca)
simspec(speca,specb)
simspec(speca,specc,plot=TRUE)
simspec(specb,specc,plot=TRUE)
#[1] 12.05652
simspec(speca,specd,plot=TRUE)
## mel scale
require(tuneR)
data(orni)
data(tico)
orni.mel <- melfcc(orni, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
orni.mel.mean <- apply(orni.mel$aspectrum, MARGIN=2, FUN=mean)
tico.mel <- melfcc(tico, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
tico.mel.mean <- apply(tico.mel$aspectrum, MARGIN=2, FUN=mean)
simspec(orni.mel.mean, tico.mel.mean, f=22050, mel=TRUE, plot=TRUE)
}

\keyword{dplot}
\keyword{ts}
