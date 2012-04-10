\name{simspec}

\alias{simspec}

\title{Similarity between two frequency spectra}

\description{This function estimates the similarity between
  two frequency spectra.}

\usage{
simspec(spec1, spec2, f = NULL, plot = FALSE, type = "l",
lty1 = 1, lty2 = 2,
lty3 = 3, col1 = 2, col2 = 4, col3 = 1, flab = "Frequency (kHz)",
alab = "Amplitude (percentage)", flim = c(0, f/2000),
alim = c(0, 100),
legend = TRUE, ...)
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
  \item{plot}{logical, if \code{TRUE} plots both spectra and similarity function 
    (by default \code{FALSE}).}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{lty1}{line type of \code{spec1} if \code{type}="l".}
  \item{lty2}{line type of \code{spec2} if \code{type}="l".}
  \item{lty3}{line type of the similarity function if \code{type}="l".}
  \item{col1}{colour of \code{spec1}.}
  \item{col2}{colour of \code{spec2}.}
  \item{col3}{colour of the similarity function.}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{the range of frequency values.}
  \item{alim}{range of amplitude axis.}
  \item{legend}{logical, if \code{TRUE} adds a legend to the plot.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  Spectra similarity is assessed according to:\cr
  \deqn{S = 100\times{\sum_{i=1}^N{\frac{\min{spec1(i),spec2(i)}}
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

\author{Jerome Sueur \email{sueur@mnhn.fr}}

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
}

\keyword{dplot}
\keyword{ts}
