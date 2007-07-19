\name{simspec}

\alias{simspec}

\title{Similarity between two frequency spectra}

\description{This function estimates the similarity between
 two frequency spectra.}

\usage{
simspec(spec1, spec2, f, plot = FALSE, type = "l", lty1 = 1, lty2 = 2,
lty3 = 3, col1 = 2, col2 = 4, col3 = 1, flab = "Frequency (kHz)",
alab = "Amplitude (percentage)", flim = c(0, f/2000), alim = c(0, 100), ...)
}

\arguments{
  \item{spec1}{a first data set resulting of a spectral analysis obtained
      with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}
  \item{spec2}{a second data set resulting of a spectral analysis obtained
      with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}
  \item{f}{sampling frequency of waves used to obtain \code{spec1}
      and \code{spec2} (in Hz).}
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
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
Spectra similarity is assessed according to:\cr
\deqn{S = 100\times{\sum_{i=1}^N{\frac{\min{spec1(i),spec2(i)}}
        {\max{spec1(i),spec2(i)}}}}}{% 
      S = 100*sum(min(spec1,spec2)/max(spec1,spec2)/N, with \emph{S} in \%.}
      with \emph{S} in \%.
}

\value{
The similarity index is returned. This value is in \%.\cr
When \code{plot} is \code{TRUE}, both spectra and the similarity function are
plotted on the same graph. The index is the mean of this function.
}

\references{Deecke, V. B. and Janik, V. M. 2006. Automated categorization of
bioacoustic signals: avoiding perceptual pitfalls. \emph{Journal of the Acoustical
Society of America}, 119: 645-653.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\seealso{\code{\link{spec}}, \code{\link{meanspec}}, \code{\link{corspec}},
\code{\link{diffspec}}, \code{\link{diffenv}}}

\examples{
a<-noise(f=8000,d=1)
b<-synth(f=8000,d=1,cf=2000)
c<-synth(f=8000,d=1,cf=1000)
d<-noise(f=8000,d=1)
speca<-spec(a,f=8000,wl=512,at=0.5,plot=FALSE)
specb<-spec(b,f=8000,wl=512,at=0.5,plot=FALSE)
specc<-spec(c,f=8000,wl=512,at=0.5,plot=FALSE)
specd<-spec(d,f=8000,wl=512,at=0.5,plot=FALSE)
simspec(speca,speca,f=8000)
#[1] 100 => similar spectra of course!
simspec(speca,specb,f=8000)
#[1] 0.3873278
simspec(speca,specc,f=8000,plot=TRUE)
#[1] 0.5920857
# two pure tone at different frequency are...
simspec(specb,specc,f=8000,plot=TRUE)
#[1] 12.05652
#... less similar than two noises
simspec(speca,specd,f=8000,plot=TRUE)
#[1] 57.40742
}

\keyword{dplot}
\keyword{ts}
