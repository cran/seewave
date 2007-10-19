\name{diffspec}

\alias{diffspec}

\title{Difference between two frequency spectra}

\description{This function estimates the surface difference between
 two frequency spectra.}

\usage{
diffspec(spec1, spec2, f, dB = FALSE, plot = FALSE, type="l", 
lty1 = 1, lty2 = 2, col1 = 2, col2 = 4, cold = 8,
flab = "Frequency (kHz)", alab = "Amplitude",
flim = c(0, f/2000), alim = c(0,max(spec1,spec2)),...)
}

\arguments{
  \item{spec1}{a first data set resulting of a spectral analysis obtained
      with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}
  \item{spec2}{a second data set resulting of a spectral analysis obtained
      with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}
  \item{f}{sampling frequency of waves used to obtain \code{spec1}
      and \code{spec2} (in Hz).}
  \item{dB}{logical, if \code{TRUE} return the spectra and their surface
      difference in dB (by default \code{FALSE}).}
  \item{plot}{logical, if \code{TRUE} plots both spectra and their surface
      difference (by default \code{FALSE}).}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
      See \code{\link{plot}} for details (by default "l" for lines).}
  \item{lty1}{line type of \code{spec1} if \code{type}="l".}
  \item{lty2}{line type of \code{spec2} if \code{type}="l".}
  \item{col1}{colour of \code{spec1}.}
  \item{col2}{colour of \code{spec2}.}
  \item{cold}{colour of the surface difference.}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{the range of frequency values.}
  \item{alim}{range of amplitude axis.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
Both spectra are first transformed as probability mass functions (PMF).\cr
Spectral difference is then computed according to:\cr
\deqn{D = \sum{|spec1-spec2|}, with D \in [0,1].}{% 
      D = sum(abs(spec1-spec2)), with 0 < D <1.}
}

\value{
The difference is returned. This value is without unit.
If \code{dB} is \code{TRUE},  the same value is returned in dB.\cr
When \code{plot} is \code{TRUE}, both spectra and their difference surface are
plotted on the same graph.
}

\author{Jérôme Sueur \email{sueur@mnhn.fr}
and Sandrine Pavoine \email{pavoine@mnhn.fr}.}

\note{This method can be used as a relative distance estimation
between different spectra.\cr
The dB value obtained can be very different from the one visually estimated
when looking at the graph (\code{plot=TRUE}).}

\seealso{\code{\link{spec}}, \code{\link{meanspec}}, \code{\link{corspec}},
\code{\link{simspec}}, \code{\link{diffenv}}}

\examples{
a<-noise(f=8000,d=1)
b<-synth(f=8000,d=1,cf=2000)
c<-synth(f=8000,d=1,cf=1000)
d<-noise(f=8000,d=1)
speca<-spec(a,f=8000,wl=512,at=0.5,plot=FALSE)
specb<-spec(b,f=8000,wl=512,at=0.5,plot=FALSE)
specc<-spec(c,f=8000,wl=512,at=0.5,plot=FALSE)
specd<-spec(d,f=8000,wl=512,at=0.5,plot=FALSE)
diffspec(speca,speca,f=8000)
#[1] 0 => similar spectra of course !
diffspec(speca,specb,f=8000)
#[1] 0.965465
diffspec(speca,specc,f=8000,plot=TRUE)
#[1] 0.95498
# two pure tone at different frequency are...
diffspec(specb,specc,f=8000,plot=TRUE)
#[1] 0.9999166
#... more different than two noises
diffspec(speca,specd,f=8000,plot=TRUE)
#[1] 0.2810392
}

\keyword{dplot}
\keyword{ts}

