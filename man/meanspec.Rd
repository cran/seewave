\name{meanspec}

\alias{meanspec}

\title{Mean frequency spectrum of a time wave}

\description{
This function returns the mean frequency spectrum (i.e. the mean relative amplitude of the
frequency distribution) of a time wave. Results can be expressed either in absolute 
or dB data.
}

\usage{
meanspec(wave, f, wl = 512, wn = "hanning", ovlp = 0, PSD = FALSE, PMF = FALSE,
dB = FALSE, from = FALSE, to = FALSE, peaks = FALSE, identify = FALSE,
col = "black", cex = 1, colpeaks = "red", cexpeaks = 1, 
fontpeaks = 1, plot = 1, flab = "Frequency (kHz)",
alab = "Amplitude", flim = c(0,f/2000),alim = NULL, type ="l",...)
}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{wl}{length of the window for the analysis (even number of points, by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{ovlp}{overlap between two successive analysis windows (in \%).}
  \item{PSD}{if \code{TRUE} return Power Spectra Density,
      \emph{i. e.} the square of the spectra.}   
  \item{PMF}{if \code{TRUE} return Probability Mass Function,
      \emph{i. e.} the probability distribution of frequencie.}  
  \item{dB}{logical, if \code{TRUE} return the spectrum in dB 
      (by default \code{FALSE}).}
  \item{from}{start mark where to compute the spectrum (in s).}
  \item{to}{end mark where to compute the spectrum (in s).}
  \item{peaks}{to return peaks values.}
  \item{identify}{to identify frequency and amplitude values on the plot
      with the help of a cursor.}
  \item{col}{colour of the spectrum.}
  \item{cex}{pitch size.}
  \item{colpeaks}{colour of peaks value plotted on the spectrum.}
  \item{cexpeaks}{character size of peaks value plotted on the spectrum.}
  \item{fontpeaks}{font of peaks value plotted on the spectrum.}
  \item{plot}{if \code{1} returns frequency on x-axis, if \code{2} returns 
  frequency on y-axis, (by default \code{1}).}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{range of frequency axis.}
  \item{alim}{range of amplitude axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
      See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
If \code{from} and \code{to} are \code{FALSE} then \code{spec} 
computes the mean spectrum of the whole signal.\cr
\code{peaks} setting corresponds to \code{dimension} of \code{\link{embed}}.\cr
When using \code{identify}, press \sQuote{stop} tools bar button to return values in the console.  
}

\value{
If \code{plot} is \code{FALSE}, \code{meanspec} returns  a vector of numeric
data corresponding to the computed spectrum.
If \code{identify} is \code{TRUE}, \code{spec} returns a list with
two elements:\cr 
  \item{freq}{the frequency of the points chosen on the spectrum}\cr 
  \item{amp}{the relative amplitude of the points chosen on the spectrum}\cr 
If \code{peaks} is \code{TRUE}, \code{spec} returns a list with
two elements:\cr 
  \item{spec}{the spectrum computed}\cr 
  \item{peaks }{the peaks values (in kHz).}\ 
}

\author{Jérôme Sueur \email{sueur@mnhn.fr} and 
Martin Maechler \email{maechler@stat.math.ethz.ch} for \code{peaks}}

\note{See examples of \code{\link{spec}}. This function is based on \code{\link{fft}}.} 

\seealso{\code{\link{spec}}, \code{\link{corspec}}, \code{\link{fft}}.}

\examples{
data(orni)
# compute the mean spectrum of the whole time wave
meanspec(orni,f=22050)
# compute the mean spectrum of a time wave section (from 0.32 s to 0.39 s)
meanspec(orni,f=22050,from=0.32,to=0.39)
# different window lengths
op<-par(mfrow=c(3,1))
meanspec(orni,f=22050,wl=256)
title("wl=256")
meanspec(orni,f=22050,wl=1024)
title("wl=1024")
meanspec(orni,f=22050,wl=4096)
title("wl=4096")
par(op)
# different overlap values (almost no effects here...)
op<-par(mfrow=c(3,1))
meanspec(orni,f=22050)
title("ovlp=0")
meanspec(orni,f=22050,ovlp=50)
title("ovlp=50")
meanspec(orni,f=22050,ovlp=95)
title("ovlp=95")
par(op)
# use of flim to zoom in
op<-par(mfrow=c(2,1))
meanspec(orni,f=22050)
title("zoom in")
meanspec(orni,f=22050,wl=512,flim=c(4,6))
par(op)
# comparaison of spectrum and mean spectrum
op<-par(mfrow=c(2,1))
spec(orni,f=22050)
title("spec()")
meanspec(orni,f=22050)
title("meanspec()")
par(op)
}

\keyword{dplot}
\keyword{ts}
