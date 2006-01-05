\name{dfreq}

\alias{dfreq}

\title{Dominant frequency of a time wave}

\description{
This function gives the dominant frequency (i. e. the frequency of highest amplitude) of a time wave.}

\usage{
dfreq(wave, f, wl, ovlp = 0, threshold = FALSE,
plot = TRUE, xlab = "Times (s)", ylab = "Frequency (kHz)",
ylim = c(0, f/2000), ...)
}

\arguments{
  \item{wave}{data describing a time wave.}
  \item{f}{sampling frequency of \code{wave}.}
  \item{wl}{length of the window for the analysis
    (even number of points).}
  \item{ovlp}{overlap between two successive analysis windows (in \% ).}
  \item{threshold}{amplitude threshold for signal detection (in \% ).}
  \item{plot}{logical, if \code{TRUE} plots the dominant frequency against time
  (by default \code{TRUE}).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{ylim}{the range of y values.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\value{
If \code{plot} is \code{FALSE}, \code{dfreq} function returns a vector of 
numeric data corresponding to the dominant frequency along the time wave.
}

\details{A Hanning function is applied to the analysis window.}

\author{Jérôme Sueur \email{jerome.sueur@ibaic.u-psud.fr}}

\note{This function is based on \code{\link{fft}}.} 

\seealso{\code{\link{spec}}, \code{\link{meanspec}},\code{\link{spectro}}.}

\examples{
data(tico)
dfreq(tico,f=22050,wl=512,ovlp=50,threshold=5)
# overlay on spectrogram
spectro(tico,f=22050,wl=512,ovlp=50,zp=16,scale=FALSE,
    collevels=seq(-40,0,1),palette=rev.terrain.colors)
par(new=TRUE)
dfreq(tico,f=22050,wl=512,ovlp=50,threshold=6,type="l",col="red",lwd=2,
    ann=FALSE,xaxs="i",yaxs="i")
}

\keyword{dplot}
\keyword{ts}