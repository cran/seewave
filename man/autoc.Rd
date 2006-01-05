\name{autoc}

\alias{autoc}

\title{Short-term autocorrelation of a time wave}

\description{
This function returns the fundamental frequency of a harmonic time wave. This is
achieved by computing a correlation of the signal with itself
after a time delay. 
}

\usage{
autoc(wave, f, wl, fmin, threshold = FALSE, plot = TRUE,
xlab = "Time (s)", ylab = "Frequency (kHz)", ylim = c(0, f/2000), ...)
}

\arguments{
  \item{wave}{data describing a time wave.}
  \item{f}{sampling frequency of \code{wave}.}
  \item{wl}{length of the window for the analysis (even number of points).}
  \item{fmin}{the minimum frequency to detect (in Hz).}  
  \item{threshold}{amplitude threshold for signal detection (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the fundamental frequency 
  against time (by default \code{TRUE}).}
  \item{xlab}{title of the x-axis.}
  \item{ylab}{title of the y-axis.}
  \item{ylim}{the range of y values.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{Autocorrelation process can be time consuming.}

\value{
If \code{plot} is \code{FALSE}, \code{autoc} returns a vector of numeric data 
describing the fundamental frequency in kHz. NA corresponds 
to pause sections in the signal (see \code{threshold}). 
}

\references{
Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg.}

\author{Jérôme Sueur \email{jerome.sueur@ibaic.u-psud.fr} and 
Thierry Aubin \email{thierry.aubin@ibaic.u-psud.fr}}

\examples{
data(sheep)
# fundamental frequency of a sheep
autoc(sheep,f=8000,wl=512,t=5,fmin=700)
# overlay on spectrogram
spectro(sheep,f=8000,wl=512,ovlp=75,zp=16,scale=FALSE)
par(new=TRUE)
autoc(sheep,f=8000,wl=512,t=5,fmin=700,type="p",col="black",pch=20,
    xaxs="i",yaxs="i",ann=FALSE,yaxt="n")
legend(0.5,3.6,"Fundamental frequency",pch=20,col="black",bty=0,cex=0.7)
}

\keyword{dplot}
\keyword{ts}
