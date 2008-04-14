\name{autoc}

\alias{autoc}

\title{Short-term autocorrelation of a time wave}

\description{
This function returns the fundamental frequency of a harmonic time wave. This is
achieved by computing a correlation of the signal with itself
after a time delay. 
}

\usage{
autoc(wave, f, wl = 512, fmax, threshold = NULL, plot = TRUE,
xlab = "Time (s)", ylab = "Frequency (kHz)", ylim = c(0, f/2000),...)
}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{wl}{length of the window for the analysis (even number of points,
  by default = 512).}
  \item{fmax}{the maximum frequency to detect (in Hz).}
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
When \code{plot} is \code{FALSE}, \code{autoc} returns a two-column matrix, the first column corresponding to time in seconds (\emph{x}-axis) and the second column corresponding to 
to fundamental frequency in kHz (\emph{y}-axis).\cr
NA corresponds to pause sections in \code{wave} (see \code{threshold}). 
}

\references{
Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg.}

\author{Jerome Sueur \email{sueur@mnhn.fr} and 
Thierry Aubin \email{thierry.aubin@u-psud.fr}}

\seealso{\code{\link{ceps}}}

\examples{
data(sheep)
# fundamental frequency of a sheep
autoc(sheep,f=8000,threshold=5,fmax=700)
# overlay on spectrogram
spectro(sheep,f=8000,ovlp=75,zp=16,scale=FALSE)
par(new=TRUE)
autoc(sheep,f=8000,wl=512,threshold=5,fmax=700,col="black",pch=20,
    xaxs="i",yaxs="i",ann=FALSE,yaxt="n")
legend(0.5,3.6,"Fundamental frequency",pch=20,col="black",bty=0,cex=0.7)
}

\keyword{dplot}
\keyword{ts}
