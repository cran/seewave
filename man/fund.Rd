\name{fund}
\alias{fund}

\title{Fundamental frequency track}

\description{This function tracks the fundamental frequency through a
short-term cepstral transform.
}

\usage{fund(wave, f, wl = 512, ovlp = 0, fmax, threshold = NULL,
plot = TRUE, xlab = "Time (s)", ylab = "Frequency (kHz)",
ylim = c(0, f/2000), ...)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{fmax}{the maximum frequency to detect (in Hz).}  
  \item{threshold}{amplitude threshold for signal detection (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the fundamental frequency
  modulations against time (by default \code{TRUE}).}
  \item{xlab}{title of the time axis (s).}
  \item{ylab}{title of the frequency axis (Hz).}
  \item{ylim}{the range of frequency values.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\value{
When \code{plot} is \code{FALSE}, \code{fund} returns a two-column matrix, the first column corresponding to time in seconds (\emph{x}-axis) and the second column corresponding to 
to fundamental frequency in kHz (\emph{y}-axis).\cr
NA corresponds to pause sections in \code{wave} (see \code{threshold}). 
}

\references{Oppenheim, A.V. and Schafer, R.W. 2004. From frequency
to quefrency: a history of the cepstrum. \emph{Signal Processing Magazine IEEE},
21: 95-106.}

\author{Jerome Sueur \email{sueur@mnhn.fr}.}

\note{This function is based on \code{\link{ceps}}.} 

\seealso{\code{\link{cepstro}}, \code{\link{ceps}}, \code{\link{autoc}}}

\examples{
data(sheep)
fund(sheep,f=8000,fmax=300,type="l")
# with 50% overlap between successive sliding windows, time zoom and 
# amplitude filter (threshold)
fund(sheep,f=8000,fmax=300,type="b",ovlp=50,threshold=5,ylim=c(0,1),cex=0.5)
# overlaid on a spectrogram
spectro(sheep,f=8000,ovlp=75,zp=16,scale=FALSE,palette=rev.gray.colors.2)
par(new=TRUE)
fund(sheep,f=8000,fmax=300,type="p",pch=24,ann=FALSE,
  xaxs="i",yaxs="i",col="black",bg="red",threshold=6)
}

\keyword{dplot}
\keyword{ts}
