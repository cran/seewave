\name{dfreq}

\alias{dfreq}

\title{Dominant frequency of a time wave}

\description{
This function gives the dominant frequency (i. e. the frequency of highest amplitude) of a time wave.}

\usage{
dfreq(wave, f, wl = 512, wn = "hanning", ovlp = 0, threshold = NULL,
plot = TRUE, xlab = "Times (s)", ylab = "Frequency (kHz)",
ylim = c(0, f/2000), type ="l",...)
}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}  \item{wl}{length of the window for the analysis
    (even number of points, by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{ovlp}{overlap between two successive analysis windows (in \% ).}
  \item{threshold}{amplitude threshold for signal detection (in \% ).}
  \item{plot}{logical, if \code{TRUE} plots the dominant frequency against time
  (by default \code{TRUE}).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{ylim}{the range of y values.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
      See \code{\link{plot}} for details (by default "l" for lines).}  
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\value{
When \code{plot} is \code{FALSE}, \code{dfreq} returns a two-column matrix, the first column corresponding to time in seconds (\emph{x}-axis) and the second column corresponding to 
to dominant frequency in kHz (\emph{y}-axis).\cr
NA corresponds to pause sections in \code{wave} (see \code{threshold}). 
}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{This function is based on \code{\link{fft}}.} 

\seealso{\code{\link{spec}}, \code{\link{meanspec}},\code{\link{spectro}}.}

\examples{
data(tico)
dfreq(tico,f=22050,ovlp=50,threshold=5)
# overlay on spectrogram
spectro(tico,f=22050,ovlp=50,zp=16,scale=FALSE,
    collevels=seq(-40,0,1),palette=rev.terrain.colors)
par(new=TRUE,las=1)
dfreq(tico,f=22050,ovlp=50,threshold=6,col="red",lwd=2,
    ann=FALSE)
}

\keyword{dplot}
\keyword{ts}
