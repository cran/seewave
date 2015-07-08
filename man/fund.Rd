\name{fund}

\alias{fund}

\title{Fundamental frequency track}

\description{This function estimates the fundamental frequency through a
  short-term cepstral transform.
}

\usage{fund(wave, f, wl = 512, ovlp = 0, fmax, threshold = NULL,
at = NULL, from = NULL, to = NULL, 
plot = TRUE, xlab = "Time (s)", ylab = "Frequency (kHz)",
ylim = c(0, f/2000), pb = FALSE, ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{fmax}{the maximum frequency to detect (in Hz).}  
  \item{threshold}{amplitude threshold for signal detection (in \%).}
  \item{at}{position where the estimate the fundamental frequency (in s)}.
  \item{from}{start position where to compute the fundamental frequency (in s).}
  \item{to}{end position to compute the fundamental frequency (in s).}
  \item{plot}{logical, if \code{TRUE} plots the fundamental frequency
  modulations against time (by default \code{TRUE}).}
  \item{xlab}{title of the time axis (s).}
  \item{ylab}{title of the frequency axis (Hz).}
  \item{ylim}{the range of frequency values.}
  \item{pb}{if \code{TRUE} returns a text progress bar in the console.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\value{
When \code{plot} is \code{FALSE}, \code{fund} returns a two-column matrix, the first column corresponding to time in seconds (\emph{x}-axis) and the second column corresponding to 
to fundamental frequency in kHz (\emph{y}-axis).\cr
NA corresponds to pause sections in \code{wave} (see
\code{threshold}). \cr
No plot is produced when using \code{at}.
}

\references{Oppenheim, A.V. and Schafer, R.W. 2004. From frequency
to quefrency: a history of the cepstrum. \emph{Signal Processing Magazine IEEE},
21: 95-106.}

\author{Jerome Sueur \email{sueur@mnhn.fr}.}

\note{This function is based on \code{\link{ceps}}.} 

\seealso{\code{\link{cepstro}}, \code{\link{ceps}}, \code{\link{autoc}}}

\examples{
data(sheep)
# estimate the fundamental frequency at a single position
fund(sheep, f=8000, fmax=300, at=1, plot=FALSE)
# track the fundamental frequency along time
fund(sheep,f=8000,fmax=300,type="l")
# with 50% overlap between successive sliding windows, time zoom and 
# amplitude filter (threshold)
fund(sheep,f=8000,fmax=300,type="b",ovlp=50,threshold=5,ylim=c(0,1),cex=0.5)
# overlaid on a spectrogram
spectro(sheep,f=8000,ovlp=75,zp=16,scale=FALSE,palette=reverse.gray.colors.2)
par(new=TRUE)
fund(sheep,f=8000,fmax=300,type="p",pch=24,ann=FALSE,
  xaxs="i",yaxs="i",col="black",bg="red",threshold=6)
}

\keyword{dplot}
\keyword{ts}
