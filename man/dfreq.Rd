\name{dfreq}

\alias{dfreq}

\title{Dominant frequency of a time wave}

\description{
  This function gives the dominant frequency (i. e. the frequency of highest amplitude) of a time wave.}

\usage{
dfreq(wave, f, wl = 512, wn = "hanning", ovlp = 0, fftw=  FALSE, at =
NULL, tlim = NULL, threshold = NULL, bandpass = NULL, clip = NULL,
plot = TRUE, xlab = "Times (s)", ylab = "Frequency (kHz)",
ylim = c(0, f/2000), ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{length of the window for the analysis
    (even number of points, by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{ovlp}{overlap between two successive analysis windows (in \% ).}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes of the \code{spectro}.}
  \item{at}{time position where the dominant frequency has to be
    computed (in s.).}
  \item{tlim}{modifications of the time X-axis limits.}
  \item{threshold}{amplitude threshold for signal detection (in \% ).}
  \item{bandpass}{a numeric vector of length two, giving the lower and
    upper limits of a frequency bandpass filter (in Hz).}
  \item{clip}{a numeric value to select dominant frequency values according to
    their amplitude in reference to a maximal value of 1 for
    the whole signal (has to be >0 & < 1).}
  \item{plot}{logical, if \code{TRUE} plots the dominant frequency against time
    (by default \code{TRUE}).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{ylim}{the range of y values.}
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
f <- 22050
# default
dfreq(tico,f)
# using the amplitude threshold and changing the graphical output
dfreq(tico, f, ovlp=50,threshold=5, type="l", col=2)
# using 'at' argument for specific positions along the time axis
dfreq(tico, f, at=c(0.25, 0.75, 1.2, 1.6))
dfreq(tico, f, at=seq(0.5, 1.4, by=0.005), threshold=5)
# a specific number of measures on a single note
dfreq(tico, f, at=seq(0.5, 0.9, len=100), threshold=5, xlim=c(0.5,0.9))
# overlap on spectrogram
# and use of 'clip' argument to better track the dominant frequency
# in noisy conditions
op <- par()
ticon <- tico@left/max(tico@left) + noisew(d=length(tico@left)/f, f)
spectro(ticon, f)
res <- dfreq(ticon, f, clip=0.3, plot=FALSE)
points(res, col=2, pch =13)
par(op)
}

\keyword{dplot}
\keyword{ts}
