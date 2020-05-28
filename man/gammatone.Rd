\name{gammatone}

\alias{gammatone}

\title{Gammatone filter}

\description{Generate gammatone filter in the time domain (impulse response).}

\usage{gammatone(f, d, cfreq, n = 4, a = 1, p = 0, output = "matrix")}

\arguments{
  \item{f}{sampling frequency (in Hz).}
  \item{d}{duration (in s).}
  \item{cfreq}{center frequency (in Hz).}
  \item{n}{filter order (no unit).}
  \item{a}{amplitude (linear scale, no unit).}
  \item{p}{initial phase (in radians).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
}

\details{The gammatone function in the time domain (impulse response) is
  obtained with:\cr
  \deqn{g(t) = a \times t^{n-1} \times e^{-2\pi \beta t} \times \cos(2 \pi cf t + \phi)}{% 
    g = a * t^{n-1} * e^{-2*pi*b*t} * cos(2*pi*f*t + p)}
    with \emph{a} the amplitude, \emph{t} time, \emph{n} the filter order, \emph{cf} the center frequency, \eqn{\phi}{p} the initial phase.\cr
    The parameter \eqn{\beta}{b} is the equivalent rectangular
  bandwidth (ERB) bandwidth which varies according to the center
  frequency \eqn{cf}{cf} following:\cr
  \deqn{\beta = 24.7 \times (4.37 \times \frac{cf}{1000} + 1)}{%
  24.7*(4.37*cfreq/1000+1)
  } 
}

\value{A wave is returned. The class of the returned object is set with the argument \code{output}.}

\references{
Holdsworth J, Nimmo-Smith I, Patterson R, Rice P (1988) Implementing a gammatone filter bank. Annex C of the SVOS Final Report: Part A: The Auditory Filterbank, 1, 1-5.
}

\author{Jerome Sueur}

\note{Use the FFT based function, as \code{\link{spec}} or
  \code{\link{meanspec}}, to get the filter in the frequency domain. See examples.}

\seealso{\code{\link{melfilterbank}}}

\examples{
## gammatone filter in the time domain (impulse response)
f <- 44100
d <- 0.05
res <- gammatone(f=f, d=d, cfreq=440, n=4)
## time display
oscillo(res, f=f)
## frequency display
spec(res, f=f)
## generate and plot a bank of 32 filters from 500 to 10000 Hz
n <- 32
cfreq <- round(seq(500, 10000, length.out=n))
res <- matrix(NA, nrow=f*d/2, ncol=n)
for(i in 1:n){
    res[,i] <- spec(gammatone(f=f, d=d, cfreq=cfreq[i]), f=f, dB="max0", plot=FALSE)[,2]
}
x <- seq(0,f/2,length.out=nrow(res))/1000
plot(x=x, y=res[,1],
     xlim=c(0,14), ylim=c(-60,0),
     type="l", col=2, las=1, 
     xlab="Frequency (kHz)", ylab="Relative amplitude (dB)")
for(i in 2:n) lines(x, res[,i], col=2)
## use the frequency domain to filter a white noise input
## here around the center frequency 2000 Hz
res <- gammatone(f=f, d=d, cfreq=2000, n=4)
gspec <- spec(res, f=f, plot=FALSE)[,2]
nw <- noisew(f=44100, d=1)
nwfilt <- fir(nw, f=44100, wl=length(gspec)*2, custom=gspec) 
spectro(nwfilt, f=f)
}

\keyword{ts}
\keyword{math}
