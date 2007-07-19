\name{fir}

\alias{fir}

\title{Finite Impulse Response filter}

\description{This function is a FIR filter that filters out a selected frequency
section of a time wave (low-pass, high-pass, low-stop, high-stop, bandpass
or bandstop frequency filter).}

\usage{fir(wave, f, from = FALSE, to = FALSE, bandpass = TRUE, custom = NULL,
wl = 512)}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[Sound]{Sample}} object created generated a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{from}{start frequency (in Hz) where to apply the filter.}
  \item{to}{end frequency (in Hz) where to apply the filter.}
  \item{bandpass}{if \code{TRUE} a band-pass filter is applied between
    \code{from} and  \code{to}, if \code{FALSE} a band-stop filter is applied
    between \code{from} and  \code{to} (by default \code{TRUE}).}
  \item{custom}{a vector describing the frequency response of a custom filter.
    This can be manually generated or obtained with \code{\link{spec}} and \code{\link{meanspec}}.
    \code{wl} is no more required. See examples.}
  \item{wl}{window length of the impulse filter (even number of points).}
}

\details{
This function is based on the reverse of the Fourier Transform
(\code{\link{fft}}) and on a convolution (\code{\link{convolve}}) between the
wave to be filtered and the impulse filter.
}

\value{\code{fir} returns a one-column matrix describing the new wave.}

\references{Stoddard, P. K. (1998). Application of Filters in Bioacoustics.
\emph{In}: Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds), \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg,pp. 105-127.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\seealso{\code{\link{ffilter}},\code{\link{lfs}}, \code{\link{afilter}},
\code{\link[signal]{fir1}}, \code{\link[Sound]{fir2}}}

\examples{
a<-noise(f=8000,d=1)
# low-pass
b<-fir(a,f=8000,to=1500)
spectro(b,f=8000)
# high-pass
c<-fir(a,f=8000,from=2500)
spectro(c,f=8000)
# band-pass
d<-fir(a,f=8000,from=1000,to=2000)
spectro(d,f=8000)
# band-stop
e<-fir(a,f=8000,from=1500,to=2500,bandpass=FALSE)
spectro(e,f=8000)
# custom filter manually generated
myfilter1<-rep(c(rep(0,64),rep(1,64)),4)
g<-fir(a,f=8000,custom=myfilter1)
spectro(g,f=8000)
# custom filter generated using spec()
data(orni)
myfilter2<-spec(orni,f=22050,at=0.21,wl=512,plot=FALSE)
b<-noise(d=1,f=22050)
h<-fir(b,f=8000,custom=myfilter2)
spectro(h,f=8000)
}

\keyword{ts}

