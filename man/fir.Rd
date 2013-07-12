\name{fir}

\alias{fir}

\title{Finite Impulse Response filter}

\description{This function is a FIR filter that filters out a selected frequency
  section of a time wave (low-pass, high-pass, low-stop, high-stop, bandpass
  or bandstop frequency filter).}

\usage{
fir(wave, f, from = NULL, to = NULL, bandpass = TRUE, custom = NULL,
wl = 512, wn = "hanning", listen = FALSE, output = "matrix")
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{from}{start frequency (in Hz) where to apply the filter.}
  \item{to}{end frequency (in Hz) where to apply the filter.}
  \item{bandpass}{if \code{TRUE} a band-pass filter is applied between
  \code{from} and  \code{to}, if not \code{NULL} a band-stop filter is applied
    between \code{from} and  \code{to} (by default \code{NULL}).}
  \item{custom}{a vector describing the frequency response of a custom filter.
    This can be manually generated or obtained with \code{\link{spec}} and \code{\link{meanspec}}.
    \code{wl} is no more required. See examples.}
  \item{wl}{window length of the impulse filter (even number of points).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
}

\details{
  This function is based on the reverse of the Fourier Transform
  (\code{\link{fft}}) and on a convolution (\code{\link{convolve}}) between the
  wave to be filtered and the impulse filter.
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\references{Stoddard, P. K. (1998). Application of filters in bioacoustics.
  \emph{In}: Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds), \emph{Animal acoustic 
    communication}. Springer, Berlin, Heidelberg,pp. 105-127.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{ffilter}}, \code{\link{lfs}}, \code{\link{afilter}}}

\examples{
a<-noisew(f=8000,d=1)
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
myfilter1<-rep(c(rep(0,32),rep(1,32)),4)
g<-fir(a,f=8000,custom=myfilter1)
spectro(g,f=8000)
# custom filter generated using spec()
data(tico)
myfilter2<-spec(tico,f=22050,at=0.7,wl=512,plot=FALSE)
b<-noisew(d=1,f=22050)
h<-fir(b,f=22050,custom=myfilter2)
spectro(h,f=22050)
}

\keyword{ts}

