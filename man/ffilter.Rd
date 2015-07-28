\name{ffilter}

\alias{ffilter}

\title{Frequency filter}

\description{This function filters out a selected frequency section
  of a time wave (low-pass, high-pass, low-stop, high-stop, bandpass
  or bandstop frequency filter).}

\usage{
ffilter(wave, f, from = NULL, to = NULL, bandpass = TRUE,
custom = NULL, wl = 1024, ovlp = 75, wn = "hanning", fftw = FALSE, output="matrix")
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{from}{start frequency (in Hz) where to apply the filter.}
  \item{to}{end frequency (in Hz) where to apply the filter.}
  \item{bandpass}{if \code{TRUE} a band-pass filter is applied between
    \code{from} and  \code{to}, if \code{FALSE} a band-stop filter is applied
    between \code{from} and  \code{to} (by default \code{TRUE}).}
  \item{custom}{a vector describing the frequency response of a custom filter.
    This can be manually generated or obtained with \code{\link{spec}} and \code{\link{meanspec}}.
    The legnth of the vector should be half the length of \code{wl}. See examples.}
  \item{wl}{window length for the analysis (even number of points).}
  \item{ovlp}{overlap between successive FFT windows (in \%).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes of the \code{spectro}.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"audioSample"} or \code{"ts"}.}
}

\details{
  A short-term Fourier transform is first applied to the signal
  (see \code{\link{spectro}}), then the frequency filter is applied and the new
  signal is eventually generated using the reverse of the Fourier Transform
  (\code{\link{istft}}).\cr
  There is therefore neither temporal modifications nor 
  amplitude modifications.
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur}

\seealso{\code{\link{afilter}},\code{\link{lfs}},\code{\link{fir}}}

\examples{
a<-noisew(f=8000,d=1)
# low-pass
b<-ffilter(a,f=8000,to=1500)
spectro(b,f=8000,wl=512)
# high-pass
c<-ffilter(a,f=8000,from=2500)
spectro(c,f=8000,wl=512)
# band-pass
d<-ffilter(a,f=8000,from=1000,to=2000)
spectro(d,f=8000,wl=512)
# band-stop
e<-ffilter(a,f=8000,from=1500,to=2500,bandpass=FALSE)
spectro(e,f=8000,wl=512)
# custom
myfilter1<-rep(c(rep(0,64),rep(1,64)),4)
g<-ffilter(a,f=8000,custom=myfilter1)
spectro(g,f=8000)
}

\keyword{ts}

