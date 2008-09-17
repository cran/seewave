\name{ffilter}

\alias{ffilter}

\title{Frequency filter}

\description{This function filters out a selected frequency section
of a time wave (low-pass, high-pass, low-stop, high-stop, bandpass
or bandstop frequency filter).}

\usage{
ffilter(wave, f, from = FALSE, to = FALSE, bandpass = TRUE,
custom = NULL, wl = 512, wn = "hanning", Sample = FALSE)
}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object created generated a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{from}{start frequency (in Hz) where to apply the filter.}
  \item{to}{end frequency (in Hz) where to apply the filter.}
  \item{bandpass}{if \code{TRUE} a band-pass filter is applied between
    \code{from} and  \code{to}, if \code{FALSE} a band-stop filter is applied
    between \code{from} and  \code{to} (by default \code{TRUE}).}
  \item{custom}{a vector describing the frequency response of a custom filter.
    This can be manually generated or obtained with \code{\link{spec}} and \code{\link{meanspec}}.
    \code{wl} is no more required. See examples.}
  \item{wl}{window length for the analysis (even number of points).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}.}  
}

\details{
A short-term Fourier transform is first applied to the signal
(see \code{\link{spectro}}), then the frequency filter is applied and the new
signal is eventually generated using the reverse of the Fourier Transform
(\code{\link{fft}}).\cr
There is therefore neither temporal modifications nor 
amplitude modifications.
}

\value{A new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{afilter}},\code{\link{lfs}},\code{\link{fir}}}

\examples{
a<-noise(f=8000,d=1)
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
myfilter1<-rep(c(rep(0,32),rep(1,32)),4)
g<-fir(a,f=8000,custom=myfilter1)
spectro(g,f=8000)
}

\keyword{ts}

