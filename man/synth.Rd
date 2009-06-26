\name{synth}

\alias{synth}

\title{Synthesis of time wave}

\description{
This functions synthesize pure tone sound with amplitude modulation (am)
and/or frequency modulation (fm).
}

\usage{
synth(f, d, cf, a = 1, shape = NULL, p = 0,
am = c(0, 0), fm = c(0, 0, 0),
plot = FALSE, listen = FALSE, Sample = FALSE,...)
}

\arguments{
  \item{f}{sampling frequency (in Hz).}
  \item{d}{duration (in s).}
  \item{cf}{carrier frequency (in Hz).}
  \item{a}{amplitude (linear scale, relative when adding different waves).}
  \item{shape}{modification of the whole amplitude shape of the wave.
  See \code{details}).}
  \item{p}{phase (in radians).}
  \item{am}{a vector of length 2 describing amplitude modulation parameters,
      see \code{details}.}
  \item{fm}{a vector of length 3 describing frequency modulation parameters,
      see \code{details}.}
  \item{plot}{if \code{TRUE} returns the spectrogram
      of the synthezised sound (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
  \item{\dots}{other \code{\link{spectro}} graphical parameters.}
}
                                   
\details{
\code{shape} allows to modify the whole amplitude shape of the wave. There
are four options to be given as as character string:
(i) "incr" : linear increase  
(ii) "decr" : linear decrease  
(iii) "sine" : sinusoid-like shape
(iv) "tria" : triangular shape

\code{am} is a vector of length 2 including:\cr
(1) the amplitude modulation
depth (in \%),\cr
(2) the frequency of the amplitude modulation.\cr

\code{fm} is a vector of length 3 including:\cr 
(1) the maximum excursion of a sinusoidal frequency modulation (in Hz),\cr
(2) the frequency of a sinusoidal frequency modulation,\cr
(3) the maximum excursion of a linear frequency modulation (in Hz).
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\references{
Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jerome Sueur \email{sueur@mnhn.fr},
Thierry Aubin \email{thierry.aubin@u-psud.fr} and 
Caroline Simonis \email{csimonis@mnhn.fr}.}

\seealso{\code{\link{noise}}, \code{\link{pulse}}, \code{\link{echo}}}

\examples{
# pure tone
synth(f=22050,d=1,cf=4000,plot=TRUE)
# pure tone with sinusoid-like overall shape
synth(f=22050,d=1,cf=4000,shape="sine",plot=TRUE,osc=TRUE)
# pure tones with am
synth(f=22050,d=1,cf=4000,am=c(50,10),plot=TRUE,osc=TRUE)
# pure tone with +2000 Hz linear fm 
synth(f=22050,d=1,cf=4000,fm=c(0,0,2000),plot=TRUE)
# pure tone with sinusoidal fm
# (maximum excursion of 1000 Hz, frequency of 10 Hz)
synth(f=22050,d=1,cf=4000,fm=c(1000,10,0),plot=TRUE,wl=256,ovlp=75)
# pure tone with sinusoidal am
# (maximum excursion of 1000 Hz, frequency of 10 Hz)
# and linear fm (maximum excursion of 1000 Hz)
synth(f=22050,d=1,cf=4000,fm=c(1000,10,1000),plot=TRUE,wl=256,ovlp=75)
# the same with am
synth(f=22050,d=1,cf=4000,am=c(50,10),
    fm=c(1000,10,1000),plot=TRUE,wl=256,ovlp=75,osc=TRUE)
# the same with am and a triangular overall shape 
synth(f=22050,d=1,cf=4000,shape="tria",am=c(50,10),
    fm=c(1000,10,1000),plot=TRUE,wl=256,ovlp=75,osc=TRUE)   
# more complex sound
F1<-synth(f=22050,cf=2000,d=1,fm=c(500,5,0))
F2<-synth(f=22050,a=0.8,cf=4000,d=1,fm=c(500,5,0))
F3<-synth(f=22050,a=0.6,cf=6000,d=1,fm=c(500,5,0))
F4<-synth(f=22050,a=0.4,cf=8000,d=1,fm=c(500,5,0))
final1<-F1+F2+F3+F4
spectro(final1,f=22050,wl=512,ovlp=75,scale=FALSE)
}

\keyword{datagen}
\keyword{ts}
