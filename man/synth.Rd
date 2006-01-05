\name{synth}

\alias{synth}

\title{Synthesis of time wave}

\description{
This functions synthesize pure tone sound with amplitude modulation (am)
and/or frequency modulation (fm).
}

\usage{
synth(f, d, cf, a = 1, p = 0, am = c(0, 0), fm = c(0, 0, 0),
plot = TRUE, wl = 512, ovlp = 50, ...)
}

\arguments{
  \item{f}{sampling frequency (in Hz).}
  \item{d}{duration (in s).}
  \item{cf}{carrier frequency (in Hz).}
  \item{a}{amplitude (linear scale, relative when adding different waves).}
  \item{p}{phase (in radians).}
  \item{am}{a vector of length 2 describing amplitude modulation parameters,
      see \code{details}.}
  \item{fm}{a vector of length 3 describing frequency modulation parameters,
      see \code{details}.}
  \item{plot}{if \code{TRUE} returns the spectrogram
      of the synthezised sound (by default \code{TRUE}).}
  \item{wl}{if \code{plot} is \code{TRUE},
      length of the window for the spectrographic analysis (even number of points).}
  \item{ovlp}{if \code{plot} is \code{TRUE}, overlap between two successive
      windows of the spectrographic analysis (in \%).}
  \item{\dots}{other \code{\link{spectro}} graphical parameters.}
}
                                   
\details{
\code{am} is a vector of length 2 including: (1) the amplitude modulation
depth (in \%), (2) the frequency of the amplitude modulation.

\code{fm} is a vector of length 3 including:\cr 
(1) the maximum excursion of a sinusoidal frequency modulation (in Hz),\cr
(2) the frequency of a sinusoidal frequency modulation,\cr
(3) the maximum excursion of a linear frequency modulation (in Hz).
}

\value{
If \code{plot} is \code{FALSE}, \code{synth} returns a one-column matrix
describing the sound synthesized.
}

\references{
Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jérôme Sueur \email{jerome.sueur@univ-tours.fr},
Thierry Aubin \email{thierry.aubin@ibaic.u-psud.fr} and 
Caroline Simonis-Sueur \email{csimonis@mnhn.fr}.}

\seealso{\code{\link{spectro}}}

\examples{
# pure tone
synth(f=22050,d=1,cf=4000)
# pure tones with am
synth(f=22050,d=1,cf=4000,am=c(50,10),osc=TRUE)
# pure tone with +2000 Hz linear fm 
synth(f=22050,d=1,cf=4000,fm=c(0,0,2000))
# pure tone with sinusoidal fm
# (maximum excursion of 1000 Hz, frequency of 10 Hz)
synth(f=22050,d=1,cf=4000,fm=c(1000,10,0),wl=256,ovlp=75)
# pure tone with sinusoidal am
# (maximum excursion of 1000 Hz, frequency of 10 Hz)
# and linear fm (maximum excurion of 1000 Hz
synth(f=22050,d=1,cf=4000,fm=c(1000,10,1000),wl=256,ovlp=75)
# the same with am
synth(f=22050,d=1,cf=4000,am=c(50,10),
    fm=c(1000,10,1000),wl=256,ovlp=75,osc=TRUE)
# more complex sound
F1<-synth(f=22050,cf=2000,d=1,fm=c(500,5,0),plot=FALSE)
F2<-synth(f=22050,a=0.8,cf=4000,d=1,fm=c(500,5,0),plot=FALSE)
F3<-synth(f=22050,a=0.6,cf=6000,d=1,fm=c(500,5,0),plot=FALSE)
F4<-synth(f=22050,a=0.4,cf=8000,d=1,fm=c(500,5,0),plot=FALSE)
final1<-F1+F2+F3+F4
spectro(final1,f=22050,wl=512,ovlp=75,scale=FALSE)
# insert 0.25 s silence
final1.1<-as.matrix(final1[1:(nrow(final1)/2),])
final1.2<-as.matrix(rep(0,0.25*22050))
final1.3<-as.matrix(final1[(nrow(final1)/2):nrow(final1),])
final2<-rbind(final1.1,final1.2,final1.3)
spectro(final2,f=22050,wl=512,ovlp=75,scale=FALSE)
}

\keyword{datagen}
\keyword{ts}
