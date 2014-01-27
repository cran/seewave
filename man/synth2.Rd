\name{synth2}

\alias{synth2}

\title{Synthesis of time wave (tonal model)}

\description{
  This functions synthesizes pure tone sound based on an amplitude
  envelope and and instantaneous frequency. The function can also be
  used to modified a reference sound.
}

\usage{
synth2(env = NULL, ifreq, f, plot = FALSE, listen = FALSE, output = "matrix", ...)
}

\arguments{
  \item{env}{a numeric vector describing the amplitude envelope
  (i.e. the amplitude modulation). By default \code{NULL}, generating a
  squared enveloppe.}
  \item{ifreq}{a numeric vector describing the instantaneous frequency
  (in Hz).}
  \item{f}{a numeric vector for the sampling frequency (in Hz)}
  \item{plot}{if \code{TRUE} returns the spectrogram
    of the synthezised sound (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{spectro}} graphical parameters.}
}

\details{\code{env} and \code{ifreq} must have exactly the same
  length.\cr
  The amplitude envelope can be obtained with the Hilbert envelope
  (function \code{\link{env}}) and the instantaneous frequency can be
  obtained with the Hilbert transform (function
  \code{\link{ifreq}}). This opens a great variety of signal
  modidications as shown in the example section.
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\references{
  Beeman, K. 1998 Digital signal analysis, editing and synthesis, in
  \emph{Animal acoustic communication} edited by Hopp SL, Owren MJ, Evans CS,
  Springer, 59-103.
}

\author{Jérôme Sueur and Laurent Lellouch}

\seealso{\code{\link{synth2}},  \code{\link{noisew}}, \code{\link{pulse}}, \code{\link{echo}}}

\examples{
## You can use plot=TRUE and spectro() options
## to directly 'see' the new-built sounds
## MODIFICATION OF A REFERENCE SIGNAL
data(tico)
env.tico <- env(tico, f=22050, plot=FALSE)
ifreq.tico <- ifreq(tico, f=22050, plot=FALSE)$f[,2]
# recover the original signal
s <- synth2(env=env.tico, ifreq=ifreq.tico*1000, f=22050)
# original signal with instantaneous frequency reversed
s <- synth2(env=env.tico, ifreq=rev(ifreq.tico)*1000, f=22050)
# original signal with a +1000 Hz linear frequency shift
s <- synth2(env=env.tico, ifreq=ifreq.tico*1000+1000, f=22050)
# original signal with instantaneous frequency multiplied by 2
s <- synth2(env=env.tico, ifreq=ifreq.tico*1000*2, f=22050)
# original signal with a linear  instantaneous frequency at 2000 Hz
s <- synth2(env=env.tico, ifreq=rep(2000, times=length(tico@left)), f=22050)

## DE NOVO SYNTHESIS
# instantaneous frequency increasing by step of 500 Hz
s <- synth2(ifreq=rep(c(500,1000,1500,2000,2500,3000,3500,4000), each=2000), f=16000)
# square function of the instantaenous frequency
s <- synth2(ifreq=500+seq(-50,50, length.out=8000)^2, f=8000)
# linear increase of the amplitude envelope
s <- synth2(env=seq(0,1,length=8000), ifreq=rep(2000,8000), f=8000)
# square-root increase of the amplitude envelope
s <- synth2(env=sqrt(seq(0,1,length=8000)), ifreq=rep(2000,8000), f=8000)
# square-root increase and decrease of the amplitude envelope
s <- synth2(env=c(sqrt(seq(0,1,length=4000)), sqrt(seq(1,0,length=4000))),
    ifreq=rep(2000,8000), f=8000)
# amplitude envelope and instantaneous frequency following a normal density shape
norm <- rep(dnorm(-4000:3999, sd=1000), 2)
s <- synth2(env=norm, ifreq=500+(norm/max(norm))*1000, f=8000)
}

\keyword{datagen}
\keyword{ts}