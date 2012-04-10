\name{th}

\alias{th}

\title{Temporal entropy}

\description{Compute the entropy of a temporal envelope.}

\usage{th(env)}

\arguments{\item{env}{a data set resulting of an envelope obtained
    using \code{\link{env}}}}

\details{
  Temporal entropy is calculated according to:\cr
  \deqn{S = -\frac{\sum_{i = 1}^{N} y_i log_2(y_i)}{log_2(N)}}{% 
    S = -sum(ylogy)/log(N)}
  with:\cr
  \emph{y} = relative amplitude of the \emph{i} envelope point,\cr
  and
  \deqn{\sum_{i = 1}^{N} y_i = 1}{% 
    sum(y) = 1}\cr 
  and \emph{N} = number of envelope points.
}

\value{A single value varying between 0 and 1 is returned.
  The value has no unit.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{
  The temporal entropy of a noisy signal with many amplitude modulations
  will tend towards 1 whereas the temporal entropy of quiet signal
  will tend towards 0. \cr
  Note, however, that a sustained sound with an almost
  flat envelope will also show a very high temporal entropy. See examples.}

\seealso{\code{\link{sh}}, \code{\link{csh}}, \code{\link{H}}}

\examples{
# Temporal entropy of a cicada song
data(orni)
envorni<-env(orni,f=22050,plot=FALSE)
th(envorni)
# Smoothing the envelope might slightly change the result.
envorniS<-env(orni,f=22050,smooth=c(50,0),plot=FALSE)
th(envorniS)
# If we mute a part of the cicada song, the temporal entropy decreases
orni2<-mutew(orni,f=22050,from=0.3,to=0.55,plot=FALSE)
envorni2<-env(orni2,f=22050,plot=FALSE)
th(envorni2)
# The temporal entropy of noise tends towards 1
a<-noisew(d=1,f=8000)
enva<-env(a,f=8000,plot=FALSE)
th(enva)
# But be aware that the temporal entropy
# of a sustained sound also tends towards 1
b<-synth(f=8000,d=1,cf=2000,plot=FALSE)
envb<-env(b,f=8000,plot=FALSE)
th(envb)
}

\keyword{ts}

