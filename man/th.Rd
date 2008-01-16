\name{th}

\alias{th}

\title{Temporal entropy}

\description{Compute the entropy of a temporal envelope.}

\usage{th(env)}

\arguments{\item{env}{a data set resulting of an envelope obtained
      using \code{\link{oscillo}} with the argument env = TRUE.}}

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

\author{Jérôme Sueur \email{sueur@mnhn.fr}}

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
envorni<-oscillo(orni,f=22050,env=TRUE,plot=FALSE)
th(envorni)
# [1] 0.8853535
# Smoothing the envelope might slightly change the result.
envorniS<-oscillo(orni,f=22050,env=TRUE,smooth=50,plot=FALSE)
th(envorniS)
# [1] 0.8515566
# If we mute a part of the cicada song, the temporal entropy decreases
orni2<-mutew(orni,f=22050,from=0.3,to=0.55,plot=FALSE)
envorni2<-oscillo(orni2,f=22050,env=TRUE,plot=FALSE)
th(envorni2)
# [1] 0.837402
# The temporal entropy of noise tends towards 1
a<-noise(d=1,f=8000)
enva<-oscillo(a,f=8000,env=TRUE,plot=FALSE)
th(enva)
# [1] 0.9683144
# But be aware that the temporal entropy of a sustained sound also tends towards 1
b<-synth(f=8000,d=1,cf=2000,plot=FALSE)
envb<-oscillo(b,f=8000,env=TRUE,plot=FALSE)
th(envb)
# [1] 0.983882
}

\keyword{ts}

