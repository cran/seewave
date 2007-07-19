\name{diffenv}

\alias{diffenv}

\title{Difference between two amplitude envelopes}

\description{This function estimates the surface difference between
two amplitude envelopes.}

\usage{
diffenv(wave1, wave2, f, smooth = NULL, plot = FALSE,
lty1 = 1, lty2 = 2, col1 = 2, col2 = 4, cold = 8, xlab = "Time (s)",
ylab = "Amplitude",...)
}

\arguments{
  \item{wave1}{data or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}) describing a first time wave.}
  \item{wave2}{data or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}) describing a second time wave.}
  \item{f}{sampling frequency of \code{wave1} and \code{wave1} (in Hz).
          Does not need to be specified if \code{wave1} and/or \code{wave2} are/is (a)
          \code{\link[Sound]{Sample}} object(s).}
  \item{smooth}{smoothes both envelopes by floating average. See examples.}
  \item{plot}{logical, if \code{TRUE} plots both envelopes and their surface
      difference (by default \code{FALSE}).}
  \item{lty1}{line type of the first envelope (envelope of \code{wave1}).}
  \item{lty2}{line type of the second envelope (envelope of \code{wave2}).}
  \item{col1}{colour of the first envelope (envelope of \code{wave1}).}
  \item{col2}{colour of the second envelope (envelope of \code{wave2}).}
  \item{cold}{colour of the surface difference.}
  \item{xlab}{title of the time axis.}
  \item{ylab}{title of the amplitude axis.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
Envelopes of both waves are first transformed as probability mass functions (PMF).\cr
Envelope difference is then computed according to:\cr
\deqn{D = \frac{\sum{|env1-env2|}}{2}, with D \in [0,1].}{% 
      D = sum(abs(env1-env2))/2, with 0 < D <1.}
}

\value{
The difference is returned. This value is without unit.
When \code{plot} is \code{TRUE}, both envelopes and their difference surface are
plotted on the same graph.
}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}.}

\note{This method can be used as a relative distance estimation
between different envelopes.}

\seealso{\code{\link{oscillo}}, \code{\link{corenv}}, \code{\link{diffspec}},
\code{\link{diffwave}}}

\examples{
data(tico)
data(orni)
# selection in tico to have two waves of similar duration (length)
tico<-as.matrix(tico[1:15842,1])
diffenv(tico,orni,f=22050,plot=TRUE)
# [1] 0.7362514
# smoothing the envelope gives a better graph but changes the result
diffenv(tico,orni,f=22050,smooth=20,plot=TRUE)
# [1] 0.6871243
# smoothing can have actually dramatic effects
# so always use the same smooth settings for comparisons
a<-synth(d=0.5,f=8000,cf=1000,plot=FALSE)
b<-synth(d=0.5,f=8000,cf=1100,plot=FALSE)
diffenv(a,b,f=8000)
# [1] 0.2727784
diffenv(a,b,f=8000,smooth=50)
# [1] 0.01012826
}

\keyword{dplot}
\keyword{ts}

