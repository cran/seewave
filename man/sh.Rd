\name{sh}

\alias{sh}

\title{Spectral entropy}

\description{Compute the entropy (H) of a spectrum.}

\usage{sh(spec)}

\arguments{\item{spec}{a data set resulting of a spectral analysis obtained
      with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}}

\details{
Spectral entropy is calculated according to:\cr
\deqn{S = -\frac{sum(y[i], i==1, N)y[i]log[2]y[i]}{log(N)}}{% 
      S = -sum(ylogy)/log(N)}
with:\cr
\emph{y} = relative amplitude of the \emph{i} frequency,\cr
and
\deqn{sum(y[i], i==1, N)=1}{% 
      sum(y) = 1}\cr 
and \emph{N} = number of frequencies.
}


\value{Spectral entropy is returned.}

\references{Nunes, R. R., Almeida de, M. P. & Sleigh, J. W. 2004
Spectral entropy: a new method for anesthetic adequacy.
\emph{Revista Brasileira de Anestesiologia}, \bold{54}, 413-422.}

\author{Jérôme Sueur \email{jerome.sueur@univ-tours.fr}}

\note{
The spectral entropy of a noisy signal will tend towards 1 whereas
the spectral entropy of a pure tone signal will tend towards 0.}

\seealso{\code{\link{csh}}, \code{\link{spec}}, \code{\link{meanspec}}}

\examples{
a<-synth(f=8000,d=1,cf=2000,plot=FALSE)
speca<-spec(a,f=8000,wl=512,at=0.5,plot=FALSE)
sh(speca)
# [1] 0.2336412
b<-noise(d=1,f=8000)
specb<-spec(b,f=8000,wl=512,at=0.5,plot=FALSE)
sh(specb)
# [1] 0.9740975
}

\keyword{ts}

