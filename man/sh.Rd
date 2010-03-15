\name{sh}

\alias{sh}

\title{Spectral entropy}

\description{This function computes the entropy of a frequency spectrum.}

\usage{sh(spec)}

\arguments{\item{spec}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}}

\details{
  Spectral entropy is calculated according to:\cr
  \deqn{S = -\frac{\sum_{i = 1}^{N} y_i log_2(y_i)}{log_2(N)}}{% 
    S = -sum(ylogy)/log(N)}
  with:\cr
  \emph{y} = relative amplitude of the \emph{i} frequency,\cr
  and
  \deqn{\sum_{i = 1}^{N} y_i = 1}{% 
    sum(y) = 1}\cr 
  and \emph{N} = number of frequencies.
}


\value{A single value varying between 0 and 1 is returned.
  The value has no unit.}

\references{Nunes RR, Almeida de MP, Sleigh JW (2004) -
  Spectral entropy: a new method for anesthetic adequacy.
  \emph{Revista Brasileira de Anestesiologia}, \bold{54}, 413-422.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{
  The spectral entropy of a noisy signal will tend towards 1 whereas
  the spectral entropy of a pure tone signal will tend towards 0.}

\seealso{\code{\link{csh}},\code{\link{th}}, \code{\link{H}}, \code{\link{sfm}}}

\examples{
a<-synth(f=8000,d=1,cf=2000,plot=FALSE)
speca<-spec(a,f=8000,at=0.5,plot=FALSE)
sh(speca)
# [1] 0.2336412
b<-noisew(d=1,f=8000)
specb<-spec(b,f=8000,at=0.5,plot=FALSE)
sh(specb)
# close to 1
}

\keyword{ts}

