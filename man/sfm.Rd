\name{sfm}

\alias{sfm}

\title{Spectral Flatness Measure}

\description{This function estimates the flatness of a frequency spectrum.}

\usage{sfm(spec)}

\arguments{\item{spec}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}}

\details{
  SFM is calculated as the ratio between the geometric mean and the
  arithmetic mean :\cr
  \deqn{F = N\times{\frac{\sqrt[N]{\prod_{i=1}^N{y_i}}}{\sum_{i=1}^N{y_i}}}}{% 
    F = N*(prod(y_i)^(1/N) / sum(y_i))}
  with:\cr
  \emph{y} = relative amplitude of the \emph{i} frequency,\cr
  and \emph{N} = number of frequencies.
}

\value{A single value varying between 0 and 1 is returned.
  The value has no unit.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{
  The SFM of a noisy signal will tend towards 1 whereas
  the SFM of a pure tone signal will tend towards 0.\cr
  See \code{\link{sh}} for another measure of signal noisiness/pureness.}

\seealso{\code{\link{sh}}, \code{\link{csh}}}

\examples{
a<-synth(f=8000,d=1,cf=2000,plot=FALSE)
speca<-spec(a,f=8000,at=0.5,plot=FALSE)
sfm(speca)
# [1] 0
b<-noisew(d=1,f=8000)
specb<-spec(b,f=8000,at=0.5,plot=FALSE)
sfm(specb)
# [1] 0.8233202
}

\keyword{ts}