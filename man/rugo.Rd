\name{rugo}

\alias{rugo}

\title{Rugosity of a time wave}

\description{This function computes the rugosity of a time wave or time series}

\usage{rugo(x, ...)}

\arguments{
  \item{x}{a vector}
  \item{\dots}{other \code{\link{mean}} parameters.}
}

\details{
The formula has been slightly modified from Mezquida & Martinez (2009:
826) to fit with the classical definition of the root-mean-square
(see \code{\link{rms}}).\cr
The rugosity is then computed as following:
  \deqn{rugo = \sqrt{\sum_{i=1}^{n-1} \frac{(x_{i+1}-x_{i})^2}{n}}}{%
      rugo = sqrt(mean((x[i+1]-x[i])^2))}
    for a vector \code{x} of length \emph{n}.
  }

\value{A vector of length 1.}

\references{Mezquida DA, Martinez JL (2009) -
  Platform for bee-hives monitoring based on sound analysis.
  A perpetual warehouse for swarm's daily acitvity.
  \emph{Spanish Journal of Agricultural Research} \bold{7}, 824-828.
}

\author{Jerome Sueur}

\note{The rugosity of a noisy signal will tend to be higher than that of
  a pure tone signal, all other things being equal.\cr
}

\seealso{\code{\link{roughness}}, \code{\link{rms}}, \code{\link{sh}}, \code{\link{th}}, \code{\link{H}}.}

\examples{
data(tico) ; tico <-tico@left
# rugosity of the original recording normalised
rugo(tico/max(tico))
# synthesis of white noise with the same duration as tico
noise <- noisew(d=length(tico)/22050, f=22050)
# tico is normalised to get similar amplitude with the noise
tico.norm <- tico/max(tico)
# addition of noise to tico
tico.noisy <- tico.norm + 0.5*noise
# new rugosity (higher) on normalised signal
rugo(tico.noisy/max(tico.noisy))
}

\keyword{ts}
