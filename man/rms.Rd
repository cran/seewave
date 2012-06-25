\name{rms}

\Rdversion{1.1}

\alias{rms}

\title{
  Root Mean Square
}

\description{This function computes the root mean square or quadratic mean.
}

\usage{
rms(x, ...)
}

\arguments{
  \item{x}{an R object}
  \item{\dots}{further arguments passed to mean}
}

\details{
  The Root Mean Square or quadratic mean is computed according to:\cr
  \deqn{RMS = \sqrt{\frac{1}{n} \times \sum_{i=1}^{N} x{_i}{^2}}}{% 
    RMS = sqrt(mean(x^2))}
}

\value{A numeric vector of length 1}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{mean}}}

\examples{
# simple rms
rms(1:10)
# rms of a normalized envelope
data(sheep)
env <- env(sheep, f=8000)
rms(env)
}

\keyword{ts}
