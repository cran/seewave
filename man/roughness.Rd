\name{roughness}

\alias{roughness}

\title{Roughness or total curvature}

\description{This function computes the roughness or total curvature of
  a curve, i.e. of a time wave or of a spectrum}

\usage{roughness(x, std = FALSE)}

\arguments{
  \item{x}{a vector}
  \item{std}{a logical, if set to \code{TRUE} then \code{x} is standardized by its maximum.}
}

\details{
  Roughness or total curvature is the integrated squared second
  derivative :\cr
  \deqn{roughness = \int [D^2x(t)]^2 \, \mathrm{d} t}{%
      roughness = sum(second derivative^2)}.
  }

\value{A vector of length 1.}

\references{
Ramsay JO, Silverman BW (2005) \emph{Functional data analysis}. Springer, Berlin.
}

\author{Jerome Sueur}

\note{The value has not unit.}

\seealso{\code{\link{rugo}}, \code{\link{rms}}, \code{\link{sh}}, \code{\link{th}}, \code{\link{H}}.}

\examples{
data(tico)
spec <- meanspec(tico, plot=FALSE)[,2]
roughness(spec) 
}

\keyword{ts}
