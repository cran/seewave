\name{sddB}

\alias{sddB}

\title{Standard deviation of dB values}

\description{This function estimates the standard deviation of dB values}

\usage{sddB(x, level = "IL")}

\arguments{
  \item{x}{a numeric vector.}
  \item{level}{intensity level (\code{"IL"}) or sound pressure level (\code{"SPL"})}
}

\details{The standard deviation of dB values is not linear. The function
is an estimation not an exact computation which is not possible.}

\value{A numeric vector of length 1.}

\references{Wikipedia, \url{http://en.wikipedia.org/wiki/Propagation_of_uncertainty}}

\author{Jérôme Sueur}

\seealso{\code{\link{meandB}}, \code{\link{moredB}},
  \code{\link{convSPL}}, \code{\link{dBweight}}}

\examples{
sddB(c(89,90,95))
sddB(c(89,90,95), level="SPL")
}

\keyword{math}
