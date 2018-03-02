\name{moredB}

\alias{moredB}

\title{Addition of dB values}

\description{This functions calculates the sum of dB values}

\usage{
moredB(x, level="IL")
}

\arguments{
  \item{x}{a numeric vector or numeric matrix.}
  \item{level}{intensity level (\code{"IL"}) or sound pressure level (\code{"SPL"})}
}

\details{The addition of dB values is not linear. See examples.}

\value{A numeric vector of length 1.}

\references{Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jerome Sueur}

\seealso{\code{meandB}, \code{sddB}, \code{\link{convSPL}}, \code{\link{dBweight}}}

\examples{
# two sources of 60 dB give an intensity level of 63 dB
moredB(c(60,60))
# addition of three sources
moredB(c(89,90,95))
}

\keyword{math}

