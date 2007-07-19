\name{moredB}

\alias{moredB}

\title{Addition of dB values}

\description{This functions calculates the sum of dB values}

\usage{
moredB(x)
}

\arguments{
  \item{x}{a numeric vector or a matrix.}
}

\details{
The addition of dB values is not linear. See examples. 
}

\value{
A numeric vector or a matrix is returned.
}

\references{
Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\examples{
# two sources of 60 dB give an intensity or pressure level of 63 dB
moredB(c(60,60))
# addition of three sources
moredB(c(89,90,95))
}

\keyword{math}

