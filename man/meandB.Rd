\name{meandB}

\alias{meandB}

\title{mean of dB values}

\description{This functions calculates the mean of dB values}

\usage{
meandB(x)
}

\arguments{
  \item{x}{a numeric vector or a matrix.}
}

\details{
  The mean of dB values is not linear. See examples. 
}

\value{
  A numeric vector of length 1 is returned.
}

\references{
  Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jerome Sueur \email{sueur@mnhn.fr} and Zev Ross}

\seealso{\code{moredB}, \code{sddB}, \code{\link{convSPL}}, \code{\link{dBweight}}}

\examples{
meandB(c(89,90,95))
}

\keyword{math}

