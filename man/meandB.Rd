\name{meandB}

\alias{meandB}

\title{Mean of dB values}

\description{This function calculates the mean of dB values}

\usage{
meandB(x, level="IL")
}

\arguments{
  \item{x}{a numeric vector or a numeric matrix.}
  \item{level}{intensity level (\code{"IL"}) or sound pressure level (\code{"SPL"})}
}

\details{The mean of dB values is not linear. See examples.}

\value{A numeric vector of length 1 is returned.}

\references{Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jerome Sueur and Zev Ross}

\seealso{\code{\link{sddB}}, \code{\link{moredB}}, \code{\link{convSPL}}, \code{\link{dBweight}}}

\examples{
meandB(c(89,90,95))
}

\keyword{math}

