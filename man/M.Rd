\name{M}

\alias{M}

\title{Median of the amplitude envelope}

\description{This function computes an acoustic index based on the median of the amplitude envelope.}

\usage{M(wave, f, envt = "hil", plot = FALSE, ...)}

\arguments{
  \item{wave}{an \code{R} object.}
  \item{f}{sampling frequency of wave (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{envt}{the type of envelope to be used: either \code{"abs"} for absolute amplitude envelope or \code{"hil"} for Hilbert amplitude envelope. See \code{\link{env}}.}
  \item{plot}{logical, if TRUE returns a plot of the amplitude envelope of wave (by default \code{FALSE}).}
  \item{\dots}{other env parameters, in particular smoothing parameters. See \code{\link{env}}.}
}

\details{This amplitude index M is computed according to:
  \deqn{M = \bar{A}(t) \times 2^{1-depth}}{%
M = median(A(t)) * 2^{1-depth}
    }
    with
    \deqn{0 \leq M \leq 1}{%
0 <= M <= 1
    }
where A(t) is the amplitude envelope and depth is the signal
  digitization depth in number of bits.
}

\value{A numeric vector of length 1 between 0 and 1, without unit.}

\references{Depraetere M, Pavoine S, Jiguet F, Gasc A, Duvail S, Sueur J (2012) Monitoring animal diversity
using acoustic indices: implementation in a temperate woodland. \emph{Ecological Indicators}, \bold{13}, 46-54.}

\author{Jerome Sueur and Marion Depraetere}

\seealso{\code{\link{env}}, \code{\link{AR}}}

\examples{
data(tico)
M(tico)
# smoothing the amplitude may change slightly the result
M(tico, msmooth=c(500,50), plot=TRUE)
}

\keyword{ts}
\keyword{index}
