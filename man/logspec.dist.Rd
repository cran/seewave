\name{logspec.dist}

\alias{logspec.dist}

\title{Log-spectral distance}

\description{Compare two distributions (e.g. two frequency spectra) by
  computing the log-spectral distance}

\usage{logspec.dist(spec1, spec2)}

\arguments{
  \item{spec1}{any distribution, especially a spectrum obtained with \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This can be either a two-column matrix (col1 = frequency, col2 = amplitude) or a vector (amplitude).}
  \item{spec2}{any distribution, especially a spectrum obtained with
    \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This can be
    either a two-column matrix (col1 = frequency, col2 = amplitude) or a
    vector (amplitude).} 
}


\details{The distance is computed according to:\cr
  \deqn{D_{LS}(spec1 \Vert spec2) = D_{LS}(spec2 \Vert spec1) = \sqrt{\sum{10
  \times log_{10}(\frac{spec1}{spec2})^{2}}}}{% 
    D(spec1 || spec2) = D(spec2 || spec1) = sqrt(sum((10*log10(spec1/spec2))^2))}
}

\value{A numeric vector of length 1 returning the D distance.}

\author{Jerome Sueur}

\note{The distance is symmetric.}

\seealso{\code{\link{ks.dist}}, \code{\link{kl.dist}},
  \code{\link{itakura.dist}}, \code{\link{simspec}}, \code{\link{diffspec}}}

\examples{
# Comparison of two spectra
data(tico)
tico1 <- spec(tico, at=0.65, plot=FALSE)
tico2 <- spec(tico, at=1.1, plot=FALSE)
logspec.dist(tico1, tico2)
}

\keyword{distribution}
\keyword{ts}
