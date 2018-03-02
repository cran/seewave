\name{itakura.dist}

\alias{itakura.dist}

\title{Itakuro-Saito distance}

\description{Compare two distributions (e.g. two frequency spectra) by
  computing the Itakuro-Saito distance}

\usage{itakura.dist(spec1, spec2, scale=FALSE)}

\arguments{
    \item{spec1}{any distribution, especially a spectrum obtained with \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This can be either a two-column matrix (col1 = frequency, col2 = amplitude) or a vector (amplitude).}
    \item{spec2}{any distribution, especially a spectrum obtained with
    \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This can be
    either a two-column matrix (col1 = frequency, col2 = amplitude) or a
    vector (amplitude).}
    \item{scale}{a logical, if \code{TRUE} the distance is scaled by dividing the
      distance by the length of \code{spec1} (or \code{spec2}).}
    }

\details{The Itakura-Saito (I-S) distance is a
  non-symmetric measure of the difference between two probability
  distributions. It is here adapted for frequency spectra. The distance
  is asymmetric, ie computing the I-S distance between spec1 and spec2 is
  not the same as computing it between spec2 and spec1. A symmetry can be
  obtained by calculating the mean between the two directions.\cr
  The distance is obtained following:\cr
  \deqn{D_{I-S}(spec1 \Vert spec2) = \sum{\frac{spec1}{spec2} -
  log(\frac{spec1}{spec2}) - 1}}{% 
    D(spec1 || spec2) = sum(spec1/spec2 - log(spec1/spec2) - 1)}
}

\value{The function returns a list of three items:
    \item{D1}{The I-S distance of 'spec2' with respect to 'spec1'
    (\emph{i.e.} D(spec1 || spec2))}
    \item{D2}{The I-S distance of 'spec1' with respect to 'spec2'
    (\emph{i.e.} D(spec2 || spec1))}
  \item{D}{The symmetric distance (\emph{i.e.} D = 0.5*(D1+D2))}
  If \code{scale = TRUE} the distance is divided by the length of \code{spec1}
     (or \code{spec2}).
  }

\note{The function works for both Hz and (htk-)mel scales.}
  
\author{Jerome Sueur, improved by Laurent Lellouch}

\seealso{\code{\link{kl.dist}}, \code{\link{ks.dist}}, \code{\link{logspec.dist}}, \code{\link{simspec}}, \code{\link{diffspec}}}

\examples{
# Comparison of two spectra
data(tico)
tico1 <- spec(tico, at=0.65, plot=FALSE)
tico2 <- spec(tico, at=1.1, plot=FALSE)
itakura.dist(tico1, tico2)
itakura.dist(tico1, tico2, scale=TRUE) 
}


\keyword{distribution}
\keyword{ts}
