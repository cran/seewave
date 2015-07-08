\name{kl.dist}

\alias{kl.dist}

\title{Kullback-Leibler distance}

\description{Compare two distributions (e.g. two frequency spectra) by
  computing the Kullback-Leibler distance}

\usage{kl.dist(spec1, spec2, base = 2)}

\arguments{
  \item{spec1}{any distribution, especially a spectrum obtained with \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This can be either a two-column matrix (col1 = frequency, col2 = amplitude) or a vector (amplitude).}
  \item{spec2}{any distribution, especially a spectrum obtained with
    \code{\link{spec}} or \code{\link{meanspec}} (not in dB). This can be
    either a two-column matrix (col1 = frequency, col2 = amplitude) or a
    vector (amplitude).} 
  \item{base}{the logarithm base used to compute the distance. See \code{\link{log}}.}
}

\details{The Kullback-Leibler distance or relative entropy is a
  non-symmetric measure of the difference between two probability
  distributions. It is here adapted for frequency spectra. The distance
  is asymmetric, ie computing the K-L distance between spec1 and spec2 is
  not the same as computing it between spec2 and spec1. A symmetry can be
  obtained by calculating the mean between the two directions.\cr
  The distance is obtained following:\cr
  \deqn{D_{K-L}(spec1 \Vert spec2) = \sum{spec1 \times log(\frac{spec1}{spec2})}}{% 
    D(spec1 || spec2) = sum(spec1 * log(spec1 / spec2))}
}
  
\note{
    The base of the logarithm can be changed using the argument
    \code{base}. When sets to base 2, the information is measured in units of
    bits. When sets to base \emph{e}, the information is measured in
    nats.\cr
    The function works for both Hz and (htk-)mel scales.
}

\value{The function returns a list of three items:
    \item{D1}{The K-L distance of 'spec2' with respect to 'spec1'
    (\emph{i.e.} D(spec1 || spec2))}
    \item{D2}{The K-L distance of 'spec1' with respect to 'spec2'
    (\emph{i.e.} D(spec2 || spec1))}
    \item{D}{The symmetric K-L distance (\emph{i.e.} D = 0.5*(D1+D2))}
  }

\references{Kullback, S., Leibler, R.A. (1951). On information and sufficiency. \emph{Annals of Mathematical Statistics}, 22: 79-86}

\author{Jerome Sueur, improved by Laurent Lellouch}

\seealso{\code{\link{ks.dist}}, \code{\link{logspec.dist}}, \code{\link{simspec}}, \code{\link{diffspec}}}

\examples{
# Comparison of two spectra
data(tico)
tico1 <- spec(tico, at=0.65, plot=FALSE)
tico2 <- spec(tico, at=1.1, plot=FALSE)
kl.dist(tico1, tico2)    # log2 (binary logarithm)
kl.dist(tico1, tico2, base=exp(1))  # ln (natural logarithm)
}

\keyword{distribution}
\keyword{ts}