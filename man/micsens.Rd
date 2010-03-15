\name{micsens}

\alias{micsens}

\title{Microphone sensitivity and conversion}

\description{
  This function converts microphone sensitivity from mV/Pa to dB. 
}

\usage{micsens(x, sref = 1, inverse = FALSE)}

\arguments{
  \item{x}{a measured sensitivity in mV/Pa (or in dB if \code{inverse} is \code{TRUE})}
  \item{sref}{the sensitivity reference (by default equals to 1 V/Pa)}
  \item{inverse}{logical, if \code{TRUE}, the inverse conversion from dB
    to mV/Pa is computed.}
}

\details{
  The sensitivity \emph{S} in dB is calculated according to:\cr
  \deqn{S_{dB} = 20 \times log_{10}(\frac{s}{s_{ref}})}{% 
    S = 20*log(s/sref)} with \emph{s} the measured sensitivity in mv/Pa 
    and \emph{sref} the reference sensitivity (by default 1 mV/Pa).
}

\value{
  A numeric value in dB \emph{re} 1V/Pa with default settings,
  in mV/Pa if \code{inverse} is set to \code{FALSE}.
}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{convSPL}}}

\examples{
# conversion of a sensitivity of 2 mV/Pa
micsens(2)
# conversion of a sensitivity of -54 dB re 1V/Pa
micsens(-54,inverse=TRUE)
}

\keyword{math}
