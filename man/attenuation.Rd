\name{attenuation}

\alias{attenuation}

\title{Generate sound intensity attenuation data}

\description{
This function generates dB data following theoretical spherical
attenuation of sound. 
}

\usage{
attenuation(lref, dref = 1, dstop, n, plot = TRUE,
xlab = "Distance (m)", ylab = "dB", type = "l", ...)
}

\arguments{
  \item{lref}{reference intensity or pressure level (in dB).}
  \item{dref}{reference distance corresponding to \code{lref} (in m.)
      (by default = 1).}
  \item{dstop}{maximal distance of propagation (in m.).}
  \item{n}{number of points generated between \code{dref} and \code{dstop}.}
  \item{plot}{logical, if \code{TRUE} plots attenuation against distance
      of propagation (by default \code{TRUE}).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
      See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}ss
}

\value{
If \code{plot} is \code{FALSE} return a numeric vector with the data generated.
}

\note{Sound attenuation in a free, unbounded medium behaves in accordance with
the inverse square law. \code{attenuation} generates data following this rule
from a reference point where sound intensity level (SIL) or 
sound pressure level (SPL) is known. Such theoretical data can be compared with
experimental data collected in a real environment.}

\references{
Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{convSPL}}, \code{\link{moredB}}}

\examples{
# theoretical attenuation up to 150 m of a 100 dB/1m sound source
attenuation(lref=100, dref=1, dstop=150, n=200)
}

\keyword{data}
