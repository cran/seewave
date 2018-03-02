\name{convSPL}

\alias{convSPL}

\title{Convert sound pressure level in other units}

\description{
This function converts sound pressure level (in dB)
in sound power (Watt), intensity (Watt/m2) and pressure (Pa). By default, these
conversions are applied to air-borne sound.}

\usage{
convSPL(x, d = 1, Iref = 10^-12, pref = 2*10^-5)
}

\arguments{
  \item{x}{a numeric vector or a matrix describind SPL values (in dB).}
  \item{d}{the distance from the sound source where SPL values have been
        measured (in meter) (by default = 1m)}
  \item{Iref}{reference intensity (in Watt/m2) (by default = 10^-12)}
  \item{pref}{reference pressure (in Pa) (by default = 2*10^-5)}
}


\value{
  \code{convSPL} returns a list containing three components:
    \item{P}{data converted in sound power (in Watt).}
    \item{I}{data converted in sound intensity (in Watt/m2).}
    \item{p}{data converted in sound pressure (in Pa).}
}

\references{
Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}
                
\note{
\code{Iref} and \code{pref} correspond to a 1 kHz sound in air. 
}

\seealso{\code{\link{moredB}}, \code{\link{dBweight}}, \code{\link{attenuation}}}

\examples{
# conversion of two SPL measurements taken at 0.5 m from the source
convSPL(c(80,85),d=0.5) 
}

\keyword{math}
