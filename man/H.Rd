\name{H}

\alias{H}

\title{Total entropy}

\description{
This functions estimates the difference between two time waves by computing their
.
}

\usage{H(wave, f, wl = 512, smooth = NULL)}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{wl}{window length for spectral entropy analysis
    (even number of points). See \code{\link{sh}}.}
  \item{smooth}{envelope smoothing by floating average. See \code{\link{th}}.}
}

\details{
This function computes the product between the values obtained with
\code{\link{sh}} and \code{\link{th}} functions. \cr
This then gives a global (time and frequency) estimation of signal entropy.\cr
The frequency mean spectrum and the amplitude envelope needed for computing
respectively \code{\link{sh}} and \code{\link{th}} are automatically generated.
They can be controlled through \code{wl} and  \code{smooth} arguments respectively.
See examples below and examples in \code{\link{sh}} and \code{\link{th}}
for implications on the results. 
}

\value{A single value varying between 0 and 1 is returned. The value has no unit.}

\references{Sueur, J., Pavoine, S., Hamerlynck, O. & Duvail, S., in preparation.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\note{
The  entropy of a noisy signal will tend towards 1 whereas
the  entropy of a pure tone signal will tend towards 0.}

\seealso{\code{\link{sh}}, \code{\link{th}}, \code{\link{csh}}}

\examples{
data(orni)
H(orni,f=22050)
# [1] 0.766589
# changing the spectral parameter (wl)
H(orni,f=22050,wl=1024)
# [1] 0.7785917
# changing the temporal parameter (smooth)
H(orni,f=22050,smooth=20)
# [1] 0.7491108
}

\keyword{ts}
