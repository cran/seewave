\name{diffwave}

\alias{diffwave}

\title{Difference between two time waves}

\description{
This function estimates the difference between two waves by computing
the product between envelope surface difference and frequency surface difference.
}

\usage{
diffwave(wave1, wave2, f, wl = 512, smooth = NULL)
}

\arguments{
  \item{wave1}{data or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}) describing a first time wave.}
  \item{wave2}{data or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}) describing a second time wave.}
  \item{f}{sampling frequency of \code{wave1} and \code{wave1} (in Hz).
          Does not need to be specified if \code{wave1} and/or \code{wave2} are/is (a)
          \code{\link[Sound]{Sample}} object(s).}
  \item{wl}{window length for spectral analysis (even number of points).}  
  \item{smooth}{envelope smoothing by floating average.}
}

\details{
This function computes the product between the values obtained with
\code{\link{diffspec}} and \code{\link{diffenv}} functions. \cr
This then gives a global (time and frequency) estimation of dissimilarity.\cr
The frequency mean spectrum and the amplitude envelope needed for computing
respectively \code{\link{diffspec}} and \code{\link{diffenv}} are automatically generated.
They can be controlled through \code{wl} and  \code{smooth} arguments respectively.
See examples below and examples in \code{\link{diffspec}} and \code{\link{diffenv}}
for implications on the results.}

\value{A single value varying between 0 and 1 is returned.
The value has no unit.}

\references{Sueur, J., Pavoine, S., Hamerlynck, O. & Duvail, S., in preparation.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\note{This method can be used as a relative distance estimation
between different waves.}

\seealso{\code{\link{diffspec}}, \code{\link{diffenv}}}

\examples{
data(tico)
data(orni)
# selection in tico to have two waves of similar duration (length)
tico<-as.matrix(tico[1:15842,1])
diffwave(tico,orni,f=22050)
# [1] 0.4027229
# changing the frequency parameter (wl)
diffwave(tico,orni,f=22050,wl=1024)
# [1] 0.4335703
# changing the temporal parameter (smooth)
diffwave(tico,orni,f=22050,smooth=20)
# [1] 0.3758508
}

\keyword{ts}
