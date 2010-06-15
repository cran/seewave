\name{diffwave}

\alias{diffwave}

\title{Difference between two time waves}

\description{
  This function estimates the difference between two waves by computing
  the product between envelope surface difference and frequency surface difference.
}

\usage{
diffwave(wave1, wave2, f, wl = 512, envt = "hil",
msmooth = NULL, ksmooth = NULL)
}

\arguments{
  \item{wave1}{a first R object.}     
  \item{wave2}{a second R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}

  \item{wl}{window length for spectral analysis (even number of points).}  
  \item{envt}{the type of envelope to be used: either "abs" for absolute
    amplitude envelope or "hil" for Hilbert amplitude envelope. See \code{\link{env}}.}
  \item{msmooth}{a vector of length 2 to smooth the amplitude envelope with a 
    mean sliding window. The first component is the window length
    (in number of points). The second component is the overlap between
    successive windows (in \%). See \code{\link{env}}.}
  \item{ksmooth}{kernel smooth via \code{\link{kernel}}. See \code{\link{env}}.}
}

\details{
  This function computes the product between the values obtained with
  \code{\link{diffspec}} and \code{\link{diffenv}} functions. \cr
  This then gives a global (time and frequency) estimation of dissimilarity.\cr
  The frequency mean spectrum and the amplitude envelope needed for computing
  respectively \code{\link{diffspec}} and \code{\link{diffenv}} are automatically generated.
  They can be controlled through \code{wl}, \code{msmooth} and \code{ksmooth} arguments respectively.\cr
  See examples below and examples in \code{\link{diffspec}} and \code{\link{diffenv}}
  for implications on the results.}

\value{A single value varying between 0 and 1 is returned.
  The value has no unit.}

\references{Sueur, J., Pavoine, S., Hamerlynck, O. & Duvail, S. (2008) -
  Rapid acoustic survey for biodiversity appraisal. \emph{PLoS ONE}, 3(12): e4065.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{This method can be used as a relative distance estimation
  between different waves.}

\seealso{\code{\link{diffspec}}, \code{\link{diffenv}}}

\examples{
data(tico) ; tico <- tico@left
data(orni) ; orni <- orni@left
# selection in tico to have two waves of similar duration (length)
tico <- tico[1:length(orni)]
diffwave(tico,orni,f=22050)
# changing the frequency parameter (wl)
diffwave(tico,orni,f=22050,wl=1024)
# changing the temporal parameter (msmooth)
diffwave(tico,orni,f=22050,msmooth=c(20,0))
}

\keyword{ts}
