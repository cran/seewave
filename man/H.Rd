\name{H}

\alias{H}

\title{Total entropy}

\description{
This function estimates the total entropy of a time wave.
}

\usage{H(wave, f, wl = 512, envt="hil", msmooth = NULL, ksmooth = NULL)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{window length for spectral entropy analysis
    (even number of points). See \code{\link{sh}}.}
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
\code{\link{sh}} and \code{\link{th}} functions. \cr
This then gives a global (time and frequency) estimation of signal entropy.\cr
The frequency mean spectrum and the amplitude envelope needed for computing
respectively \code{\link{sh}} and \code{\link{th}} are automatically generated.
They can be controlled through \code{wl} and  \code{smooth} arguments respectively.
See examples below and examples in \code{\link{sh}} and \code{\link{th}}
for implications on the results. 
}

\value{A single value varying between 0 and 1 is returned. The value has no unit.}

\references{Sueur, J., Pavoine, S., Hamerlynck, O. & Duvail, S. (2008) - Rapid acoustic survey for biodiversity appraisal. \emph{PLoS ONE}, 3(12): e4065.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{
The entropy of a noisy signal will tend towards 1 whereas
the entropy of a pure tone signal will tend towards 0.}

\seealso{\code{\link{sh}}, \code{\link{th}}, \code{\link{csh}}}

\examples{
data(orni)
H(orni,f=22050)
# changing the spectral parameter (wl)
H(orni,f=22050,wl=1024)
# changing the temporal parameter (msmooth)
H(orni,f=22050,msmooth=c(20,0))
}

\keyword{ts}
