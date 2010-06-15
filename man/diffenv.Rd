\name{diffenv}

\alias{diffenv}

\title{Difference between two amplitude envelopes}

\description{This function estimates the surface difference between
  two amplitude envelopes.}

\usage{
diffenv(wave1, wave2, f, envt = "hil", msmooth = NULL, ksmooth = NULL,
plot = FALSE, lty1 = 1, lty2 = 2, col1 = 2, col2 = 4, cold = 8,
xlab = "Time (s)", ylab = "Amplitude", ylim = NULL, legend = TRUE, ...)
}

\arguments{
  \item{wave1}{a first R object.}     
  \item{wave2}{a second R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{envt}{the type of envelope to be used: either "abs" for absolute
    amplitude envelope or "hil" for Hilbert amplitude envelope. See \code{\link{env}}.}
  \item{msmooth}{a vector of length 2 to smooth the amplitude envelope with a 
    mean sliding window. The first component is the window length
    (in number of points). The second component is the overlap between
    successive windows (in \%). See \code{\link{env}}.}
  \item{ksmooth}{kernel smooth via \code{\link{kernel}}. See \code{\link{env}}.}
  \item{plot}{logical, if \code{TRUE} plots both envelopes and their surface
    difference (by default \code{FALSE}).}
  \item{lty1}{line type of the first envelope (envelope of \code{wave1}).}
  \item{lty2}{line type of the second envelope (envelope of \code{wave2}).}
  \item{col1}{colour of the first envelope (envelope of \code{wave1}).}
  \item{col2}{colour of the second envelope (envelope of \code{wave2}).}
  \item{cold}{colour of the surface difference.}
  \item{xlab}{title of the time axis.}
  \item{ylab}{title of the amplitude axis.}
  \item{ylim}{range of amplitude axis.}
  \item{legend}{logical, if \code{TRUE} adds a legend to the plot.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  Envelopes of both waves are first transformed as probability mass functions (PMF).\cr
  Envelope difference is then computed according to:\cr
  \deqn{D = \frac{\sum{|env1-env2|}}{2}, with D \in [0,1].}{% 
    D = sum(abs(env1-env2))/2, with 0 < D <1.}
}

\value{
  The difference is returned. This value is without unit.
  When \code{plot} is \code{TRUE}, both envelopes and their difference surface are
  plotted on the same graph.
}

\references{Sueur, J., Pavoine, S., Hamerlynck, O. & Duvail, S. (2008) -
  Rapid acoustic survey for biodiversity appraisal. \emph{PLoS ONE},
  3(12): e4065.}

\author{Jerome Sueur \email{sueur@mnhn.fr}.}

\note{This method can be used as a relative distance estimation
  between different envelopes.}

\seealso{\code{\link{env}}, \code{\link{corenv}}, \code{\link{diffspec}},
  \code{\link{diffwave}}}

\examples{
data(tico) ; tico <- tico@left
data(orni) ; orni <- orni@left
# selection in tico of two waves with similar duration
tico2<-tico[1:length(orni)]
diffenv(tico2,orni,f=22050,plot=TRUE)
# smoothing the envelope gives a better graph but slightly changes the result
diffenv(tico2,orni,f=22050,msmooth=c(20,0),plot=TRUE)
}

\keyword{dplot}
\keyword{ts}