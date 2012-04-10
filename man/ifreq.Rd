\name{ifreq}

\alias{ifreq}

\title{Instantaneous frequency}

\description{
  This function returns the instantaneous frequency (and/or phase) of a time wave
  through the computation of the analytic signal (Hilbert transform).
}

\usage{
ifreq(wave, f, phase = FALSE, threshold = NULL,
plot = TRUE, xlab = "Time (s)", ylab = NULL,
ylim = NULL, type = "l", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{phase}{if \code{TRUE} and  \code{plot} is also \code{TRUE}
    plots the instantaneous phase instead of the instantaneous frequency.}
  \item{threshold}{amplitude threshold for signal detection (in \% ).}
  \item{plot}{logical, if \code{TRUE} plots the instantaneous frequency or phase
    against time (by default \code{TRUE}).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{ylim}{the range of y values.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}  
  \item{\dots}{other \code{\link{plot}} graphical parameters.} 
}

\details{
  The instantaneous phase is the argument of the 
  analytic signal obtained throught the Hilbert transform.\cr
  The instantaneous phase is then unwrapped and derived against time to
  get the instantaneous frequency.\cr
  There may be some edge effects at both start and end of the time wave. 
}

\value{If \code{plot} is \code{FALSE}, \code{ifreq}
  returns a list of two components:
  \item{p}{a two-column matrix, the first column corresponding to time in seconds
    (\emph{x}-axis) and the second column corresponding to wrapped
    instantaneous phase in radians (\emph{y}-axis).}
  \item{f}{a two-column matrix, the first column corresponding to time in seconds
    (\emph{x}-axis) and the second column corresponding to instantaneous
    frequency in kHz (\emph{y}-axis).}
}

\references{Mbu Nyamsi, R. G., Aubin, T. & Bremond, J. C. 1994
  On the extraction of some time dependent parameters of an acoustic signal
  by means of the analytic signal concept.
  Its application to animal sound study. \emph{Bioacoustics}, 5: 187-203.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{This function is based on the analytic signal obtained with the 
  Hilbert transform (see \code{\link{hilbert}}).\cr
  The function requires the package \pkg{signal}.\cr
  The matrix describing the instantaneous phase has one more row than the
  one describing the instantaneous frequency.
}

\seealso{\code{\link{hilbert}}, \code{\link{zc}}}

\examples{
# generate a sound with sine and linear frequency modulations
a<-synth(d=1, f=8000, cf=1500, fm=c(200,10,1000))
# plot on a single graphical device the instantaneous frequency and phase
op<-par(mfrow=c(2,1))
ifreq(a,f=8000,main="Instantaneous frequency")
ifreq(a,f=8000,phase=TRUE,main="Instantaneous phase")
par(op)
}

\keyword{ts}
\keyword{dplot}