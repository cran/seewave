\name{afilter}

\alias{afilter}

\title{Amplitude filter}

\description{This function deletes all signal which amplitude is below a selected threshold.}

\usage{afilter(wave, f, threshold = 5, plot = TRUE, ...)}

\arguments{
  \item{wave}{data describing the time wave
  or a \code{\link[Sound]{Sample}} object created loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).}
  \item{threshold}{amplitude threshold (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the new oscillogram
      (by default \code{TRUE}).}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{
The \code{threshold} value is in \% relative to the maximal value
of \code{wave}. Signal inferior to this value is clipped. 
}

\value{If \code{plot} is \code{FALSE}, \code{afilter} returns a one column
matrix with the new wave.}

\author{Jérôme Sueur \email{jerome.sueur@univ-tours.fr}}

\note{This function is used as an argument (\code{threshold}) in the 
following functions: \code{\link{autoc}}, \code{\link{csh}},
\code{\link{dfreq}}, \code{\link{timer}} and \code{\link{zc}}.   
}

\seealso{\code{\link{ffilter}}, \code{\link{oscillo}}}

\examples{
data(orni)
op<-par(mfrow=c(2,1))
afilter(orni,f=22050)
title(main = "threshold level = 5")
afilter(orni,f=22050,threshold=0.5,colwave="blue")
title(main = "threshold level = 0.5")
par(op)
}

\keyword{ts}
