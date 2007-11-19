\name{zapsilw}

\alias{zapsilw}

\title{Zap silence periods of a time wave}

\description{
This function simply delete the silence periods of a time wave. 
}

\usage{zapsilw(wave, f, threshold = 5, plot = TRUE, Sample = FALSE, ...)}

\arguments{
  \item{wave}{data describing the time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{threshold}{amplitude threshold (in \%) between silence and signal.}
  \item{plot}{logical, if \code{TRUE} plots the new oscillogram
      (by default \code{TRUE}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.      
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jérôme Sueur \email{sueur@mnhn.fr}}

\note{Use the argument \code{threshold} to set the level of silence. See
the \code{examples}.} 

\seealso{\code{\link{afilter}}, \code{\link{oscillo}}}

\examples{
data(orni)
op<-par(mfrow=c(3,1))
oscillo(orni,f=22050)
title(main = "original signal")
zapsilw(orni,f=22050,colwave="red")
title(main = "threshold level = 5")
zapsilw(orni,f=22050,threshold=1,colwave="blue")
title(main = "threshold level = 1")
par(op)
}

\keyword{ts}
