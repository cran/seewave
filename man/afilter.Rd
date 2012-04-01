\name{afilter}

\alias{afilter}

\title{Amplitude filter}

\description{This function deletes all signal which amplitude
  is below a selected threshold.}

\usage{afilter(wave, f, threshold = 5, plot = TRUE,
listen = FALSE, output = "matrix", ...)}

\arguments{
\item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{threshold}{amplitude threshold (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the new oscillogram
      (by default \code{TRUE}).}
    \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return,
    either \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{
The \code{threshold} value is in \% relative to the maximal value
of \code{wave}. Signal inferior to this value is clipped. 
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

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
