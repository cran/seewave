\name{afilter}

\alias{afilter}

\title{Amplitude filter}

\description{This function deletes all signal which amplitude
is below a selected threshold.}

\usage{afilter(wave, f, threshold = 5, plot = TRUE,
listen = FALSE, Sample = FALSE,...)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{threshold}{amplitude threshold (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the new oscillogram
      (by default \code{TRUE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.      
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{
The \code{threshold} value is in \% relative to the maximal value
of \code{wave}. Signal inferior to this value is clipped. 
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

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
