\name{pastew}

\alias{pastew}

\title{Paste a time wave to another one}

\description{
This function pastes a first time wave to a second one. The time wave to be 
pasted, the time wave to be completed and the resulting time wave can be 
displayed in a three-frame oscillographic plot.
}

\usage{
pastew(wave1, wave2, f, at = FALSE, plot = FALSE, marks = TRUE, ...)
}

\arguments{
  \item{wave1}{data or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound})
  describing the time wave to be pasted.}
  \item{wave2}{data or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound})
  describing the time wave into which \code{wave1} will be pasted.}
  \item{f}{sampling frequency of \code{wave1} and \code{wave1} (in Hz).
          Does not need to be specified if \code{wave1} and/or \code{wave2} are/is (a)
          \code{\link[Sound]{Sample}} object(s).}
  \item{at}{\code{wave2} position in seconds
    where \code{wave1} will be pasted into (by default at the end of \code{wave2}).}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of 
  \code{wave1}, \code{wave2} and \code{wave1} + \code{wave2} (by default \code{FALSE}).}
  \item{marks}{logical, if \code{TRUE} shows where \code{wave1} has been pasted (by default \code{TRUE}).}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
  }
}

\details{If \code{plot} is \code{TRUE} returns a two-frame plot with three waves:\cr 
(1) the wave to be pasted (\code{wave1}),\cr
(2) the wave to be completed (\code{wave2}),\cr
(3) the resulting wave.}
 
\value{A one-column matrix describing a wave with the same sampling
frequency as original data.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\seealso{\code{\link{cutw}}, \code{\link{deletew}}, \code{\link{mute}},
\code{\link{oscillo}}, \code{\link{timer}}.}

\examples{
data(tico)
# double a data set describing a bird song
a<-pastew(tico,tico,f=22050)
oscillo(a,f=22050,k=1,j=1)
# a direct way to see what has been pasted
pastew(tico,tico,f=22050,plot=TRUE)
# cut a section and then paste it at the beginning
a<-cutw(tico, f=22050, from=0.5, to=0.9)
pastew(a,tico,f=22050,plot=TRUE)
# or paste it at the end
pastew(tico,a,f=22050,plot=TRUE)
# similar to
pastew(a,tico,f=22050,at=nrow(tico)/22050,plot=TRUE)
# or paste it at a specific location
pastew(a,tico,f=22050,at=1.4,plot=TRUE)
} 
\keyword{dplot}
\keyword{ts}
