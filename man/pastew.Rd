\name{pastew}

\alias{pastew}

\title{Paste a time wave to another one}

\description{
This function pastes a first time wave to a second one. The time wave to be 
pasted, the time wave to be completed and the resulting time wave can be 
displayed in a three-frame oscillographic plot.
}

\usage{
pastew(wave1, wave2, f, at = "end", choose = FALSE, plot = FALSE,
marks = TRUE, Sample = FALSE,...)
}

\arguments{
  \item{wave1}{data or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound})
  describing the time wave to be pasted.}
  \item{wave2}{data or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound})
  describing the time wave into which \code{wave1} will be pasted.}
  \item{f}{sampling frequency of \code{wave1} and \code{wave2} (in Hz).
          Does not need to be specified if \code{wave1} and/or \code{wave2} are/is (a)
          \code{\link[sound]{Sample}} object(s).}
  \item{at}{\code{wave2} position in seconds where \code{wave1} will be pasted into.
    Can be also specified as \code{"start"}, \code{"middle"} or \code{"end"}.}
  \item{choose}{logical, if \code{TRUE} the point where \code{wave1} will
	be pasted into \code{wave2} (=\code{at}) can be graphically chosen with a cursor.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of 
  \code{wave1}, \code{wave2} and \code{wave1} + \code{wave2} (by default \code{FALSE}).}
  \item{marks}{logical, if \code{TRUE} shows where \code{wave1} has been pasted (by default \code{TRUE}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
  }

\details{If \code{plot} is \code{TRUE} returns a two-frame plot with three waves:\cr 
(1) the wave to be pasted (\code{wave1}),\cr
(2) the wave to be completed (\code{wave2}),\cr
(3) the resulting wave.}
 
\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}}, \code{\link{cutw}},
\code{\link{deletew}}, \code{\link{fadew}}, \code{\link{mutew}}
, \code{\link{revw}}, \code{\link{repw}}, \code{\link{zapsilw}}}

\examples{
data(tico)
# double a data set describing a bird song
a<-pastew(tico,tico,f=22050)
oscillo(a,f=22050)
# a direct way to see what has been pasted
pastew(tico,tico,f=22050,plot=TRUE)
# cut a section and then paste it at the beginning
a<-cutw(tico, f=22050, from=0.5, to=0.9)
pastew(a,tico,f=22050,at="start",plot=TRUE)
# or paste it at a specific location
pastew(a,tico,f=22050,at=1.4,plot=TRUE)
} 
\keyword{dplot}
\keyword{ts}
