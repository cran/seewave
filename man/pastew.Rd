\name{pastew}

\alias{pastew}

\title{Paste a time wave to another one}

\description{
  This function pastes a first time wave to a second one. The time wave to be 
  pasted, the time wave to be completed and the resulting time wave can be 
  displayed in a three-frame oscillographic plot.
}

\usage{
pastew(wave1, wave2, f, at = "end",
join = FALSE, tjunction=0,
choose = FALSE, plot = FALSE,
marks = TRUE, output = "matrix", ...)
}

\arguments{
  \item{wave1}{a first R object.}
  \item{wave2}{a second R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{at}{\code{wave2} position in seconds where \code{wave1} will be pasted into.
    Can be also specified as \code{"start"}, \code{"middle"} or
    \code{"end"}.}
  \item{join}{if \code{TRUE} the two waves will be pasted and jointed by
    removing the last point of wave2. See examples.}
  \item{tjunction}{a numeric vector to remove clicks at the junction of ‘wave1’ and ‘wave2’. The value specifies the duration in seconds where the real vales will be replaced by a linear interpolation. This duration should be a few milliseconds.}
  \item{choose}{logical, if \code{TRUE} the point where \code{wave1} will
    be pasted into \code{wave2} (=\code{at}) can be graphically chosen with a cursor.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of 
    \code{wave1}, \code{wave2} and \code{wave1} + \code{wave2} (by default \code{FALSE}).}
  \item{marks}{logical, if \code{TRUE} shows where \code{wave1} has been
    pasted (by default \code{TRUE}).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{If \code{plot} is \code{TRUE} returns a two-frame plot with three waves:\cr 
  (1) the wave to be pasted (\code{wave1}),\cr
  (2) the wave to be completed (\code{wave2}),\cr
  (3) the resulting wave.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur, improved by Laurent Lellouch}

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
# setting the argument 'join' to TRUE might be useful
# to smooth pasting when some phase problem occur
# generate two sine waves
a <- synth(cf=50, f=400, d=0.1)
b <- synth(cf=100, f=400, d=0.1)
# paste it with 'join' turned to FALSE
# there is a click at the junction between the two waves
pastew(a, b, f=400, plot=TRUE)
# that can be removed by setting 'join' to TRUE
pastew(a, b, f=400, join=TRUE, plot=TRUE)
# or by using the argument 'tjunction'
pastew(a, b, f=400, tjunction=0.01, plot=TRUE)
} 
\keyword{dplot}
\keyword{ts}
