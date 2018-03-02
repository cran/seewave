\name{drawenv}

\alias{drawenv}

\title{Draw the amplitude envelope of a time wave}

\description{This function lets the user modifying the amplitude envelope
  of a time wave by drawing it with the graphics device }

\usage{
drawenv(wave, f, n = 20, plot = FALSE, listen = FALSE, output = "matrix")
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{n}{the maximum number of points to draw the new envelope. Valid values start at 1.}
  \item{plot}{if \code{TRUE} returns the oscillogram
    of the new time wave (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
}

\details{
  The function first plots an oscillogram view of \code{wave}.\cr
  The user has then to choose points on the positive side of the y-axis (amplitude). 
  The junction of these points will draw a new amplitude envelope.\cr
  The order of points along the x-axis (time) is not important
  but points cannot be cancelled. When this process is finished
  the new time wave is returned in the console or as an oscillogram
  in a second graphics device if \code{plot} is \code{TRUE}.
  \cr
  The function uses \code{\link{locator}}.
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{setenv}}, \code{\link{env}}, \code{\link{synth}}}

\examples{
\dontrun{
a<-synth(d=1,f=22050,cf=1000)
# drawenv(a,f=22050,plot=TRUE)
# choose points on the oscillogram view to draw a new enveloppe
# stop (ESC on Windows; right mouse button on Linux)
# check the result on the second graphics device opened thanks to plot=TRUE
}
}

\keyword{datagen}
\keyword{ts}
