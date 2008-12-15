\name{drawenv}

\alias{drawenv}

\title{Draw the amplitude envelope of a time wave}

\description{This function lets the user modifying the amplitude envelope
of a time wave by drawing it with the graphics device }

\usage{
drawenv(wave, f, n = 20, plot = FALSE, listen = FALSE, Sample = FALSE)
}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{n}{the maximum number of points to draw the new envelope. Valid values start at 1.}
  \item{plot}{if \code{TRUE} returns the oscillogram
      of the new time wave (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}.}
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

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{setenv}}, \code{\link{env}}, \code{\link{synth}}}

\examples{
a<-synth(d=1,f=22050,cf=1000)
# drawenv(a,f=22050,plot=TRUE)
# choose points on the oscillogram view to draw a new enveloppe
# stop (ESC on Windows; right mouse button on Linux)
# check the result on the second graphics device opened thanks to plot=TRUE
}

\keyword{datagen}
\keyword{ts}
