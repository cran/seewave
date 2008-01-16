\name{fadew}

\alias{fadew}

\title{Fade in and fade out of a time wave}

\description{
This function applies a \dQuote{fade in} and/or a \dQuote{fade out} to a time wave following
a linear, exponential or cosinus-like shape.}

\usage{
fadew(wave, f, din = 0, dout = 0, shape = "linear", plot = FALSE,
listen = FALSE, Sample = FALSE, ...)}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{din}{fade in duration}
  \item{dout}{fade out duration}
  \item{shape}{fade shape, \code{"linear"}, \code{"exp"} for exponential,
  \code{"cos"} for cosinus-like, (by default \code{"linear"})}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of the wave
  modified (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jérôme Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}}, \code{\link{cutw}},
\code{\link{deletew}},\code{\link{mutew}}, \code{\link{pastew}}, \code{\link{revw}},
\code{\link{zapsilw}}
}

\examples{
a<-noise(d=5,f=4000)
op<-par(mfrow=c(3,1))
fadew(a,f=4000,din=1,dout=2,plot=TRUE,title="Linear",cexlab=0.8)
fadew(a,f=4000,din=1,dout=2,shape="exp",plot=TRUE,title="Exponential shape",
    colwave="blue",coltitle="blue",cexlab=0.8)
fadew(a,f=4000,din=1,dout=2,shape="cos",plot=TRUE,title="Cosinus-like shape",
    colwave="red",coltitle="red",cexlab=0.8)
par(op)
}

\keyword{dplot}
\keyword{ts}
