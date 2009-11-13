\name{setenv}

\alias{setenv}

\title{Set the amplitude envelope of a time wave to another one}

\description{This function sets the amplitude envelope of a time wave
  to another one}

\usage{setenv(wave1, wave2, f, envt="hil", msmooth = NULL, ksmooth = NULL,
plot = FALSE, listen = FALSE, output = "matrix", ...)
}

\arguments{
  \item{wave1}{a first R object.}     
  \item{wave2}{a second R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{envt}{the type of envelope to be used for \code{wave2}: either "abs" for absolute
    amplitude envelope or "hil" for Hilbert amplitude envelope. See \code{\link{env}}.}
  \item{msmooth}{a vector of length 2 to smooth the amplitude envelope of \code{wave2} 
    with a mean sliding window. The first component is the window length
    (in number of points). The second component is the overlap between
    successive windows (in \%). See \code{\link{env}}.}
  \item{ksmooth}{kernel smooth via \code{\link{kernel}} to apply
    to the amplitude envelope of\code{wave2}. See \code{\link{env}}.}
  \item{plot}{if \code{TRUE} returns the oscillogram
    of the new time wave (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{\code{wave1} and \code{wave2} can have different duration (length)\cr
  Smoothing the envelope with \code{smooth} or \code{ksmooth} can significantly
  change the value returned.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{drawenv}}, \code{\link{env}}, \code{\link{synth}}}

\examples{
data(tico)
a<-synth(d=1,f=22050,cf=1000)
# apply 'tico' ammplitude envelope to 'a' that has a square amplitude envelope
setenv(a,tico,f=22050,plot=TRUE)
# the same but with smoothing the envelope
setenv(a,tico,f=22050,ksmooth=kernel("daniell",50),plot=TRUE)
}

\keyword{datagen}
\keyword{ts}
