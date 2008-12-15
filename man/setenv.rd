\name{setenv}

\alias{setenv}

\title{Set the amplitude envelope of a time wave to another one}

\description{This function sets the amplitude envelope of a time wave
to another one}

\usage{setenv(wave1, wave2, f, envt="hil", msmooth = NULL, ksmooth = NULL,
plot = FALSE, listen = FALSE, Sample = FALSE, ...)
}

\arguments{
	\item{wave1}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
	\item{wave2}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel) describing the time wave
	which envelope will be used to set \code{wave1} envelope.}
  \item{f}{sampling frequency of \code{wave1} and \code{wave2} (in Hz).
  Does not need to be specified if \code{wave1} and/or \code{wave2} are/is
	of class \code{ts}, \code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
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
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{\code{wave1} and \code{wave2} can have different duration (length)\cr
Smoothing the envelope with \code{smooth} or \code{ksmooth} can significantly
change the value returned.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

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
