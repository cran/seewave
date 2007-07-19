\name{mute}

\alias{mute}

\title{Replace time wave data by 0 values}

\description{
This functions replaces a time wave or a section of a time wave by 0 values.
For a time wave describing a sound, this corresponds in muting the sound or a section of it. 
}

\usage{
mute(wave, f, from = FALSE, to = FALSE, plot = TRUE, ...)
}

\arguments{
  \item{wave}{data describing the time wave
  or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{from}{start of the silence section (in s).}
  \item{to}{end of the silence section (in s).}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot
      of \code{wave} with the new silence section (by default \code{TRUE}).}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{
By default, \code{from} and \code{from} are \code{FALSE},
this results in completely muting \code{wave}.  
}

\value{New wave is returned as one-column matrix if \code{plot} is \code{FALSE}.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{cutw}},
\code{\link{pastew}}, \code{\link{deletew}}, \code{\link{addsilw}}}

\examples{
data(tico)
mute(tico,f=22050,from=0.5,to=0.9)
}

\keyword{dplot}

\keyword{ts}
