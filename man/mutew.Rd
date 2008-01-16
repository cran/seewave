\name{mutew}

\alias{mutew}

\title{Replace time wave data by 0 values}

\description{
This functions replaces a time wave or a section of a time wave by 0 values.
For a time wave describing a sound, this corresponds in muting the sound or a section of it. 
}

\usage{
mutew(wave, f, from = FALSE, to = FALSE, plot = TRUE, Sample = FALSE,...)
}

\arguments{
  \item{wave}{data describing the time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{from}{start of the silence section (in s).}
  \item{to}{end of the silence section (in s).}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot
      of \code{wave} with the new silence section (by default \code{TRUE}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{
By default, \code{from} and \code{from} are \code{FALSE},
this results in completely muting \code{wave}.  
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jérôme Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}}, \code{\link{cutw}},
\code{\link{deletew}},\code{\link{fadew}}, \code{\link{pastew}},
\code{\link{revw}}, \code{\link{zapsilw}}}

\examples{
data(tico)
mutew(tico,f=22050,from=0.5,to=0.9)
}

\keyword{dplot}

\keyword{ts}
