\name{listen}

\alias{listen}

\title{Play a sound wave}

\description{Play a sound wave}

\usage{listen(wave, f, from, to)}

\arguments{
  \item{wave}{data describing a time wave or a \code{\link[sound]{Sample}}
  object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need
  to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{from}{start of play (in s).}
  \item{to}{end of play (in s).}
}

\author{Jérôme Sueur \email{sueur@mnhn.fr} but the original
\code{\link[sound]{play}} function is by Matthias Heymann (package \pkg{sound}).}

\note{This function is based on \code{\link[sound]{play}} but allows to read
one-colum matrix, data.frame and Sample objects.}

\seealso{\code{\link[sound]{play}}}

\examples{
data(tico)
listen(tico,f=22050)
listen(tico,f=22050,from=0.5,to=1.5)
listen(noise(d=1,f=8000,Sample=TRUE))
# change f to play the sound a different speed
data(sheep)
# normal
listen(sheep,f=8000)
# two times faster
listen(sheep,f=8000*2)
# two times slower
listen(sheep,f=8000/2)
}

\keyword{ts}
