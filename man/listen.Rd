\name{listen}

\alias{listen}

\title{Play a sound wave}

\description{Play a sound wave}

\usage{listen(wave, f = NULL, from = NULL, to = NULL, choose = FALSE)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{from}{start of play (in s).}
  \item{to}{end of play (in s).}
  \item{choose}{logical, if \code{TRUE} start (=\code{from}) and end (=\code{to})
    points can be graphically chosen with a cursor on the oscillogram.}
}

\author{Jerome Sueur \email{sueur@mnhn.fr} but the original
  \code{\link[tuneR]{play}} function is by Uwe Ligges (package \pkg{tuneR}).}

\note{This function is based on \code{\link[tuneR]{play}} but allows to read
  one-colum matrix, data.frame, time-series and Sample objects.}

\seealso{\code{\link[tuneR]{play}}}

\examples{
## NOT RUN
# data(tico)
# listen(tico,f=22050)
# listen(tico,f=22050,from=0.5,to=1.5)
# listen(noise(d=1,f=8000,Wave=TRUE))
## change f to play the sound a different speed
# data(sheep)
## normal
# listen(sheep,f=8000)
## two times faster
# listen(sheep,f=8000*2)
## two times slower
# listen(sheep,f=8000/2)
}

\keyword{ts}
