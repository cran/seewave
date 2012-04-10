\name{mutew}

\alias{mutew}

\title{Replace time wave data by 0 values}

\description{
  This functions replaces a time wave or a section of a time wave by 0 values.
  For a time wave describing a sound, this corresponds in muting the sound or a section of it. 
}

\usage{
mutew(wave, f, from = NULL, to = NULL, choose = FALSE, plot = TRUE,
output = "matrix", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{from}{start of the silence section (in s).}
  \item{to}{end of the silence section (in s).}
  \item{choose}{logical, if \code{TRUE} start (=\code{from}) and end (=\code{to})
    points can be graphically chosen with a cursor on the oscillogram.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot
    of \code{wave} with the new silence section (by default \code{TRUE}).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{
  By default, \code{from} and \code{from} are \code{NULL},
  this results in completely muting \code{wave}.  
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}}, \code{\link{cutw}},
  \code{\link{deletew}},\code{\link{fadew}}, \code{\link{pastew}},
  \code{\link{revw}}, \code{\link{zapsilw}}}

\examples{
data(tico)
mutew(tico,f=22050,from=0.5,to=0.9)
}

\keyword{dplot}

\keyword{ts}
