\name{rmnoise}

\alias{rmnoise}

\title{Remove noise}

\description{
  This function removes background noise by smoothing 
}

\usage{
rmnoise(wave, f, output = "matrix", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link[stats]{smooth.spline}} arguments.}
}

\details{
  This function is based on \code{\link[stats]{smooth.spline}}. You can
  use the arguments of the later to modify the smoothing.
}


\value{A new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{Low frequency noise might not be removed out properly.}

\seealso{\code{\link{afilter}},  \code{\link{noisew}}}

\examples{
# synthesis of a 440 Hz sound with background noise
n <- noisew(d=1,f=8000)
s <- synth(d=1,f=8000,cf=440)
ns <- n+s
# remove noise (but low frequency content still there)
a <- rmnoise(ns,f=8000)
}

\keyword{ts}

