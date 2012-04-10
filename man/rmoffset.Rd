\name{rmoffset}

\alias{rmoffset}

\title{Remove the offset of a time wave}

\description{This function removes the offset of a time wave.}

\usage{rmoffset(wave, f, plot = FALSE, output = "matrix", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of the wave
    after removing the offset (by default \code{FALSE}).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}}

\examples{
data(tico)
# artifically generates an offset
tico2<-tico+0.1
# see the wave with an offset
oscillo(tico2,f=22050)
# remove the offset
rmoffset(tico2,f=22050,plot=TRUE)
}

\keyword{dplot}
\keyword{ts}
