\name{rmoffset}

\alias{rmoffset}

\title{Remove the offset of a time wave}

\description{This function removes the offset of a time wave.}

\usage{rmoffset(wave, f, FUN = mean, plot = FALSE, output = "matrix", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{FUN}{a function used to apply the offset correction. See Details.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of the wave
    after removing the offset (by default \code{FALSE}).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{The offset is removed by substracting the wave by its mean
  (argument \code{FUN}). But other function can be used. For instance, it
  can be more approriate to use the median to remove the
  offtset and transients. See Examples. \cr
   If \code{plot} is \code{FALSE}, a new wave is returned. The class
   of the returned object is set with the argument \code{output}.
 }

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}}

\examples{
data(tico)
# artifically generates an offset
tico2<-tico+0.1
# see the wave with an offset
oscillo(tico2, f=22050)
# remove the offset with the mean (by default)
rmoffset(tico2, f=22050, plot=TRUE)
# remove the offset with the median
rmoffset(tico2, f=22050, FUN=median, plot=TRUE)
}

\keyword{dplot}
\keyword{ts}
