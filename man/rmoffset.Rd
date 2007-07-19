\name{rmoffset}

\alias{rmoffset}

\title{Remove the offset of a time wave}

\description{This function removes the offset of a time wave.}

\usage{rmoffset(wave, f, plot = FALSE, ...)}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of the wave
  after removing the offset (by default \code{FALSE}).}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{
If \code{plot} is \code{FALSE}, \code{rmoffset} returns a one-column matrix
describing the new wave with the same sampling frequency as the original wave.
}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

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
