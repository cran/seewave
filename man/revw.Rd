\name{revw}

\alias{revw}

\title{Time reverse of a time wave}

\description{Reverse the wave along the time axis.}

\usage{revw(wave, f, plot = FALSE, Sample = FALSE, ...)}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of
  the reversed wave (by default \code{FALSE}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.  
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{If \code{plot} is \code{TRUE} returns an oscillogram of the reversed
wave.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}}, \code{\link{deletew}},
\code{\link{fadew}}, \code{\link{pastew}}, \code{\link{mutew}}}

\examples{
data(tico)
revw(tico,f=22050,plot=TRUE)
}

\keyword{dplot}
\keyword{ts}