\name{pulsew}

\alias{pulsew}

\title{Generate rectangle pulse}

\description{
  This function generates a rectangle pulse.
}

\usage{
pulsew(dbefore, dpulse, dafter, f, plot = FALSE, output = "matrix", ...)
}

\arguments{
  \item{dbefore}{duration of the silent period before the pulse}
  \item{dpulse}{duration of the pulse to generate}
  \item{dafter}{duration of silent period after the pulse}
  \item{f}{sampling frequency of the signal to be generated (in Hz)}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of the 
    pulse generated (by default \code{FALSE}).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{plot}} parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{synth}}, \code{\link{noisew}}}

\examples{
pulsew(dbefore=0.5,dpulse=0.1,dafter=0.3,f=8000,plot=TRUE)
}

\keyword{datagen}
\keyword{ts}
