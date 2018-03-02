\name{noisew}

\alias{noisew}

\title{Generate noise}

\description{
  This function generates noise.
}

\usage{
noisew(f, d, type="unif", listen = FALSE, output = "matrix")
}

\arguments{
  \item{f}{sampling frequency of the signal to be generated (in Hz)}
  \item{d}{duration of the signal to be generated.}
  \item{type}{a character string to specify the type of noise, either "unif" or "gaussian".}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
}

\details{
  Uniform noise is generated using \code{\link{runif}} and gaussian noise is based on \code{\link{rnorm}} 
}

\value{A new wave is returned. The class of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{synth}}, \code{\link{pulse}}}

\examples{
# add noise to a synthetic signal
a<-noisew(d=1,f=8000)
b<-synth(f=8000,d=1,cf=2000,plot=FALSE)
c<-a+b
spectro(c,f=8000)
}

\keyword{datagen}
\keyword{ts}
