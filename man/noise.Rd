\name{noise}

\alias{noise}

\title{Generate noise}

\description{
This function generates noise.
}

\usage{
noise(f, d, type="unif", listen = FALSE, Sample = FALSE)
}

\arguments{
  \item{f}{sampling frequency of the signal to be generated (in Hz)}
  \item{d}{duration of the signal to be generated.}
  \item{type}{a character string to specify the type of noise, either "unif" or "gaussian".}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
}

\details{
Uniform noise is generated using \code{\link{runif}} and gaussian noise is based on \code{\link{rnorm}} 
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{synth}}, \code{\link{pulse}}}

\examples{
# add noise to a synthetic signal
a<-noise(d=1,f=8000)
b<-synth(f=8000,d=1,cf=2000,plot=FALSE)
c<-a+b
spectro(c,f=8000)
}

\keyword{datagen}
\keyword{ts}
