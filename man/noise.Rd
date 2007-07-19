\name{noise}

\alias{noise}

\title{Generate noise}

\description{
This function generates noise.
}

\usage{
noise(f, d)
}

\arguments{
  \item{f}{sampling frequency of the signal to be generated (in Hz)}
  \item{d}{duration of the signal to be generated}
}

\details{
\code{noise} is based on \code{\link{runif}}
}

\value{
Returns a one column matrix with random data (by default: mean = 0, sd = 1). 
}


\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

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
