\name{noise}

\alias{noise}

\title{Generate noise}

\description{
This function generates noise.
}
\usage{
noise(f, d, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{sampling frequency of the signal to be generated}
  \item{d}{duration of the signal to be generated}
  \item{\dots}{other \code{\link{rnorm}} parameters.}
}

\details{
\code{am} is based on \code{\link{rnorm}}
}

\value{
Returns a one column matrix with random data (by defalut: mean = 0, sd = 1). 
}


\author{Jérôme Sueur \email{jerome.sueur@ibaic.u-psud.fr}}

\seealso{\code{\link{synth}}}

\examples{
# add noise to a synthetic signal
a<-noise(d=1,f=8000)
b<-synth(f=8000,d=1,cf=2000,plot=FALSE)
c<-a+b
spectro(c,f=8000,wl=512)
}
\keyword{datagen}
\keyword{ts}
