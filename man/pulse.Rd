\name{pulse}

\alias{pulse}

\title{Generate rectangle pulse}

\description{
This function generates a rectangle pulse.
}

\usage{
pulse(dbefore, dpulse, dafter, f, plot = FALSE, ...)
}

\arguments{
  \item{dbefore}{duration of the silent period before the pulse}
  \item{dpulse}{duration of the pulse to generate}
  \item{dafter}{duration of silent period after the pulse}
  \item{f}{sampling frequency of the signal to be generated (in Hz)}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of the 
  pulse generated (by default \code{FALSE}).}
  \item{\dots}{other \code{\link{plot}} parameters.}
}

\value{
Returns a one column matrix corresponding to the signal generated. 
}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\seealso{\code{\link{synth}}, \code{\link{noise}}}

\examples{
pulse(dbefore=0.5,dpulse=0.1,dafter=0.3,f=8000,plot=TRUE)
}

\keyword{datagen}
\keyword{ts}
