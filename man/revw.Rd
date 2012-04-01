\name{revw}

\alias{revw}

\title{Time reverse of a time wave}

\description{Reverse the wave along the time axis.}

\usage{revw(wave, f, env = TRUE, ifreq = TRUE,
plot = FALSE, output = "matrix", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{env}{logical, if \code{TRUE} the amplitude envelope is reversed.}
  \item{ifreq}{logical, if \code{TRUE} the instantaneous frequency is reversed.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of
    the reversed wave (by default \code{FALSE}).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{If \code{plot} is \code{TRUE} returns an oscillogram of the reversed
  wave. The amplitude and the instantaneous frequency can be independently reversed
  thanks to the arguments \code{env} and \code{ifreq}. See the examples.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\references{
  Beeman, K. 1998. Digital signal analysis, editing and synthesis \emph{in} 
  Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
    communication}, pp. 59-103. Springer, Berlin, Heidelberg.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}}, \code{\link{deletew}},
  \code{\link{fadew}}, \code{\link{pastew}}, \code{\link{mutew}}}

\examples{
data(tico)
# simple reverse
revw(tico,f=22050,plot=TRUE)
# envelope reverse only
revw(tico,f=22050,ifreq=FALSE, plot=TRUE)
# instantaneous frequency reverse only
revw(tico,f=22050,env=FALSE, plot=TRUE)
}

\keyword{dplot}
\keyword{ts}