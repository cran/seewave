\name{revw}

\alias{revw}

\title{Time reverse of a time wave}

\description{Reverse the wave along the time axis.}

\usage{revw(wave, f, env = TRUE, ifreq = TRUE, plot = FALSE, Sample = FALSE, ...)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{env}{logical, if \code{TRUE} the amplitude envelope is reversed.}
	\item{ifreq}{logical, if \code{TRUE} the instantaneous frequency is reversed.}
	\item{plot}{logical, if \code{TRUE} returns an oscillographic plot of
  the reversed wave (by default \code{FALSE}).}
  \item{Sample}{logical, if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.  
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{If \code{plot} is \code{TRUE} returns an oscillogram of the reversed
wave. The amplitude and the instantaneous frequency can be independently reversed
thanks to the arguments \code{env} and \code{ifreq}. See the examples.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

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