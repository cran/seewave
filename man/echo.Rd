\name{echo}

\alias{echo}

\title{Echo generator}
\description{This function generate echoes of a time wave.}

\usage{echo(wave, f, amp, delay, plot = FALSE,
listen = FALSE, output = "matrix", ...)}

\arguments{
\item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{amp}{a vector describing the relative amplitude
    of the successive echoes. Each value of the vector should be in [0,1]}
  \item{delay}{a vector describing the time delays of the successive echoes
    from the beginning of \code{wave} (in s.)}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of the wave
    modified (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
  \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{
\code{amp} and \code{delay} should strictly have the same length corresponding
to the number of desired echoes.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
of the returned object is set with the argument \code{output}.}

\references{Stoddard, P. K. (1998). Application of filters in bioacoustics.
\emph{In}: Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds), \emph{Animal acoustic
communication}. Springer, Berlin, Heidelberg,pp. 105-127.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{This function is based on a convolution (\code{\link{convolve}}) between the
input wave and a pulse echo filter.}

\seealso{\code{\link{synth}}}

\examples{
# generation of the input wave
a <- synth(f=11025,d=1,cf=2000,shape="tria",am=c(50,10),fm=c(1000,10,1000))
# generation of three echoes
# with respectively a relative amplitude of 0.8, 0.4, and 0.2
# and with a delay of 1s, 2s, and 3s  from the beginning of the input wave
aecho <- echo(a,f=11025,amp=c(0.8,0.4,0.2),delay=c(1,2,3))
# another echo with time delays overlapping with the input wave
aecho <- echo(a,f=11025,amp=c(0.4,0.2,0.4),delay=c(0.6,0.8,1.5))
}

\keyword{datagen}
\keyword{ts}
