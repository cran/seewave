\name{lfs}

\alias{lfs}

\title{Linear Frequency Shift}

\description{
This function linearly shifts all the frequency content of a time wave. 
}

\usage{lfs(wave, f, shift, wl = 128, wn = "hanning", Sample = FALSE)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{shift}{positive or negative frequency shift to apply (in Hz).)}
  \item{wl}{window length for the analysis (even number of points, by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
}

\details{
A short-term Fourier transform is first applied to the signal
(see \code{\link{spectro}}), then the frequency shift is applied and the new
signal is eventually generated using the reverse of the Fourier Transform
(\code{\link{fft}}).\cr
There is therefore neither temporal modifications nor 
amplitude modifications.
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\references{
Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg.}

\author{Jerome Sueur \email{sueur@mnhn.fr} and 
Thierry Aubin \email{thierry.aubin@u-psud.fr}}

\seealso{\code{\link{ffilter}}, \code{\link{spectro}}}

\examples{
data(orni)
a<-lfs(orni,f=22050,shift=1000)
spectro(a,f=22050)
# to be compared with the original signal
spectro(orni,f=22050)
}

\keyword{ts}
