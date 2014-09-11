\name{lfs}

\alias{lfs}

\title{Linear Frequency Shift}

\description{
  This function linearly shifts all the frequency content of a time wave. 
}

\usage{lfs(wave, f, shift, wl = 1024, ovlp = 75, wn = "hanning", fftw = FALSE, output = "matrix")}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{shift}{positive or negative frequency shift to apply (in Hz).)}
  \item{wl}{window length for the analysis (even number of points, by
    default = 1024).}
  \item{ovlp}{overlap between successive FFT windows (in \%, by default 75\%).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default
    \code{"hanning"}).}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes of the \code{spectro}.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"},
    \code{"audioSample"} or \code{"ts"}.}
}

\details{
  A short-term Fourier transform is first applied to the signal
  (see \code{\link{spectro}}), then the frequency shift is applied and the new
  signal is eventually generated using the reverse of the Fourier Transform
  (\code{\link{istft}}).\cr
  There is therefore neither temporal modifications nor 
  amplitude modifications.
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

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
