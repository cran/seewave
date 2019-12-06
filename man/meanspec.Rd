\name{meanspec}

\alias{meanspec}

\title{Mean frequency spectrum of a time wave}

\description{
  This function returns the mean frequency spectrum (i.e. the mean relative amplitude of the
  frequency distribution) of a time wave. Results can be expressed either in absolute 
  or dB data.
}

\usage{
meanspec(wave, f, channel = 1, wl = 512, wn = "hanning", ovlp = 0, fftw = FALSE,
norm = TRUE, PSD = FALSE, PMF = FALSE, FUN = mean, correction = "none", dB = NULL,
dBref = NULL, from = NULL, to = NULL, identify = FALSE,
col = "black", cex = 1, plot = 1, flab = "Frequency (kHz)",
alab = "Amplitude", flim = NULL, alim = NULL, type ="l", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{channel}{channel of the R object, by default left channel (1).}
  \item{wl}{length of the window for the analysis (even number of points, by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{ovlp}{overlap between two successive analysis windows (in \%).}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes of \code{spectro}.}
  \item{norm}{if \code{TRUE} the mean spectrum is normalised (i. e. scaled) by  its maximum.}
  \item{PSD}{if \code{TRUE} return Power Spectra Density,
    \emph{i. e.} the square of the spectra.}   
  \item{PMF}{if \code{TRUE} return Probability Mass Function,
    \emph{i. e.} the probability distribution of frequencie.}
  \item{FUN}{the function to apply on the rows of the STFT matrix, by
  default \code{mean} for a mean spectrum but could be other as
  \code{median} or \code{var} for a median spectrum or variance spectrum.}
  \item{correction}{a character vector of length 1 to apply an
  amplitude ("amplitude") or an energy ("energy") correction
  to the FT window. This argument is useful only when one wish to obtain
  absolute values that is when \code{norm=FALSE} and \code{PMF=FALSE}.
  By default no correction is applied ("none").}    
  \item{dB}{a character string specifying the type dB to return: "max0" for a
    maximum dB value at 0, "A", "B", "C", "D", and "ITU" for common dB weights.}
  \item{dBref}{a dB reference value when \code{dB} is not \code{NULL}. \code{NULL} by default
    but should be set to 2*10e-5 for a 20 microPa reference (SPL).}
  \item{from}{start mark where to compute the spectrum (in s).}
  \item{to}{end mark where to compute the spectrum (in s).}
  \item{identify}{to identify frequency and amplitude values on the plot
    with the help of a cursor.}
  \item{col}{colour of the spectrum.}
  \item{cex}{pitch size.}
  \item{plot}{if \code{1} returns frequency on x-axis, if \code{2} returns 
    frequency on y-axis, (by default \code{1}).}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{range of frequency axis (in kHz).}
  \item{alim}{range of amplitude axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  See examples of \code{\link{spec}}. This function is based on \code{\link{fft}}.
}

\value{
  If \code{plot} is \code{FALSE}, \code{meanspec} returns  a two columns matrix,
  the first column corresponding to the frequency axis, the second column
  corresponding to the amplitude axis.\cr 
  If \code{identify} is \code{TRUE}, \code{spec} returns a list with
  two elements:\cr 
  \item{freq}{the frequency of the points chosen on the spectrum}\cr 
  \item{amp}{the relative amplitude of the points chosen on the spectrum}\cr 
}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{The argument \code{fftw} can be used to try to speed up process
  time. When set to \code{TRUE}, the Fourier transform is computed
  through the function \code{FFT} of the package {fftw}. This pacakge is a
  wrapper around the fastest Fourier transform of the free C subroutine
  library FFTW (\url{http://www.fftw.org/}). FFT should be then installed on your OS.} 

\section{Warning}{The argument \code{peaks} is no more available
  (version > 1.5.6). See the function \code{\link{fpeaks}}
  for peak(s) detection.} 

\seealso{\code{\link{spec}},\code{\link{fpeaks}},
  \code{\link{localpeaks}}, \code{\link{dynspec}}, 
  \code{\link{corspec}}, \code{\link{diffspec}}, \code{\link{simspec}}, \code{\link{fft}}.}

\examples{
data(orni)
# compute the mean spectrum of the whole time wave
meanspec(orni,f=22050)
# compute the mean spectrum of a time wave section (from 0.32 s to 0.39 s)
meanspec(orni,f=22050,from=0.32,to=0.39)
# different window lengths
op<-par(mfrow=c(3,1))
meanspec(orni,f=22050,wl=256)
title("wl=256")
meanspec(orni,f=22050,wl=1024)
title("wl=1024")
meanspec(orni,f=22050,wl=4096)
title("wl=4096")
par(op)
# different overlap values (almost no effects here...)
op<-par(mfrow=c(3,1))
meanspec(orni,f=22050)
title("ovlp=0")
meanspec(orni,f=22050,ovlp=50)
title("ovlp=50")
meanspec(orni,f=22050,ovlp=95)
title("ovlp=95")
par(op)
# use of flim to zoom in
op<-par(mfrow=c(2,1))
meanspec(orni,f=22050)
title("zoom in")
meanspec(orni,f=22050,wl=512,flim=c(4,6))
par(op)
# comparaison of spectrum and mean spectrum
op<-par(mfrow=c(2,1))
spec(orni,f=22050)
title("spec()")
meanspec(orni,f=22050)
title("meanspec()")
par(op)
# log scale on frequency axis
meanspec(orni, f=22050, log="x")
# median spectrum
meanspec(orni,f=22050, FUN=median)
# variance spectrum
meanspec(orni,f=22050, FUN=var)
}

\keyword{dplot}
\keyword{ts}
