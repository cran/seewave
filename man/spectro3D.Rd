\name{spectro3D}

\alias{spectro3D}

\title{3D-spectrogram of a time wave}

\description{
  This function returns a three-dimension spectrographic representation of a
  time wave. The function corresponds to short-term Fourier transform.
}

\usage{
spectro3D(wave, f, wl = 512, wn = "hanning", zp = 0,
ovlp = 0, norm = TRUE, fftw = FALSE, dB = "max0", dBref = NULL, plot = TRUE,
magt = 10, magf = 10, maga = 2,
palette = reverse.terrain.colors)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{length of the window for the analysis
    (even number of points).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{zp}{zero-padding (even number of points), see \code{Details}.}
  \item{ovlp}{overlap between two successive windows (in \% ).}
  \item{norm}{if \code{TRUE} the STFT is normalised (i. e. scaled) by its maximum.}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes of the \code{spectro}.}
  \item{dB}{a character string specifying the type dB to return:
    "max0" for a maximum dB value at 0, "A", "B", "C" and "D" for common dB weights.}
  \item{dBref}{a dB reference value when \code{dB} is \code{TRUE}. \code{NULL} by default
    but should be set to 2*10e-5 for a 20 microPa reference.}
  \item{plot}{logical, if \code{TRUE} plots the spectrogram
    (by default \code{TRUE}).}
  \item{magt}{magnification of the time axis.}
  \item{magf}{magnification of the frequency axis.}
  \item{maga}{magnification of the amplitude axis.}
  \item{palette}{a color palette function to be used to assign colors in
    the plot, see \code{Details}.}
}

\details{
  Following Heisenberg uncertainty principle, the short-term Fourier transform
  cannot be precised in both time and frequency. The temporal and frequency
  precisions of the function are actually dependent of the \code{wl} value.
  Choosing a high \code{wl} value will increase the frequency resolution but
  reduce the temporal one, and \emph{vice versa}. The frequency precision is
  obtained by calculating the ratio \code{f}/\code{wl},
  and the temporal precision is obtained by calculating the reverse ratio
  \code{wl}/\code{f}. This problem can be reduced in some way with \code{zp} that
  adds 0 values on both sides of the analysis window. This increases frequency
  resolution without altering time resolution.\cr
  Any colour palette can be used. In particular, it is possible to use other palettes coming with
  \pkg{seewave}: \code{reverse.gray.colors.1}, \code{reverse.gray.colors.2},
  \code{spectro.colors}, \code{temp.colors}, \code{reverse.heat.colors}, 
  \code{reverse.cm.colors}, \code{reverse.topo.colors}, corresponding to the reverse
  of \code{heat.colors},\code{topo.colors},
  \code{cm.colors}.\cr
  Use \code{magt}, \code{magf} and \code{maga} to resize the plot.   
}

\value{
  This function returns a list of three items:
  \item{time}{a numeric vector corresponding to the time axis.}
  \item{freq}{a numeric vector corresponding to the frequency axis.}
  \item{amp}{a numeric matrix corresponding to the amplitude
  values. Each column is a Fourier transform of length \code{wl}/2.}
}

\author{Jerome Sueur \email{sueur@mnhn.fr} and 
  Caroline Simonis \email{csimonis@mnhn.fr}.}

\note{This function requires \pkg{rgl} and is based on \code{\link{fft}}.
  See examples of \code{\link{spectro}} for analysis arguments
  (\code{wl},\code{zp}, \code{ovlp}).}

\seealso{\code{\link{spectro}}, \code{\link{dynspec}},
  \code{\link{wf}}, \code{\link{fft}}.}

\examples{
\dontrun{
require(rgl)
data(tico)
spectro3D(tico,f=22050,wl=512,ovlp=75,zp=16,maga=4,palette=reverse.terrain.colors)
# linear amplitude scale without a normisation of the STFT matrix
# time and frequency scales need to be dramatically amplified
spectro3D(tico, norm=FALSE, dB=NULL, magt=100000, magf=100000)
}}

\keyword{dplot}
\keyword{ts}
