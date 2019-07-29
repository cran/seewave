\name{dynspectro}

\alias{dynspectro}

\title{Dynamic sliding spectrogramn}

\description{This function plots dynamically a sliding spectrogram along a time wave.}

\usage{
dynspectro(wave, f, channel = 1, slidframe = 10,
wl = 512, wn = "hanning", zp = 0, ovlp = 75,
fftw = FALSE, dB = TRUE, plot = TRUE,
title = TRUE, osc = FALSE,
tlab = "Time (s)", flab = "Frequency (kHz)", alab = "Amplitude",
from = NULL, to = NULL,
collevels = NULL, palette = spectro.colors,
envt = NULL, msmooth = NULL, ksmooth = NULL,
coltitle = "black", colbg = "white", colline = "black",
colaxis = "black", collab = "black", cexlab = 1,
fontlab = 1, colwave = "black",
coly0 = "lightgrey", colcursor = "red", bty = "l")
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{channel}{channel of the R object, by default left channel (1).}  
  \item{slidframe}{size of the sliding frame (in percent of the wave duration).}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{zp}{zero-padding (even number of points), see \code{Details}.}
  \item{ovlp}{overlap between two successive windows (in \% ).}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes of the \code{spectro}.}
  \item{dB}{a logical, if \code{TRUE} then uses dB values}
  \item{plot}{logical, if \code{TRUE} plots in an ew graphics device the successive
    spectrograms sliding along the time wave (by default \code{TRUE}).}
  \item{title}{logical, if \code{TRUE} adds a title with the time position of the current
    spectrogram along the time wave.}
  \item{osc}{logical, if \code{TRUE} plots an oscillogram beneath
    the sliding spectrogram with a cursor showing the position of the 
    current spectrum (by default \code{FALSE}).}
  \item{tlab}{title of the time axis.}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{from}{start mark where  to compute the sliding spectrogram (in s).}
  \item{to}{end mark where to compute the sliding spectrogram (in s).}
  \item{collevels}{a set of levels which are used to partition
    the amplitude range of the spectrogram.}
  \item{palette}{a color palette function to be used to assign colors in
    the plot.}
  \item{envt}{the type of envelope to be plooted:
    either "abs" for absolute amplitude envelope or "hil" for Hilbert amplitude envelope.
    See \code{\link{env}}.}
  \item{msmooth}{when \code{env} is not \code{NULL},
    a vector of length 2 to smooth the amplitude envelope with a 
    mean sliding window. The first component is the window length
    (in number of points). The second component is the overlap between
    successive windows (in \%). See \code{\link{env}}.}
  \item{ksmooth}{when \code{env} is not \code{NULL},
    kernel smooth via \code{\link{kernel}}. See \code{\link{env}}.}
  \item{coltitle}{if \code{title} is \code{TRUE}, colour of the title.}
  \item{colbg}{background colour.}
  \item{colline}{colour of axes line.}
  \item{colaxis}{colour of the axes.}
  \item{collab}{colour of axes title.}  
  \item{cexlab}{character size for axes title.}
  \item{fontlab}{font for axes title.}
  \item{colwave}{colour of the oscillogram or of the envelope (only when \code{osc} is \code{TRUE}).}
  \item{coly0}{colour of the y=0 line (only when \code{osc} is \code{TRUE}).}
  \item{colcursor}{colour of oscillogram cursor (only when \code{osc} is \code{TRUE}).}
  \item{bty}{the type of box to be drawn around the oscillogram (only
    when \code{osc} is \code{TRUE}).} 
}

\details{
  Use the slider panel to move along the time wave. \cr
  The function requires the package \pkg{rpanel} that is based on the
  package \pkg{tcltk}.\cr
  The function is mainly written for inspecting long sounds.\cr
  The function is based on \code{\link{image}} for fast display when
  \code{\link{spectro}} is based on \code{\link{filled.contour}}.
  Displaying the amplitude envelope with the argument \code{envt} can
  slow down significantly the display.
}


\value{
  This function returns a list of three items:
  \item{time}{a numeric vector corresponding to the time axis.}
  \item{freq}{a numeric vector corresponding to the frequency axis.}
  \item{amp}{a numeric matrix corresponding to the amplitude values.
   Each column is a Fourier transform of length \code{wl/2}.}
}

\author{David Pinaud and Jerome Sueur}

\note{This function is very similar to a spectrogram. See the \code{Details} of
  \code{\link{spectro}} for some information regarding the short term Fourier 
  transform.}

\seealso{\code{\link{spectro}}, \code{\link{spectro3D}},
  \code{\link{wf}}, \code{\link{spec}}, \code{\link{dynspec}},
  \code{\link{fft}}, \code{\link{oscillo}}, \code{\link{env}}.}

\examples{
\dontrun{
data(sheep)
require(rpanel)
dynspectro(sheep, ovlp=95, osc=TRUE)
}}

\keyword{dplot}
\keyword{ts}
