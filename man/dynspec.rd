\name{dynspec}

\alias{dynspec}

\title{Dynamic sliding spectrum}

\description{
This function plots dynamically a sliding spectrum along a time wave. 
This basically corresponds to a short-term Fourier transform.
}

\usage{
dynspec(wave, f, wl = 512, wn = "hanning", zp = 0, ovlp = 0, norm = FALSE,
dB = FALSE, plot = TRUE, title = TRUE, osc = FALSE, flab = "Frequency (kHz)",
alab = "Amplitude", alim = NULL, flim = c(0, f/2000), type = "l", from = FALSE,
to = FALSE, env = FALSE, smooth = NULL, ksmooth = NULL, colspec = "black",
coltitle = "black", colbg = "white", colline = "black", colaxis = "black",
collab = "black", cexlab = 1, fontlab = 1, colwave = "black",
coly0 = "lightgrey", colcursor = "red", bty = "l")
}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave}
  is a \code{\link[sound]{Sample}} object.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{zp}{zero-padding (even number of points), see \code{Details}.}
  \item{ovlp}{overlap between two successive windows (in \% ).}
  \item{norm}{logical, if \code{TRUE} compute a normalised sliding spectrum.}
  \item{dB}{logical, if \code{TRUE} returns the sliding spectrum in dB 
      (by default \code{FALSE}).}
  \item{plot}{logical, if \code{TRUE} plots in an ew graphics device the successive
      spectra sliding along the time wave (by default \code{TRUE}).}
  \item{title}{logical, if \code{TRUE} adds a title with the time position of the current
      spectrum along the time wave.}
  \item{osc}{logical, if \code{TRUE} plots an oscillogram beneath
      the sliding spectrum with a cursor showing the position of the 
      current spectrum (by default \code{FALSE}).}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{range of frequency axis.}
  \item{alim}{range of amplitude axis.}
  \item{type}{type of plot that should be drawn for the sliding spectrum.
      See \code{\link{plot}} for details (by default "l" for lines).}
  \item{from}{start mark where  to compute the sliding spectrum (in s).}
  \item{to}{end mark where to compute the sliding spectrum (in s).}
  \item{env}{logical, if \code{TRUE} plots the envelope of the time wave instead 
      (by default \code{FALSE}) (only when \code{osc} is \code{TRUE}).}
  \item{smooth}{if \code{env} is \code{TRUE} then smoothes the envelope 
      by floating average (only when \code{osc} is \code{TRUE}).}
  \item{ksmooth}{if \code{env} is \code{TRUE} then smoothes the envelope 
      by applying a \code{\link{kernel}} smoothing (only when \code{osc} is \code{TRUE}).}
  \item{colspec}{colour of the sliding spectrum.}
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
  \item{bty}{the type of box to be drawn around the oscillogram (only when \code{osc} is \code{TRUE}).} 
}

\details{
Use the slider panel to move along the time wave.\cr
Use the argument \code{norm} if you wish to have each spectrum normalised, \emph{i.e.}
with values between 0 and 1 or maximised to 0 dB when \code{dB} is \code{TRUE}.\cr
The function requires the package \pkg{rpanel} that is based on the package \pkg{tcltk}.
}

\value{
If \code{plot} is \code{FALSE}, this function returns a matrix which columns
correspond to the spectra computed along the time wave.}

\author{Jérôme Sueur \email{sueur@mnhn.fr} and 
Caroline Simonis \email{csimonis@mnhn.fr}.}

\note{This function is very similar to a spectrogram. See the \code{Details} of
\code{\link{spectro}} for some information regarding the short term Fourier 
transform.}

\seealso{\code{\link{spectro}}, \code{\link{spectro3D}}, \code{\link{spec}},
\code{\link{fft}}, \code{\link{oscillo}}.}

\examples{
data(sheep)
dynspec(sheep,f=8000,wl=1024,ovlp=50,osc=TRUE)
dev.off()
}

\keyword{dplot}
\keyword{ts}