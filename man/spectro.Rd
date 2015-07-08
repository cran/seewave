\name{spectro}

\alias{spectro}

\title{2D-spectrogram of a time wave}

\description{
  This function returns a two-dimension spectrographic representation of a time
  wave. The function corresponds to short-term Fourier transform.
  An amplitude contour plot can be overlaid.}

\usage{
spectro(wave, f, wl = 512, wn = "hanning", zp = 0,
ovlp = 0, complex = FALSE, norm = TRUE, fftw = FALSE,
dB = "max0", dBref = NULL, plot = TRUE,
grid = TRUE, osc = FALSE, scale = TRUE, cont = FALSE,
collevels = NULL, palette = spectro.colors,
contlevels = NULL, colcont = "black",
colbg = "white", colgrid = "black",
colaxis = "black", collab="black",
cexlab = 1, cexaxis = 1, 
tlab = "Time (s)",
flab = "Frequency (kHz)",
alab = "Amplitude",
scalelab = "Amplitude\n(dB)",
main=NULL,
scalefontlab = 1, scalecexlab =0.75,
axisX = TRUE, axisY = TRUE, tlim = NULL, trel = TRUE,
flim = NULL, flimd = NULL,
widths = c(6,1), heights = c(3,1),
oma = rep(0,4),
listen=FALSE,
...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}  
  \item{wl}{window length for the analysis
    (even number of points) (by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default
    \code{"hanning"}).}
  \item{zp}{zero-padding (even number of points), see \code{Details}.}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{complex}{if \code{TRUE} the STFT will be returned as complex
  numbers.}  
  \item{norm}{if \code{TRUE} the STFT is normalised (i. e. scaled) by its maximum.}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes.}
  \item{dB}{a character string specifying the type dB to return: "max0"
    (default) for a maximum dB value at 0, "A", "B", "C" and "D" for
  common dB weights. If set to \code{NULL}, then a linear scale is used.}
  \item{dBref}{a dB reference value. \code{NULL} by default
    but should be set to 2*10e-5 for a 20 microPa reference.}
  \item{plot}{logical, if \code{TRUE} plots the spectrogram
    (by default \code{TRUE}).}
  \item{grid}{logical, if \code{TRUE} plots a y-axis grid
    (by default \code{TRUE}).}
  \item{osc}{logical, if \code{TRUE} plots an oscillogram beneath
    the spectrogram (by default \code{FALSE}).}
  \item{scale}{logical, if \code{TRUE} plots a dB colour scale on the right
    side of the spectrogram (by default \code{TRUE}).}
  \item{cont}{logical, if \code{TRUE} overplots contour lines on the spectrogram
    (by default \code{FALSE}).}
  \item{collevels}{a set of levels which are used to partition
    the amplitude range of the spectrogram (in dB).}
  \item{palette}{a color palette function to be used to assign colors in
    the plot, see \code{Details}.}
  \item{contlevels}{a set of levels which are used to partition
    the amplitude range for contour overplot (in dB).}
  \item{colcont}{colour for \code{cont} plotting.}
  \item{colbg}{background colour.}
  \item{colgrid}{colour for \code{grid} plotting.}
  \item{colaxis}{color of the axes.}
  \item{collab}{color of the labels.}
  \item{cexlab}{size of the labels.}
  \item{cexaxis}{size of the axes.}
  \item{tlab}{label of the time axis.}
  \item{flab}{label of the frequency axis.}
  \item{alab}{label of the amplitude axis.}
  \item{scalelab}{amplitude scale label.}
  \item{main}{label of the main title.}
  \item{scalefontlab}{font of the amplitude scale label.}
  \item{scalecexlab}{cex of the amplitude scale label.}
  \item{axisX}{logical, if \code{TRUE} plots time X-axis (by default \code{TRUE}).}
  \item{axisY}{logical, if \code{TRUE} plots frequency Y-axis (by default \code{TRUE}).}
  \item{tlim}{modifications of the time X-axis limits.}
  \item{trel}{time X-axis with a relative scale when \code{tlim} is not null,
    \emph{i.e.} relative to \code{wave}.} 
  \item{flim}{modifications of the frequency Y-axis limits (in kHz).} 
  \item{flimd}{dynamic modifications of the frequency Y-axis limits. New \code{wl}
    and \code{ovlp} arguments are applied to increase time/frequency resolution.}
  \item{widths}{a vector of length 2 to control the relative widths of columns on
  the device when \code{scale} is \code{TRUE}.}
  \item{heights}{a vector of length 2 to control the relative heights of rows on
    the device when \code{osc} is \code{TRUE}.}
  \item{oma}{a vector of length 4 to control the size of outer margins
  when either \code{scale} or \code{osc} is \code{TRUE}.}
  \item{listen}{if \code{TRUE} the sound is played back (by default
  \code{FALSE}).}
  \item{\dots}{other \code{\link{contour}} and \code{\link{oscillo}}
    graphical parameters.}
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
  Any colour palette can be used. In particular, it is possible to use other
  palettes coming with \pkg{seewave}: \code{temp.colors}, 
  \code{reverse.gray.colors.1},
  \code{reverse.gray.colors.2}, \code{reverse.heat.colors},
  \code{reverse.terrain.colors},
  \code{reverse.topo.colors},
  \code{reverse.cm.colors} corresponding to the reverse of \code{heat.colors}, 
  \code{terrain.colors}, \code{topo.colors}, \code{cm.colors}.\cr
  Use \code{\link{locator}} to identify  points.
}

\note{The argument \code{fftw} can be used to try to speed up process
  time. When set to \code{TRUE}, the Fourier transform is computed
  through the function \code{FFT} of the package {fftw}. This pacakge is a
  wrapper around the fastest Fourier transform of the free C subroutine
  library FFTW (\url{http://www.fftw.org/}). FFT should be then installed on your OS.} 

\value{
  This function returns a list of three items:
  \item{time}{a numeric vector corresponding to the time axis.}
  \item{freq}{a numeric vector corresponding to the frequency axis.}
  \item{amp}{a numeric or a complex matrix corresponding to the amplitude values.
   Each column is a Fourier transform of length \code{wl/2}.}
}

\references{
  Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
    communication}. Springer, Berlin, Heidelberg.}

\author{Jerome Sueur and Caroline Simonis.}

\note{This function is based on \code{\link{fft}}, \code{\link{contour}} and
  \code{\link{filled.contour}}}

\seealso{\code{\link{ggspectro}}, \code{\link{spectro3D}}, \code{\link{dynspec}}, \code{\link{wf}}, 
  \code{\link{oscillo}}, \code{\link{dBscale}}, \code{\link{fft}}.}

\examples{
\dontrun{
data(tico)
data(pellucens)
# simple plots
spectro(tico,f=22050)
spectro(tico,f=22050,osc=TRUE)
spectro(tico,f=22050,scale=FALSE)
spectro(tico,f=22050,osc=TRUE,scale=FALSE)
# change the dB scale by setting a different dB reference value (20microPa)
spectro(tico,f=22050, dBref=2*10e-5)
# unnormalised spectrogram with a linear amplitude scale
spectro(tico, dB=NULL, norm=FALSE, scale=FALSE)
# manipulating wl
op<-par(mfrow=c(2,2))
spectro(tico,f=22050,wl=256,scale=FALSE)
title("wl = 256")
spectro(tico,f=22050,wl=512,scale=FALSE)
title("wl = 512")
spectro(tico,f=22050,wl=1024,scale=FALSE)
title("wl = 1024")
spectro(tico,f=22050,wl=4096,scale=FALSE)
title("wl = 4096")
par(op)
# vertical zoom using flim
spectro(tico,f=22050, flim=c(2,6))
spectro(tico,f=22050, flimd=c(2,6))
# a full plot
pellu2<-cutw(pellucens,f=22050,from=1,plot=FALSE)
spectro(pellu2,f=22050,ovlp=85,zp=16,osc=TRUE,
    cont=TRUE,contlevels=seq(-30,0,20),colcont="red",
    lwd=1.5,lty=2,palette=reverse.terrain.colors)
# black and white spectrogram 
spectro(pellu2,f=22050,ovlp=85,zp=16,
    palette=reverse.gray.colors.1)
# colour modifications
data(sheep)
spectro(sheep,f=8000,palette=temp.colors,collevels=seq(-115,0,1))
spectro(pellu2,f=22050,ovlp=85,zp=16,
palette=reverse.cm.colors,osc=TRUE,colwave="orchid1") 
spectro(pellu2,f=22050,ovlp=85,zp=16,osc=TRUE,palette=reverse.heat.colors,
colbg="black",colgrid="white", colwave="white",colaxis="white",collab="white")
}
}

\keyword{dplot}
\keyword{ts}
