\name{spectro}

\alias{spectro}

\title{2D-spectrogram of a time wave}

\description{
This function returns a two-dimension spectrographic representation of a time
wave. The function corresponds to short-term Fourier transform.
An amplitude contour plot can be overlaid.}

\usage{
spectro(wave, f, wl = 512, wn = "hanning", zp = 0, ovlp = 0, plot = TRUE,
grid = TRUE, osc = FALSE, scale = TRUE, cont = FALSE,
collevels = seq(-30, 0, 1), palette = spectro.colors,
contlevels = seq(-30, 0, 10), colcont = "black", colbg = "white", colgrid = "black",
colaxis = "black", collab="black",
plot.title = title(main = "", xlab = "Time (s)", 
ylab = "Frequency (kHz)"), scalelab = "Amplitude\n(dB)",
scalefontlab = 1, scalecexlab =0.75, axisX = TRUE, axisY = TRUE, tlim = NULL,
flim = NULL, flimd = NULL, listen=FALSE,
...)
}

\arguments{
  \item{wave}{data describing a time wave
     or a \code{\link[sound]{Sample}} object generated loading a wav file
     with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
     Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{wl}{window length for the analysis
    (even number of points). (by default = 512)}
  \item{zp}{zero-padding (even number of points), see \code{Details}.}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
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
  \item{plot.title}{statements which add titles to the plot.}
  \item{scalelab}{amplitude scale label.}
  \item{scalefontlab}{font of the amplitude scale label.}
  \item{scalecexlab}{cex of the amplitude scale label.}
  \item{axisX}{logical, if \code{TRUE} plots time X-axis (by default \code{TRUE}).}  
  \item{axisY}{logical, if \code{TRUE} plots frequency Y-axis (by default \code{TRUE}).}
  \item{tlim}{modifications of the time X-axis limits.} 
  \item{flim}{modifications of the frequency Y-axis limits.} 
  \item{flimd}{dynamic modifications of the frequency Y-axis limits. New \code{wl}
      and \code{ovlp} arguments are applied to increase time/frequency resolution.}    
  \item{listen}{if \code{TRUE} the sound is played back (by default \code{FALSE}).}
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
palettes coming with \pkg{seewave}: \code{rev.gray.colors.1}, \code{rev.gray.colors.2}, 
\code{rev.heat.colors}, \code{rev.terrain.colors}, \code{rev.topo.colors},
\code{rev.cm.colors} corresponding to the reverse of \code{heat.colors}, 
\code{terrain.colors}, \code{topo.colors}, \code{cm.colors}.\cr
Use \code{\link{locator}} to identify  points.
}

\value{
If \code{plot} is \code{FALSE}, this function returns a matrix. Each column
corresponds to a Fourier transform of length \code{wl}/2.}

\references{
Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg.}

\author{Jérôme Sueur \email{sueur@mnhn.fr} and 
Caroline Simonis \email{csimonis@mnhn.fr}.}

\note{This function is based on \code{\link{fft}}, \code{\link{contour}} and
\code{\link{filled.contour}}}

\seealso{\code{\link{spectro3D}}, \code{\link{dynspec}}, 
\code{\link{oscillo}}, \code{\link{dBscale}}, \code{\link{fft}}.}

\examples{
data(tico)
data(pellucens)
# simple plots
spectro(tico,f=22050)
spectro(tico,f=22050,osc=TRUE)
spectro(tico,f=22050,scale=FALSE)
spectro(tico,f=22050,osc=TRUE, scale=FALSE)
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
spectro(tico,f=22050, ylim=c(2,6))
spectro(tico,f=22050, ylimd=c(2,6))
# a full plot
pellu2<-cutw(pellucens,f=22050,from=1,plot=FALSE)
spectro(pellu2,f=22050,ovlp=85,zp=16,osc=TRUE,
    cont=TRUE,contlevels=seq(-30,0,20),colcont="red",
    lwd=1.5,lty=2,palette=rev.terrain.colors)
# black and white spectrogram 
spectro(pellu2,f=22050,ovlp=85,zp=16,
    palette=rev.gray.colors.1)
# colour modifications
spectro(pellu2,f=22050,ovlp=85,zp=16,
palette=rev.cm.colors,osc=TRUE,colwave="orchid1") 
spectro(pellu2,f=22050,ovlp=85,zp=16,osc=TRUE,palette=rev.heat.colors,
colbg="black",colgrid="white", colwave="white",colaxis="white",collab="white")
}

\keyword{dplot}
\keyword{ts}
