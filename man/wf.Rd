\name{wf}

\alias{wf}

\title{Waterfall display}

\description{This function returns a waterfall display of
  a short-term Fourier transform or of any matrix.}

\usage{
wf(wave, f = NULL, wl = 512, zp = 0, ovlp = 0, fftw= FALSE, dB = "max0",
dBref = NULL, wn = "hanning", x = NULL,
hoff = 1, voff = 1, col = heat.colors,
xlab = "Frequency (kHz)", ylab = "Amplitude (dB)",
xaxis = TRUE, yaxis = TRUE,
density = NULL, border = NULL, lines = FALSE, lwd=NULL, ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{window length for the analysis
    (even number of points). (by default = 512)}
  \item{zp}{zero-padding (even number of points), see \code{Details}.}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
  library \code{fftw}. See Notes of the \code{spectro}.}
  \item{dB}{a character string specifying the type dB to return: "max0" for a
    maximum dB value at 0, "A", "B", "C" and "D" for common dB weights.}
  \item{dBref}{a dB reference value when \code{dB} is \code{TRUE}. \code{NULL} by default
    but should be set to 2*10e-5 for a 20 microPa reference.}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{x}{a matrix if \code{wave} is not provided.}
  \item{hoff}{horizontal 'offset' which shifts actual x-values slightly
    per row for visibility. Fractional parts will be removed.}
  \item{voff}{vertical 'offset' which separates traces.}
  \item{col}{a color or a color palette function to be used to assign colors in
    the plot}
  \item{xlab}{title of the frequency x-axis.}
  \item{ylab}{title of the amplitude y-axis.}
  \item{xaxis}{a logical, if \code{TRUE} adds the frequency x-axis
    according to \code{f}.}
  \item{yaxis}{a logical, if \code{TRUE} adds the amplitude y-axis
    according.}
  \item{density}{argument of \code{\link[graphics]{polygon}}:
    the density of shading lines, in lines per inch.  The default
    value of 'NULL' means that no shading lines are drawn. A zero
    value of 'density' means no shading nor filling whereas
    negative values (and 'NA') suppress shading (and so allow
    color filling).}
  \item{border}{argument of \code{\link[graphics]{polygon}}:
    the color to draw the border.  The default, 'NULL', means to
    use 'par("fg")'.  Use 'border = NA' to omit borders.}
  \item{lines}{a logical, if \code{TRUE} plots lines instead of surfaces
    (polygons).}
  \item{lwd}{line width.}
  \item{\dots}{other graphical arguments to passed to \code{\link{plot}}}
}

\details{Data input can be either a time wave (\code{wave}) or a matrix
  (\code{x}). In that case, if \code{xaxis} is set to \code{TRUE} the
  x-axis will follow the row index. To change it, turn \code{xaxis} to
  \code{FALSE} and use \code{\link[graphics]{axis}} afterwards. See examples.}

\note{The function is well adapted to display short-term Fourier
  transform. However, any matrix can be called using the argument
  \code{x} instead of \code{wave}.} 

\author{Carl G. Witthoft and Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{spectro}}, \code{\link{spectro3D}}, \code{\link{dynspec}}}

\examples{
data(tico)
wf(tico,f=22050)
# changing the display parameters
jet.colors <- colorRampPalette(c("blue", "green"))
wf(tico,f=22050, hoff=0, voff=2, col=jet.colors, border = NA)
# matrix input instead of a time wave and transparent lines display
m <- numeric()
for(i in seq(-pi,pi,len=40)) {m <- cbind(m,10*(sin(seq(0,2*pi,len=100)+i)))}
wf(x=m, lines=TRUE, col="#0000FF50",xlab="Time", ylab="Amplitude",
main="waterfall display")
}

\keyword{dplot}
\keyword{ts}