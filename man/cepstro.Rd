\name{cepstro}

\alias{cepstro}

\title{2D-cepstrogram of a time wave}

\description{
This function returns a two-dimension cepstrographic representation of a time
wave. The function corresponds to a short-term cepstral transform.
An amplitude contour plot can be overlaid.}

\usage{
cepstro(wave, f, wl = 512, ovlp = 0, plot = TRUE, grid = TRUE,
scale = TRUE, cont = FALSE, collevels = seq(0, 1, 0.01),
palette = rev.heat.colors, contlevels = seq(0, 1, 0.01),
colcont = "black", colbg="white", colgrid = "black",
colaxis = "black", collab = "black",
plot.title = title(main = "", xlab = "Time (s)", ylab = "Quefrency (kHz)"),
scalelab = "Amplitude", scalefontlab = 1, scalecexlab = 0.75,
axisX = TRUE, axisY = TRUE, ...)
}

\arguments{
  \item{wave}{data describing a time wave or a \code{\link[sound]{Sample}}
  object generated loading a wav file with \code{\link[sound]{loadSample}}
  (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified
  if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the cepstrogram
        (by default \code{TRUE}).}
  \item{grid}{logical, if \code{TRUE} plots a y-axis grid
        (by default \code{TRUE}).}
  \item{scale}{logical, if \code{TRUE} plots a dB colour scale on the right
      side of the cesptrogram (by default \code{TRUE}).}
  \item{cont}{logical, if \code{TRUE} overplots contour lines on the cepstrogram
      (by default \code{FALSE}).}
  \item{collevels}{a set of levels which are used to partition
      the amplitude range of the cepstrogram (in dB).}
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
  \item{axisX}{if \code{TRUE} plots time X-axis (by default \code{TRUE}).}  
  \item{axisY}{if \code{TRUE} plots frequency Y-axis (by default \code{TRUE}).}  
  \item{\dots}{other \code{\link{contour}} graphical parameters.}
}

\details{
It is unfortunately not possible to turn the \emph{y}-axis to a frequency scale.\cr
See \code{\link{spectro}} for the use of the graphical arguments. 
}

\value{
When \code{plot} is \code{FALSE}, a matrix is returned with the successive
cepstral profiles computed along time.}

\references{Oppenheim, A.V. and Schafer, R.W. 2004. From frequency
to quefrency: a history of the cepstrum. \emph{Signal Processing Magazine IEEE},
21: 95-106.}

\author{Jerome Sueur \email{sueur@mnhn.fr}.}

\note{This function is based on \code{\link{ceps}}.} 

\seealso{\code{\link{ceps}}, \code{\link{fund}}, \code{\link{autoc}}}

\examples{
data(sheep)
cepstro(sheep,f=8000)
}

\keyword{dplot}
\keyword{ts}

