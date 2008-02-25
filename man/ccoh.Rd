\name{ccoh}

\alias{ccoh}

\title{Continuous coherence function between two time waves}

\description{
This function returns a two-dimension coherence representation between two time
waves. The function corresponds to a sliding coherence function along the 
two signals. This produces a 2-D density plot.
An amplitude contour plot can be overlaid.}

\usage{
ccoh(wave1, wave2, f, wl = 512, ovlp = 0, plot = TRUE,
grid = TRUE, scale = TRUE, cont = FALSE,
collevels = seq(0, 1, 0.01), palette = rev.heat.colors,
contlevels = seq(0, 1, 0.01), colcont = "black",
colbg="white", colgrid = "black",
colaxis = "black", collab="black",
plot.title = title(main = "", xlab = "Time (s)", 
ylab = "Frequency (kHz)"), scalelab = "Coherence",
scalefontlab = 1, scalecexlab =0.75, axisX = TRUE, axisY = TRUE,
flim = NULL, flimd = NULL,
...)
}

\arguments{
  \item{wave1}{data or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}) describing a first time wave.}
  \item{wave2}{data or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}) describing a second time wave.}
  \item{f}{sampling frequency of \code{wave1} and \code{wave1} (in Hz).
          Does not need to be specified if \code{wave1} and/or \code{wave2} are/is (a)
          \code{\link[sound]{Sample}} object(s).}
  \item{wl}{window length for the analysis
    (even number of points, by default = 512).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the continuous coherence function
        (by default \code{TRUE}).}
  \item{grid}{logical, if \code{TRUE} plots a y-axis grid
        (by default \code{TRUE}).}
  \item{scale}{logical, if \code{TRUE} plots a dB colour scale on the right
      side of the plot (by default \code{TRUE}).}
  \item{cont}{logical, if \code{TRUE} overplots contour lines on the plot
      (by default \code{FALSE}).}
  \item{collevels}{a set of levels which are used to partition
      the amplitude range of the coherence (should be between 0 and 1.}
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
  \item{flim}{modifications of the frequency Y-axis limits.} 
  \item{flimd}{dynamic modifications of the frequency Y-axis limits. New \code{wl}
      and \code{ovlp} arguments are applied to increase time/frequency resolution.}    
  \item{\dots}{other \code{\link{contour}} and \code{\link{oscillo}}
    graphical parameters.}
}

\details{Coherence is a frequency domain function computed to show the degree
of a relationship between two signals.
The value of the coherence function ranges between zero and one,
where a value of zero indicates there is no causal relationship between the signals.
A value of one indicates the existence of linear frequency response between the
two signals. This can be used, for instance, to compare the input and output
signals of a system.\cr
Any colour palette can be used. In particular, it is possible to use other
palettes coming with \pkg{seewave}: \code{rev.gray.colors.1},
\code{rev.gray.colors.2}, 
\code{spectro.colors}, \code{rev.terrain.colors}, \code{rev.topo.colors},
\code{rev.cm.colors} corresponding to the reverse of  
\code{terrain.colors}, \code{topo.colors}, \code{cm.colors}.\cr
Use \code{\link{locator}} to identify  points.
}

\value{
If \code{plot} is \code{FALSE}, this function returns a matrix. Each column
corresponds to a coherence function of length \code{wl}.}

\author{Jérôme Sueur \email{sueur@mnhn.fr} but this function is
mainly based on \code{\link{spec.pgram}} by Martyn Plummer, Adrian Trapletti
and B.D. Ripley}

\note{This function is based on \code{\link{spec.pgram}}, \code{\link{contour}} and
\code{\link{filled.contour}}. See \code{\link{spectro}} for graphical changes.}

\seealso{\code{\link{coh}}, \code{\link{spectro}}, \code{\link{spec.pgram}}.}

\examples{
wave1<-synth(d=1,f=4000,cf=500)
wave2<-synth(d=1,f=4000,cf=800)
ccoh(wave1,wave2,f=4000)
}

\keyword{dplot}
\keyword{ts}
