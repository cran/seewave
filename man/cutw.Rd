\name{cutw}

\alias{cutw}

\title{Cut a section of a time wave}

\description{
This function selects and cuts a section of data describing a time wave.
Original and cut sections can be plotted as oscillograms for comparison. 
}

\usage{
cutw(wave, f, from, to, plot = FALSE, marks = TRUE, ...)
}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{from}{start mark (in s).}
  \item{to}{end mark (in s).}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of original
   and cut sections (by default \code{FALSE}).}
  \item{marks}{logical, if \code{TRUE} shows the start and end mark on the plot (by default \code{TRUE}).}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
  }

\details{If \code{plot} is \code{TRUE} returns a two-frame plot with both 
original and cut sections.}

\value{
If \code{plot} is \code{FALSE}, \code{cutw} returns a one-column matrix
describing the cut wave with the same sampling frequency as the original wave.
}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{pastew}}, \code{\link{deletew}},
\code{\link{mute}}, \code{\link{timer}}}

\examples{
# a 0.4 s section in a bird song
data(tico)
a<-cutw(tico,f=22050,from=0.5,to=0.9)
oscillo(a,22050)
# a direct way to see what has been cut
cutw(tico,f=22050,from=0.5,to=0.9,plot=TRUE)
}

\keyword{dplot}
\keyword{ts}
