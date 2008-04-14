\name{deletew}
\alias{deletew}

\title{Delete a section of a time wave}

\description{
This function selects and delete a section of data describing a time wave.
Original section and section after deletion can be plotted as oscillograms for comparison. 
}
\usage{
deletew(wave, f, from = NULL, to = NULL, choose = FALSE, plot = FALSE,
marks = TRUE, Sample = FALSE,...)
}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{from}{start position (in s).}
  \item{to}{end position (in s).}
  \item{choose}{logical, if \code{TRUE} start (=\code{from}) and end (=\code{to})
		points can be graphically chosen with a cursor on the oscillogram.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of original
   and cut sections (by default \code{FALSE}).}
  \item{marks}{logical, if \code{TRUE} shows the start and end mark on the plot
(by default \code{TRUE}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.  
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
  }
  
\details{If \code{plot} is \code{TRUE} returns a two-frame plot with both 
original and resulting sections.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}},\code{\link{cutw}},
\code{\link{fadew}}, \code{\link{mutew}}, \code{\link{pastew}},
\code{\link{revw}}, \code{\link{zapsilw}}}

\examples{
# deletion a 0.4 s section in a bird song
data(tico)
a<-deletew(tico,f=22050,from=0.5,to=0.9)
oscillo(a,22050)
# a direct way to see what has been cut
deletew(tico,f=22050,from=0.5,to=0.9,plot=TRUE)
}

\keyword{dplot}
\keyword{ts}
