\name{deletew}
\alias{deletew}

\title{Delete a section of a time wave}

\description{
This function selects and delete a section of data describing a time wave.
Original section and section after deletion can be plotted as oscillograms for comparison. 
}
\usage{
deletew(wave, f, from = FALSE, to = FALSE, plot = FALSE, marks = TRUE, ...)
}

\arguments{
  \item{wave}{data describing a time wave.}
  \item{f}{sampling frequency of \code{wave}.}
  \item{from}{start mark (in s).}
  \item{to}{end mark (in s).}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of original
   and cut sections (by default \code{FALSE}).}
  \item{marks}{logical, if \code{TRUE} shows the start and end mark on the plot (by default \code{TRUE}).}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
  }
  
\details{If \code{plot} is \code{TRUE} returns a two-frame plot with both 
original and resulting sections.}

\value{
If \code{plot} is \code{FALSE}, \code{deletew} returns a one-column matrix
describing the resulting wave with the same sampling frequency as original data.
}

\author{Jérôme Sueur \email{jerome.sueur@univ-tours.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{cutw}}, \code{\link{pastew}},
\code{\link{timer}}.}

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
