\name{cutw}

\alias{cutw}

\title{Cut a section of a time wave}

\description{
  This function selects and cuts a section of data describing a time wave.
  Original and cut sections can be plotted as oscillograms for comparison. 
}

\usage{
cutw(wave, f, from = NULL, to = NULL, choose = FALSE,
plot = FALSE, marks = TRUE, output="matrix", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{from}{start mark (in s).}
  \item{to}{end mark (in s).}
  \item{choose}{logical, if \code{TRUE} start (=\code{from}) and end (=\code{to})
    points can be graphically chosen with a cursor on the oscillogram.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot of original
    and cut sections (by default \code{FALSE}).}
  \item{marks}{logical, if \code{TRUE} shows the start and end mark on the plot
    (by default \code{TRUE}).}
  \item{output}{character string, the class of the object to return,
    either \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{If \code{plot} is \code{TRUE} returns a two-frame plot with both 
  original and cut sections.}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}},\code{\link{deletew}},
  \code{\link{fadew}},\code{\link{mutew}},\code{\link{pastew}},\code{\link{revw}},
  \code{\link{zapsilw}}}

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
