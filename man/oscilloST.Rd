\name{oscilloST}

\alias{oscilloST}

\title{Show a stereo time wave as oscillograms}

\description{
  This graphical function displays a stereo (2 channels) time wave as an
  oscillogram in a two-frame plot. The envelope of the wave can also be shown.
}

\usage{
oscilloST(wave1, wave2 = NULL, f, from = NULL, to = NULL,
identify = FALSE, plot = TRUE, colwave1 = "black",
colwave2 = "blue", coltitle = "black",
collab = "black", cexlab = 1, fontlab = 1, colaxis = "black",
cexaxis = 1, coly01 = "grey47", coly02 = "black", title = FALSE,
bty = "l")
}

\arguments{
  \item{wave1}{a first R object.}     
  \item{wave2}{a second R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{from}{start of the oscillogram (in s).}
  \item{to}{end of the oscillogram (in s).}
  \item{identify}{returns the time coordinate of points chosen  with a cursor 
    on the bottom oscillogram.}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic or envelope plot
    of \code{wave}(by default \code{TRUE}).}
  \item{colwave1}{colour of the oscillogram or of the envelope of \code{wave1}.}
  \item{colwave2}{colour of the oscillogram or of the envelope of \code{wave2}.}
  \item{coltitle}{if \code{title} is \code{TRUE}, colour of the title.}
  \item{collab}{colour of axes title.}  
  \item{cexlab}{character size for axes title.}
  \item{fontlab}{font for axes title.}
  \item{colaxis}{colour of the axes}
  \item{cexaxis}{mangification for axes annotation.}
  \item{coly01}{colour of the y=0 line of \code{wave1}.}
  \item{coly02}{colour of the y=0 line of \code{wave1}.}
  \item{title}{logical, if \code{TRUE} plots the title 
    with information on time and \code{f} (by default \code{FALSE}).}
  \item{bty}{the type of box to be drawn around the oscillogram.} 
}

\value{
  Data are returned as two-column matrix if \code{plot} is \code{FALSE}.
  \code{identify} returns a numeric object with the time coordinate 
  of points successively chosen on the bottom oscillogram.}

\author{Jerome Sueur and Caroline Simonis.}

\seealso{\code{\link{oscillo}}, \code{\link{dynoscillo}}}

\examples{
a<-synth(f=8000,d=1,cf=2000,am=c(50,10),plot=FALSE)
b<-synth(f=8000,d=1,cf=1000,fm=c(0,0,2000),plot=FALSE)
oscilloST(a,b,f=8000)
}

\keyword{dplot}
\keyword{ts}
