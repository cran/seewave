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
coly01 = "grey47", coly02 = "black", title = FALSE,
bty = "l")
}

\arguments{
	\item{wave1}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
	\item{wave2}{a \code{vector}, a \code{matrix} (second column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (right channel),
	or \code{\link[tuneR]{Wave}} (right channel).}
  \item{f}{sampling frequency of \code{wave1} and \code{wave2} (in Hz).
  Does not need to be specified if \code{wave1} and/or \code{wave2} are/is
	of class \code{ts}, \code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
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

\author{Jerome Sueur \email{sueur@mnhn.fr} and 
Caroline Simonis \email{csimonis@mnhn.fr}.}

\seealso{\code{\link{oscillo}}}

\examples{
a<-synth(f=8000,d=1,cf=2000,am=c(50,10),plot=FALSE)
b<-synth(f=8000,d=1,cf=1000,fm=c(0,0,2000),plot=FALSE)
oscilloST(a,b,f=8000)
}

\keyword{dplot}
\keyword{ts}
