\name{dBscale}

\alias{dBscale}

\title{dB colour scale for a spectrogram display}

\description{
This function displays a vertical or horizontal dB colour scale to be used with
\code{\link{spectro}} plots.
}

\usage{
dBscale(collevels, palette = spectro.colors, side = 4,
textlab = "Amplitude\n(dB)", cexlab = 0.75,
fontlab = 1, collab = "black", colaxis = "black",...)
}


\arguments{

  \item{collevels}{a set of levels which are used to partition
      the amplitude range of the spectrogram (in dB).}
  \item{palette}{a color palette function to be used to assign colors in
      the plot, see \code{note}.}
  \item{side}{side of the axis.}
  \item{textlab}{text of the label.}  
  \item{cexlab}{character size of the label.}
  \item{fontlab}{font of the label.}
  \item{collab}{colour of the label.}
  \item{colaxis}{colour of the axis.}
  \item{\dots}{other \code{\link{axis}} arguments.}
}    

\author{Jerome Sueur \email{sueur@mnhn.fr} and 
Caroline Simonis \email{csimonis@mnhn.fr}.}

\note{This function, based on \code{\link{filled.contour}} by Ross Ihaka,
is not supposed to be used by itself but as a legend of \code{\link{spectro}}.\cr
Any colour palette can be used. In particular, it is possible
to use other palettes coming with
\pkg{seewave}: \code{rev.gray.colors.1}, \code{rev.gray.colors.2},
\code{rev.heat.colors}, \code{rev.terrain.colors}, \code{rev.topo.colors}, 
\code{rev.cm.colors} corresponding to the reverse of \code{heat.colors}, 
\code{terrain.colors}, \code{topo.colors}, \code{cm.colors}.}

\seealso{\code{\link{spectro}}.}

\examples{
data(pellucens)
# place the scale on the left and not on the right as spectro() does
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1, 2), nc = 2), widths = c(1, 5))
par(mar=c(5,3,4,2))
dBscale(collevels=seq(-30,0,1),side=2)
par(mar=c(5,4,4,2))
spectro(pellucens, f=22050,wl=512,scale=FALSE)
par(def.par)
# place the scale on the top and not on the right as spectro() does
def.par <- par(no.readonly = TRUE)
layout(matrix(c(0,1,2,2), nc = 2, byrow=TRUE),widths=c(1,2),heights=(c(1,5.5)))
par(mar=c(0.5,3,4,2))
dBscale(collevels=seq(-30,0,1), textlab = "",side=3)
mtext("Amplitude (dB)",side=2,line = 1,at=0.6,cex=0.75)
par(mar=c(5,4,0.5,2))
spectro(pellucens, f=22050,wl=512,scale=FALSE)
par(def.par)   
}

\keyword{dplot}
\keyword{ts}
