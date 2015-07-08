\name{oscillo}

\alias{oscillo}

\title{Show a time wave as an oscillogram}

\description{
  This graphical function displays a time wave as an oscillogram in a single or 
  multi-frame plot. The envelope of the wave can also be shown.
}

\usage{
oscillo(wave, f, from = NULL, to = NULL, scroll = NULL,
zoom = FALSE, k=1, j=1, cex,
labels = TRUE,  tlab = "Time (s)", alab = "Amplitude",            
byrow = TRUE, identify = FALSE, nidentify = NULL, 
plot = TRUE, colwave = "black", 
coltitle = "black", cextitle = 1.2, fonttitle = 2,
collab = "black", cexlab = 1, fontlab = 1,
colline = "black",
colaxis = "black", cexaxis = 1, fontaxis = 1,
coly0 = "lightgrey",
tcl = 0.5, title = FALSE, xaxt="s", yaxt="n", type="l", bty = "l")
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{from}{start of the oscillogram (in s).}
  \item{to}{end of the oscillogram (in s).}
  \item{scroll}{a numeric of length 1 allowing to move along the time wave using a slider panel.
    This numeric corresponds to the number of successive windows dividing the time wave.}
  \item{zoom}{time zoom in with start and end points chosen on the oscillogram with a cursor.}  
  \item{k}{number of horizontal sections (by default =1).}
  \item{j}{number of vertical sections (by default =1).}
  \item{cex}{Pitch size if \code{type = "p"}.}
  \item{labels}{if \code{TRUE} plots time and amplitude labels (by default \code{TRUE}).}  
  \item{tlab}{Label of time axis.}
  \item{alab}{Label of amplitude axis.}
  \item{byrow}{logical, if \code{TRUE}, the sections are filled by
    rows, otherwise the sections are filled by colmuns (by default \code{TRUE}).}
  \item{identify}{returns the time and amplitude coordinates of points chosen  with a cursor
    on the oscillogram.}
  \item{nidentify}{a numeric vector of length 1, specifies the number of
  points to identified on \code{wave} if \code{identify} is \code{TRUE}.}  
  \item{plot}{logical, if \code{TRUE} returns an oscillographic or envelope plot
    of \code{wave}(by default \code{TRUE}).}
  \item{colwave}{colour of the oscillogram or of the envelope.}
  \item{coltitle}{if \code{title} is \code{TRUE}, colour of the title.}
  \item{cextitle}{character size for the title.}
  \item{fonttitle}{font for the title.}
  \item{cexlab}{character size for axes labels.}
  \item{fontlab}{font for axes labels.}
  \item{collab}{colour of axes labels.}  
  \item{colline}{colour of axes line.}
  \item{colaxis}{colour of the axis annotation.}
  \item{fontaxis}{font of axis annotation.}
  \item{cexaxis}{magnification for axis annotation.}
  \item{coly0}{colour of the y=0 line.}
  \item{tcl}{length of tick marks.}
  \item{title}{\code{TRUE} to add a title
    with information on \code{wave} duration and \code{f},
    \code{FALSE} to live it blank,
    or a character string to add any desired title.}
  \item{xaxt}{equivalent to \code{xaxt} of \code{\link{par}} (by default =\code{"s"}).}    
  \item{yaxt}{equivalent to \code{yaxt} of \code{\link{par}} (by default
    =\code{"n"}).}  
  \item{type}{type of plot, by default \code{"l"}. Use \code{"n"} for no
  plot.}
  \item{bty}{the type of box to be drawn around the oscillogram.} 
}

\value{
  Data are returned as one-column matrix if \code{plot} is \code{FALSE}.
  \code{identify} returns a two-column matrix with the time and
  amplitude coordinates of points successively chosen on the oscillogram. 
}


\note{
  \code{zoom} is similar to but more visual than \code{from} and/or \code{to}.
  \code{zoom} and \code{identify} do work with a single-frame window only
  (\emph{i. e.} with \code{k} = 1 and \code{j} = 1).\cr
  Press \sQuote{Stop} button of the tools bar after choosing the appropriate
  points on the oscillogram.\cr
} 


\author{Jerome Sueur \email{sueur@mnhn.fr} and 
  Caroline Simonis \email{csimonis@mnhn.fr}.}

\seealso{\code{\link{dynoscillo}}, \code{\link{oscilloST}}, \code{\link{cutw}}, \code{\link{pastew}},
  \code{\link{timer}}}

\examples{
data(tico)
# a simple oscillogram of a bird song
oscillo(tico)
# zoom in
op<-par(mfrow=c(4,1),mar=c(4.5,4,2,2))
oscillo(tico,22050,cexlab=0.75)
oscillo(tico,22050,from=0.5,to=0.9,cexlab=0.75)
oscillo(tico,22050,from=0.65,to=0.75,cexlab=0.75)
oscillo(tico,22050,from=0.68,to=0.70,cexlab=0.75)
par(op)
# the same divided in four lines
oscillo(tico,f=22050,k=4,j=1)
# the same divided in different numbers of lines and columns
oscillo(tico,f=22050,k=4,j=4)
oscillo(tico,f=22050,k=2,j=2,byrow=TRUE)
oscillo(tico,f=22050,k=2,j=2,byrow=FALSE)
# overplot of oscillographic and envelope representations
oscillo(tico,f=22050)
par(new=TRUE)
env(tico,f=22050,colwave=2)
# full colour modifications in a two-frame oscillogram
op<-par(bg="grey")
oscillo(tico,f=22050,k=4,j=1,title=TRUE,colwave="black",
    coltitle="yellow",collab="red",colline="white",
    colaxis="blue",coly0="grey50")
par(op)
# change the title
data(orni)
oscillo(orni,f=22050,title="The song of a famous cicada")
# move along the signal using scroll
\dontrun{
require(rpanel)
oscillo(tico,f=22050,scroll=8)}
}

\keyword{dplot}
\keyword{ts}
