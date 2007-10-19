\name{oscillo}

\alias{oscillo}

\title{Show a time wave as an oscillogram}

\description{
This graphical function displays a time wave as an oscillogram in a single or 
multi-frame plot. The envelope of the wave can also be shown.
}

\usage{
oscillo(wave, f, from = FALSE, to = FALSE, zoom = FALSE,
k=1, j=1, labels = TRUE, byrow = TRUE, env = FALSE, smooth = NULL,
ksmooth = NULL, identify = FALSE, plot = TRUE, colwave = "black", 
colbg = "white", coltitle = "black", collab = "black", 
cexlab = 1, fontlab = 1, colline = "black", 
colaxis = "black", coly0 = "grey47",
title = FALSE, xaxt="s", yaxt="n", bty = "l")
}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{from}{start of the oscillogram (in s).}
  \item{to}{end of the oscillogram (in s).}
  \item{zoom}{time zoom in with start and end points chosen on the oscillogram with a cursor.}  
  \item{k}{number of horizontal sections (by default =1).}
  \item{j}{number of vertical sections (by default =1).}
  \item{labels}{if \code{TRUE} plots time and amplitude labels (by default \code{TRUE}).}  
  \item{byrow}{logical, if \code{TRUE}, the sections are filled by
          rows, otherwise the sections are filled by colmuns (by default \code{TRUE}).}
  \item{env}{logical, if \code{TRUE} plots the envelope of the time wave instead 
      (by default \code{FALSE}).}
  \item{smooth}{if \code{env} is \code{TRUE} then smoothes the envelope 
      by floating average.}
  \item{ksmooth}{if \code{env} is \code{TRUE} then smoothes the envelope 
      by applying a \code{\link{kernel}} smoothing.}
  \item{identify}{returns the time coordinate of points chosen  with a cursor
      on the oscillogram.}
   \item{plot}{logical, if \code{TRUE} returns an oscillographic or envelope plot
      of \code{wave}(by default \code{TRUE}).}
  \item{colwave}{colour of the oscillogram or of the envelope.}
  \item{colbg}{background colour.}
  \item{coltitle}{if \code{title} is \code{TRUE}, colour of the title.}
  \item{cexlab}{character size for axes title.}
  \item{fontlab}{font for axes title.}
  \item{collab}{colour of axes title.}  
  \item{colline}{colour of axes line.}
  \item{colaxis}{colour of the axes.}
  \item{coly0}{colour of the y=0 line.}
  \item{title}{logical, if \code{TRUE} plots the title 
      with information on time and \code{f} (by default \code{FALSE}).}
  \item{xaxt}{equivalent to \code{xaxt} of \code{\link{par}} (by default =\code{"s"}).}    
  \item{yaxt}{equivalent to \code{yaxt} of \code{\link{par}} (by default =\code{"n"}).}    
  \item{bty}{the type of box to be drawn around the oscillogram.} 
}

\value{
Data are returned as one-column matrix if \code{plot} is \code{FALSE}.
\code{identify} returns a numeric object with the time coordinate 
of points successively chosen on the oscillogram. 
}


\note{
\code{zoom} is similar to but more visual than \code{from} and/or \code{to}.
\code{zoom} and \code{identify} do work with a single-frame window only
(\emph{i. e.} with \code{k} = 1 and \code{j} = 1).\cr
Press \sQuote{Stop} button of the tools bar after choosing the appropriate
points on the oscillogram.\cr
} 


\author{Jérôme Sueur \email{sueur@mnhn.fr} and 
Caroline Simonis \email{csimonis@mnhn.fr}.}

\seealso{\code{\link{oscilloST}},\code{\link{cutw}}, \code{\link{pastew}},
\code{\link{timer}}}

\examples{
data(tico)
# a simple oscillogram of a bird song
oscillo(tico,f=22050)
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
# envelope representation without and with smoothing
oscillo(tico,f=22050,env=TRUE)
oscillo(tico,f=22050,env=TRUE,smooth=80)
oscillo(tico,f=22050,env=TRUE,ksmooth=kernel("daniell",10))
oscillo(tico,f=22050,env=TRUE,ksmooth=kernel("daniell",500))
# overplot of oscillographic and envelope representation
env<-oscillo(tico,f=22050,env=TRUE,plot=FALSE,smooth=40)
ticonorm<-tico/max(tico)
envnorm<-env/max(env)
oscillo(ticonorm,f=22050)
par(new=TRUE)
plot(envnorm,type="l",col="red",xaxs="i",yaxs="i",ann=FALSE,xaxt="n",yaxt="n",
    ylim=range(ticonorm),bty="l",lwd=2)
legend(x=4, y=1,"smoothed envelope", col="red",lty=1,lwd=2,bty="n",cex=0.75)
# full colour modifications in a two-frame oscillogram
oscillo(tico,f=22050,k=4,j=1,title=TRUE,colwave="black",colbg="grey",
    coltitle="yellow",collab="red",colline="white",
    colaxis="blue",coly0="grey50")
# change the title
data(orni)
oscillo(orni,f=22050,title=FALSE)
mtext("The song of a famous cicada",side=3,line=2,font=2,col="red",cex=1.5)
}

\keyword{dplot}
\keyword{ts}
