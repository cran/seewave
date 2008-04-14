\name{spec}

\alias{spec}

\title{Frequency spectrum of a time wave}

\description{
This function returns the frequency spectrum (\emph{i.e.} the relative amplitude of the
frequency content) of a time wave. Results can be obtained either as absolute 
or dB data.
}

\usage{
spec(wave, f, wl = 512, wn = "hanning", PSD = FALSE,
PMF = FALSE, dB = FALSE,
at = NULL, from = NULL, to = NULL, peaks = FALSE,
identify = FALSE, col = "black", cex = 1, colpeaks = "red",
cexpeaks = 1, fontpeaks = 1, plot = 1, flab = "Frequency (kHz)",
alab = "Amplitude", flim = c(0, f/2000),
alim = NULL, type="l",...)
}

\arguments{
  \item{wave}{data describing a time wave or a \code{\link[sound]{Sample}}
      object generated loading a wav file with \code{\link[sound]{loadSample}}
     (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
     Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{PSD}{if \code{TRUE} return Power Spectra Density,
      \emph{i. e.} the square of the spectra.}  
  \item{PMF}{if \code{TRUE} return Probability Mass Function,
      \emph{i. e.} the probability distribution of frequencie.}    
  \item{dB}{logical, if \code{TRUE} returns the spectrum in dB 
      (by default \code{FALSE}).}
  \item{at}{position where to compute the spectrum (in s).}
  \item{from}{start mark where  to compute the spectrum (in s).}
  \item{to}{end mark where to compute the spectrum (in s).}
  \item{peaks}{returns peaks value for a given span (see \code{details}).}
  \item{identify}{to identify frequency and amplitude values on the plot
      with the help of a cursor.}
  \item{col}{colour of the spectrum.}
  \item{cex}{pitch size of the spectrum.}
  \item{colpeaks}{colour of peaks value plotted on the spectrum.}
  \item{cexpeaks}{character size of peaks value plotted on the spectrum.}
  \item{fontpeaks}{font of peaks value plotted on the spectrum.}
  \item{plot}{if \code{1} returns frequency on x-axis, if \code{2} returns 
  frequency on y-axis, (by default \code{1}).}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{flim}{range of frequency axis.}
  \item{alim}{range of amplitude axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
      See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
If \code{at}, \code{from} or \code{to} are \code{FALSE} then \code{spec} 
computes the spectrum of the whole signal.\cr
\code{peaks} setting corresponds to \code{dimension} of \code{\link{embed}}.\cr 
When using \code{identify}, press \sQuote{stop} tools bar button to return values in the console.
}

\value{
This function returns  a two-column matrix,
the first column corresponding to the frequency axis, the second column
corresponding to the amplitude axis.\cr 
If \code{identify} is \code{TRUE}, \code{spec} returns a list with
two elements:\cr 
  \item{freq}{the frequency of the points chosen on the spectrum}\cr 
  \item{amp}{the relative amplitude of the points chosen on the spectrum}\cr 
If \code{peaks} is \code{TRUE}, \code{spec} returns a list with
two elements:\cr 
  \item{spec}{the spectrum computed} of length \code{wl}/2, \cr 
  \item{peaks }{the peaks value (in kHz).}
}

\author{Jerome Sueur \email{sueur@mnhn.fr} and 
Martin Maechler \email{maechler@stat.math.ethz.ch} for \code{peaks}.}

\note{This function is based on \code{\link{fft}}.} 

\seealso{\code{\link{meanspec}}, \code{\link{dynspec}}, 
\code{\link{corspec}}, \code{\link{fft}}.}

\examples{
data(tico)
# spectrum of the whole signal, in absolute or dB amplitude,
# horizontaly or vertically
op<-par(mfrow=c(2,2))
spec(tico,f=22050)
spec(tico,f=22050,col="red",plot=2)
spec(tico,f=22050,dB=TRUE,col="blue")
spec(tico,f=22050,dB=TRUE,col="green",plot=2)
par(op)
# spectra computed at specific locations with peak value
op<-par(mfrow=c(2,2))
spec(tico,f=22050,wl=512,at=0.2,peak=175)
title("Note A")
spec(tico,f=22050,wl=512,at=0.7,peak=175)
title("Note B")
spec(tico,f=22050,wl=512,at=1.1,peak=175)
title("Note C")
spec(tico,f=22050,wl=512,at=1.6,peak=165)
title("Note D")
par(op)
# an indirect way to compare spectra 
a<-spec(tico,f=22050,wl=512,at=0.2,plot=FALSE)
b<-spec(tico,f=22050,wl=512,at=0.7,plot=FALSE)
c<-spec(tico,f=22050,wl=512,at=1.1,plot=FALSE)
d<-spec(tico,f=22050,wl=512,at=1.6,plot=FALSE)
all<-cbind(a[,2],b[,2],c[,2],d[,2])
matplot(x=a[,1],y=all,yaxt="n",
    xlab="Frequency (kHz)",ylab="Amplitude",xaxs="i",type="l")
legend(8,0.8,c("Note A","Note B", "Note C", "Note D"),bty="o",
    lty=c(1:4),col=c(1:4))
# spectrum from a particular position to another one
op<-par(mfrow=c(2,1))
oscillo(tico,f=22050)
abline(v=c(0.5,0.9),col="red",lty=2)
text(c("from 0.5 s","to 0.7 s"),
    x=c(0.5,0.9),y=rep(max(tico/1.1),2),col="red",pos=4)
spec(tico,f=22050,wl=512,from=0.5,to=0.9,col="red")
title("Spectrum of the note B")
par(op)
# spectrum and spectrogram
data(orni)
orni1<-cutw(orni,f=22050,from=0.32,to=0.39)
layout(matrix(c(1,2),nc=2),widths=c(3,1))
par(mar=c(5,4,3,0.5))
spectro(orni1,f=22050,wl=128,zp=8,ovlp=85,scale=FALSE)
par(mar=c(5,1,3,0.5))
spec(orni1,f=22050,col="red",plot=2,flab="",yaxt="n")
}
\keyword{dplot}
\keyword{ts}
