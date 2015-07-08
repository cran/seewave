
\name{spec}

\alias{spec}

\title{Frequency spectrum of a time wave}

\description{
  This function returns the frequency spectrum (\emph{i.e.} the relative amplitude of the
  frequency content) of a time wave. Results can be obtained either as absolute 
  or dB data.
}

\usage{
spec(wave, f, wl = 512, wn = "hanning", fftw = FALSE, norm = TRUE, PSD = FALSE,
PMF = FALSE, dB = NULL, dBref = NULL,
at = NULL, from = NULL, to = NULL,
identify = FALSE, col = "black", cex = 1,
plot = 1, flab = "Frequency (kHz)",
alab = "Amplitude", flim = NULL,
alim = NULL, type="l",...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{if \code{at} is not null, length of the window for the analysis
    (even number of points, by defaults = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default
    \code{"hanning"}).}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
    library \code{fftw} for faster computation. See Notes of the function \code{\link{spectro}}.}
  \item{norm}{if \code{TRUE} the spectrum is normalised (i. e. scaled) by
  its maximum.}
  \item{PSD}{if \code{TRUE} return Power Spectra Density,
    \emph{i. e.} the square of the spectra.}  
  \item{PMF}{if \code{TRUE} return Probability Mass Function,
    \emph{i. e.} the probability distribution of frequencie.}    
  \item{dB}{a character string specifying the type dB to return: "max0" for a
    maximum dB value at 0, "A", "B", "C" and "D" for common dB weights.}
  \item{dBref}{a dB reference value when \code{dB} is not \code{NULL}. \code{NULL} by default
    but should be set to 2*10e-5 for a 20 microPa reference (SPL).}
  \item{at}{position where to compute the spectrum (in s).}
  \item{from}{start mark where  to compute the spectrum (in s).}
  \item{to}{end mark where to compute the spectrum (in s).}
  \item{identify}{to identify frequency and amplitude values on the plot
    with the help of a cursor.}
  \item{col}{colour of the spectrum.}
  \item{cex}{pitch size of the spectrum.}
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
  computes the spectrum of the whole signal.}

\value{
  This function returns  a two-column matrix,
  the first column corresponding to the frequency axis, the second column
  corresponding to the amplitude axis.\cr 
  If \code{identify} is \code{TRUE}, \code{spec} returns a list with
  two elements:\cr 
  \item{freq}{the frequency of the points chosen on the spectrum}\cr 
  \item{amp}{the relative amplitude of the points chosen on the spectrum}\cr 
}

\author{Jerome Sueur}

\note{This function is based on \code{\link{fft}}.} 

\section{Warning}{The argument \code{peaks} is no more available
  (version > 1.5.6). See the function \code{\link{fpeaks}}
  for peak(s) detection.} 

\seealso{\code{\link{meanspec}}, \code{\link{fpeaks}},
  \code{\link{localpeaks}}, \code{\link{dynspec}}, 
  \code{\link{corspec}}, \code{\link{fft}}.}

\examples{
data(tico)
# spectrum of the whole signal, in absolute or dB amplitude,
# horizontaly or vertically
op<-par(mfrow=c(2,2))
spec(tico,f=22050)
spec(tico,f=22050,col="red",plot=2)
spec(tico,f=22050,dB="max0",col="blue")
spec(tico,f=22050,dB="max0",col="green",plot=2)
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
