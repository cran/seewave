\name{dBweight}

\Rdversion{1.1}

\alias{dBweight}

\title{dB weightings}

\description{This function returns the four most common dB weightings.}

\usage{dBweight(f, dBref = NULL)}

\arguments{
  \item{f}{frequency (in Hz).}
  \item{dBref}{dB reference level (by default \code{NULL}).}
}

\details{By default, the function returns four weightings. When
  \code{dBref} is not \code{NULL} then the function returns the
  conversion from a dB reference level to four dB weighting levels.}

\value{\code{dBweight} returns a list of four items corresponding to four
  dB weightings.
\item{A}{dB (A)}
\item{B}{dB (B)}
\item{C}{dB (C)}
\item{D}{dB (D)}
}

\references{\url{http://en.wikipedia.org/wiki/A-weighting}.}

\author{Jerome Sueur \email{sueur@mnhn.fr} and Zev Ross}

\note{The transfer equations used here come from Wipipedia but they were originally coming
  from the appendix of an international standard on the design
  performance of sound level meters IEC 651:1979 (Neil Glenister, pers. com.).
}


\seealso{\code{\link{convSPL}}, \code{\link{moredB}}}

\examples{
# weight for a 50 Hz frequency
dBweight(f=50)
# A weight for the 1/3 Octave centre frequencies.
dBweight(f=c(20,25,31.5,40,50,63,80,100,125,160,200,250,
315,400,500,630,800,1000,1500,
1600,2000,2500,3150,4000,5000,
6300,8000,10000,12500,16000,20000))$A
# correction for a 50 Hz sound emitted at 100 dB
dBweight(f=50, dB=100)
# weighting curves plot
f <- seq(10,20000,by=10)
par(las=1)
plot(f, dBweight(f)$A, type="n", log="x",
xlim=c(10,10^5),ylim=c(-80,20),xlab="",ylab="",xaxt="n",yaxt="n")
abline(v=c(seq(10,100,by=10),seq(100,1000,by=100),
seq(1000,10000,by=1000),seq(10000,100000,by=10000),
c(100,1000,10000,100000)),col="lightgrey",lty=2)
abline(v=c(100,1000,10000,100000),col="grey")
abline(h=seq(-80, 20, 20),col="grey")
par(new=TRUE)
plot(f, dBweight(f)$A, type="l", log="x",
xlab="Frequency (Hz)", ylab="dB",lwd=2, col="blue", xlim=c(10,10^5),ylim=c(-80,20))
title(main="Acoustic weighting curves (10 Hz - 20 kHz)")
lines(x=f, y=dBweight(f)$B, col="green",lwd=2)
lines(x=f, y=dBweight(f)$C, col="red",lwd=2)
lines(x=f, y=dBweight(f)$D, col="black",lwd=2)
legend("bottomright",legend=c("dB(A)","dB(B)","dB(C)","dB(D)"),
lwd=2,col=c("blue","green","red","black"),bty="o",bg="white")
}

\keyword{ts}
\keyword{math}
