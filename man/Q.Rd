\name{Q}

\alias{Q}

\title{Resonance quality factor of a frequency spectrum}

\description{
This function estimates the frequency pureness of a time wave by returning the resonant
quality factor Q at a specific dB level.}

\usage{
Q(x, range, level = -3, plot = TRUE, colval = "red",
cexval = 1, fontval = 1, flab = "Frequency (kHz)",
alab = "Relative amplitude (dB)",...)
}

\arguments{
  \item{x}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}}, or \code{\link{meanspec}} (in dB).}
  \item{range}{range of \code{x} (in kHz).}
  \item{level}{frequency bandwidth set by an amplitude value relative
      to \code{spectrum} (in dB).}
  \item{plot}{logical, if \code{TRUE} returns the spectrum 
      with Q plotted (by default \code{TRUE}).}
  \item{colval}{colour of plotting Q.}
  \item{cexval}{character size of plotting Q.}
  \item{fontval}{font of plotting Q.}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the ampltiude axis.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{A high Q value indicates a highly resonant system.}

\value{Q is returned as a single numeric data.}

\author{Jérôme Sueur \email{jerome.sueur@ibaic.u-psud.fr}}

\note{This function is based on \code{\link{fft}}.} 

\seealso{\code{\link{spec}}, \code{\link{meanspec}}, \code{\link{corspec}},
\code{\link{fft}}.}

\examples{
# bird song
data(tico)
t<-spec(tico,f=22050,wl=512,at=1.1,plot=FALSE,dB=TRUE)
op<-par(mfrow=c(2,1),las=1)
Q(t,range=c(0,11.025),type="l")
Q(t,range=c(0,11.025),type="l",xlim=c(3.8,4.2),ylim=c(-60,0))
title("zoom in")
par(op)
# cricket, changing the dB level
data(pellucens)
p<-spec(pellucens,f=11025,wl=512,at=0.5,plot=FALSE,dB=TRUE)
op<-par(mfrow=c(3,1))
Q(p,range=c(0,11.25/2),type="l",xlim=c(1.8,2.6),ylim=c(-70,0))
title("level = - 3 (default value)",col.main="red")
Q(p,range=c(0,11.25/2),type="l",level=-6,
    xlim=c(1.8,2.6),ylim=c(-70,0),colval="blue")
title("level = - 6",col.main="blue")
Q(p,range=c(0,11.25/2),type="l",level=-9,
    xlim=c(1.8,2.6),ylim=c(-70,0),colval="green")
title("level = - 9",col.main="green")
par(op)
}

\keyword{dplot}
\keyword{ts}
