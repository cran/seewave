\name{Q}

\alias{Q}

\title{Resonance quality factor of a frequency spectrum}

\description{
  This function estimates the frequency pureness of a time wave by returning the resonant
  quality factor Q at a specific dB level.}

\usage{
Q(spec, f = NULL, level = -3, mel = FALSE, plot = TRUE, colval = "red",
cexval = 1, fontval = 1, flab = NULL,
alab = "Relative amplitude (dB)", type = "l", ...)
}

\arguments{
  \item{spec}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}}, or \code{\link{meanspec}} (in dB).
    This can be either a two-column matrix (col1 = frequency, col2 = amplitude)
    or a vector (amplitude).}
  \item{f}{sampling frequency of the wave used to obtain \code{spec}
    (in Hz). Not necessary if \code{spec} is a two columns matrix obtained
    with \code{\link{spec}} or \code{\link{meanspec}}.}
  \item{level}{frequency bandwidth set by an amplitude value relative
    to \code{spectrum} (in dB).}
  \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}
  \item{plot}{logical, if \code{TRUE} returns the spectrum 
    with Q plotted (by default \code{TRUE}).}
  \item{colval}{colour of plotting Q.}
  \item{cexval}{character size of plotting Q.}
  \item{fontval}{font of plotting Q.}
  \item{flab}{title of the frequency axis.}
  \item{alab}{title of the amplitude axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{A high Q value indicates a highly resonant system.}

\value{A list is returned with the following four items:
  \item{Q}{a numeric vector of length 1 returning the Q factor (no units)}
  \item{dfreq}{a numeric vector of length 1 the dominant frequency (kHz)}
  \item{fmin}{a numeric vector of length 1 returning the minimum
  frequency of the -dB level bandwidth (kHz)} 
  \item{fmax}{a numeric vector of length 1 returning the minimum
  frequency of the -dB level bandwidth (kHz)}
  \item{bwd}{a numeric vector of length 1 returning the bandwidth,
  i. e. \code{fmax-fmin} (kHz)}
 }

\author{Jerome Sueur, improved by Laurent Lellouch}

\note{This function is based on \code{\link{fft}}.} 

\seealso{\code{\link{spec}}, \code{\link{meanspec}}, \code{\link{corspec}},
  \code{\link{fft}}.}

\examples{
# bird song
data(tico)
t<-spec(tico,f=22050,at=1.1,plot=FALSE,dB="max0")
op<-par(mfrow=c(2,1),las=1)
Q(t,type="l")
Q(t,type="l",xlim=c(3.8,4.2),ylim=c(-60,0))
title("zoom in")
par(op)
# cricket, changing the dB level
data(pellucens)
p<-spec(pellucens,f=11025,at=0.5,plot=FALSE,dB="max0")
op<-par(mfrow=c(3,1))
Q(p,type="l",xlim=c(1.8,2.6),ylim=c(-70,0))
title("level = - 3 (default value)",col.main="red")
Q(p,type="l",level=-6,
    xlim=c(1.8,2.6),ylim=c(-70,0),colval="blue")
title("level = - 6",col.main="blue")
Q(p,type="l",level=-9,
    xlim=c(1.8,2.6),ylim=c(-70,0),colval="green")
title("level = - 9",col.main="green")
par(op)
}

\keyword{dplot}
\keyword{ts}
