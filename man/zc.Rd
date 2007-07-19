\name{zc}

\alias{zc}

\title{Instantaneous frequency of a time wave by zero-crossing}

\description{
This function measures the period of a full oscillating cycle.
}

\usage{
zc(wave, f, plot = TRUE, interpol = 1, threshold = FALSE,
xlab = "Time (s)", ylab = "Frequency (kHz)", ylim = c(0, f/2000),...)
}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{plot}{logical, if \code{TRUE} plots the dominant frequency along
  the time wave(by default \code{TRUE}).}
  \item{interpol}{interpolation factor.}
  \item{threshold}{amplitude threshold for signal detection (in \% ).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{ylim}{the range of y values.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
If \code{plot} is \code{FALSE}, \code{zc} returns a vector of numeric data with
the instantaneous frequency. 
}

\value{
If \code{plot} is \code{FALSE}, \code{zc} functions returns a vector of numeric data.
Values corresponds to the instantaneous frequency of the time wave. \sQuote{NA}s 
correspond either to pause periods (e. g. detected applying \code{threshold} or
sections of the time wave not corrsing the zero line. To remove \sQuote{NA}s with 
\code{\link{na.omit}} allows to get only instantaneous frequency values but
discards information about pause sections.}

\note{\code{interpol} adds points to the time wave by linear interpolation
(through \code{\link{approx}}). This increases measurement precision but as well
time process. Type argument of \code{plot} cannot be set to \dQuote{l}.} 

\references{Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998.
\emph{Animal acoustic communication.} Springer, Berlin, Heidelberg.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr},
Caroline Simonis-Sueur \email{csimonis@mnhn.fr} and 
Thierry Aubin \email{thierry.aubin@u-psud.fr}}

\examples{
data(pellucens)
pellu1<-cutw(pellucens,f=22050,from=0,to=1,plot=FALSE)
# without interpolation
zc(pellu1,f=22050,threshold=5,pch=20)
# with interpolation
zc(pellu1,f=22050,threshold=5,interpol=20,pch=20)
# a way to plot with a line and to filter low frequencies
pellu2<-zc(pellu1,f=22050,threshold=5,interpol=20,plot=FALSE)
pellu3<-na.omit(pellu2)
pellu4<-pellu3[pellu3>3]
plot(x=seq(0,nrow(pellu1)/22050,length.out=length(pellu4)),
    y=pellu4,type="l",xlab="Time(s)",ylab="Frequency(kHz)")
}

\keyword{dplot}
\keyword{ts}
