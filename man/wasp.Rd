\name{wasp}

\alias{wasp}

\title{WAve length and SPeed of sound}

\description{
  This function returns the wavelength and the speed of sound
  of a given frequency in air, fresh-water or sea-water.}

\usage{wasp(f, t = 20, c = NULL, s = NULL, d = NULL, medium = "air")}

\arguments{
  \item{f}{frequency (Hz).}
  \item{t}{temperature (degree Celsius).}
  \item{c}{celerity (m/s) if a wavelength is to be found at a particular speed of sound.}
  \item{s}{salinity (parts per thousand) when \code{medium} is \code{"sea"}.}
  \item{d}{depth (m) when \code{medium} is \code{"sea"}.}
  \item{medium}{medium for sound propagation,
    either "air", "fresh" for fresh, or pure, water, "sea" for sea water.}
}

\details{
  Speed of sound in air is computed according to:\cr
  \deqn{c = 331.4 + 0.6\times{t}}{% 
    c = 331.4+0.6*t}

  Speed of sound in fresh-water is computed according to Marczak equation:\cr

  \deqn{c = 1.402385.10^{3} + 5.038813\times{t} - 5.799136.10^{-2}\times{t^{2}}}{% 
    c = 1.402385e3+5.038813*t-(5.799136e-2)*t^2}
  \deqn{+ 3.287156.10^{-4}\times{t^{3}} - 1.398845.10^{-6}\times{t^{4}}}{% 
    +(3.287156e-4)*t^3-(1.398845e-6)*t^4}
  \deqn{+ 2.787860.10^{-9}\times{t^{5}}}{% 
    +(2.787860e-9)*t^5}

  with \emph{t} = temperature in degrees Celsius;
  range of validity: 0-95 degrees Celcius at atmospheric pressure.\cr

  Speed of sound in sea-water is computed according to Mackenzie equation:\cr
  \deqn{c = 1448.96 + 4.591\times{t}- 5.304.10^{-2}\times{t^{2}}}{%
    c = 1448.96+4.591*t-(5.304e-2)*t^2}      
  \deqn{+ 2.374.10^{-4}\times{t^{3}} + 1.34\times{(s-35)} + 1.63.10^{-2}\times{d}}{%
    +(2.374e-4)*t^3+1.34*(s-35)+(1.63e-2)*d}  
  \deqn{+ 1.675.10^{-7}\times{d^{2}} - 1.025.10^{-2}\times{t}\times{(s-35)}}{%
    +(1.675e-7)*d^2-(1.025e-2)*t*(s-35)}      
  \deqn{- 7.139.10^{-13}\times{t}\times{d^3}}{% 
    -(7.139e-13)*t*d^3}

  with \emph{t} = temperature in degrees Celsius;
  \emph{s} = salinity in parts per thousand;
  \emph{d} = depth in meters;
  range of validity: temperature 2 to 30 degrees Celcius, salinity 25 to 40 parts per thousand, depth 0 to 8000 m.\cr

  Wavelength is obtained following:\cr
  \deqn{\lambda = \frac{c}{f}}{% 
    lambda = c/f}
  with \emph{c} = speed of sound in meters/second;
  \emph{f} = frequency in Hertz. 
}

\value{
  A list of two values is returned:
  \item{l}{wavelength in meters}
  \item{c}{speed of sound in meters/second.}
}

\references{\url{http://resource.npl.co.uk}}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\examples{
# wavelength (m) of a 2000 Hz air-borne sound at 20 degrees Celsius
wasp(f=2000)$l
# [1] 0.1717

# sound speed in sea at 0 and -500 m
# for a respective temperature of 22 degrees Celcius and 11 degrees Celcius
wasp(f=1000,s=30,d=c(0,500),t=c(22,11),medium="sea")$c
# [1] 1521.246 1495.414

# wavelength (m) of a 1000 Hz sound in a medium unspecified where c = 1497 m/s
wasp(f=1000,c=1497)$l
# [1] 1.497

# variation of wavelength according to frequency and air temperature
op<-par(bg="lightgrey")
a<-seq(1000,20000,by=100) ; na<-length(a)
b<-seq(-20,40,by=10) ; nb<-length(b)
res<-matrix(numeric(na*nb),nrow=na)
for(i in 1:nb) res[,i]<-wasp(a,t=b[i])$l
matplot(x=a,y=res,type="l",lty=1,col= spectro.colors(nb),
  xlab="Frequency (Hz)",ylab="Wavelength (m)")
title("Wavelength of air-borne sound at different temperatures (deg. C)")
legend(x=15000,y=0.3,c("-20","-10","0","10","20","30","40"),
  lty=1,col= spectro.colors(nb),bg="grey")
par(op)
}

\keyword{math}


