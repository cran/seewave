\name{field}

\alias{field}

\title{Near field and far field limits}

\description{
This function helps in knowing whether you are working in 
the near or far field.
}

\usage{field(f, d)}

\arguments{
  \item{f}{frequency (Hz)}
  \item{d}{distance from the sound source (m)}
}

\details{
Areas very close to the sound source are in the near-field where the contribution
of particle velocity to sound energy is greater thant that of sound pressure and where
these components are not in phase. Sound propagation properties are also different
near or far from the source. It is therefore important to know where the microphone
was from the source.\cr
To know this, the product k*d is computed according to:
\deqn{k\times{d} = \frac{f}{c}\times{d}}{% 
      k*d = (f/c)*d}
with \emph{d} = distance from the source (m), \emph{f} = frequency (Hz)
and \emph{c} = sound celerity (m/s).\cr
If k*d is greatly inferior 1 then the microphone is in the near field.\cr
The decision help returned by the function follows the rule:\cr
far field: \deqn{k\times{d} > 1}{%
k*d > 1} 
between near and far field limits: \deqn{0.1 \leq k\times{d} \leq 1}{%
0.1 <= k*d <= 1}
near field: \deqn{k\times{d} < 0.1}{%
k*d < 0.1}.
}

\value{
A list of two values is returned:
  \item{kd}{the numeric value k*d used to take a decision}
  \item{d}{a character string giving the help decision.}
}

\note{This function works for air-borne sound only.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\examples{
# 1 kHz near field at 1 cm from the source
field(f=1000,d=0.01)
# playing with distance from source and sound frequency
op<-par(bg="lightgrey")
D<-seq(0.01,0.5,by=0.01); nD<-length(D)
F<-seq(100,1000,by=25); nF<-length(F)
a<-matrix(numeric(nD*nF),nrow=nD)
for(i in 1:nF) a[,i]<-field(f=F[i],d=D)$kd
matplot(x=D,y=a,type="l",lty=1,col= spectro.colors(nF),
  xlab="Distance from the source (m)", ylab="k*d")
title("Variation of the product k*d with distance and frequency")
text(x=c(0.4,0.15),y=c(0.02,1), c("Near Field","Far Field"),font=2)
legend(x=0.05,y=1.4,c("100 Hz","1000 Hz"),lty=1,
  col=c(spectro.colors(nF)[1],spectro.colors(nF)[nF]),bg="grey")
abline(h=0.1)
par(op)
}

\keyword{ts}
