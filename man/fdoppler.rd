\name{fdoppler}

\alias{fdoppler}

\title{Doppler effect}

\description{
This function computes the altered frequency of a moving source due to the Doppler effect.
}

\usage{fdoppler(f, c = 340, vs, vo = 0, movs = "toward", movo = "toward")}

\arguments{
  \item{f}{original frequency produced by the source (in Hz or kHz)}
  \item{c}{speed of sound in meters/second.}
  \item{vs}{speed of the source in meters/second.}
  \item{vo}{speed of the observer in meters/second. The observer is static by default
    \emph{i.e.} \code{vo} = 0}
  \item{movs}{movement direction of the source in relation with observer position,
  either \code{"toward"} (by default) or \code{"away"}.}
  \item{movo}{movement direction of the observer in relation with the source position,
  either \code{"toward"} (by default, but be aware that
  the observer is static by default) or \code{"away"}.}
}

\details{
The altered frequency \emph{f'} is computed according to:\cr
\deqn{f{'} = f\times{\frac{c \pm v_{o}}{c \pm v_{s}}}}{%
      f' = f*(c+/-vo/(c+/-vs))}
with \emph{f} = original frequency produced by the source (in Hz or kHz),\cr
\emph{vs} = speed of the source,\cr
\emph{vo} = speed of the observer.
}

\value{The altered frequency is returned in a vector.}

\references{\url{http://www.kettering.edu/~drussell/Demos/doppler/doppler.html}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{You can use \code{\link{wasp}} to have exact values of \code{c}.
See examples.}

\seealso{\code{\link{wasp}}}

\examples{
# a 400 Hz source moving toward or away from the observer at 85 m/s
fdoppler(f=400,vs=85)
# [1] 533.3333
fdoppler(f=400,vs=85,movs="away")
# [1] 320
# use wasp() if you wish to have exact sound speed at a specific temperature
fdoppler(f=wasp(f=400,t=25)$c, vs=85)
# [1] 461.8667
# Doppler effect at different source speeds
f<-seq(1,10,by=1); lf<-length(f)
v<-seq(10,300,by=20); lv<-length(v)
res<-matrix(numeric(lf*lv),ncol=lv)
for(i in 1:lv) res[,i]<-fdoppler(f=f,vs=v[i])
op<-par(bg="lightgrey")
matplot(x=f,y=res,type="l",lty=1,las=1,col= spectro.colors(lv),
xlab="Source frequency (kHz)", ylab="Altered frequency (kHz)")
legend("topleft",legend=paste(as.character(v),"m/s"),
lty=1,col= spectro.colors(lv))
title(main="Doppler effect at different source speeds")
par(op)
}

\keyword{math}
