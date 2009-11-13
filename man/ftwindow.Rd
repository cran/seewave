\name{ftwindow}

\alias{ftwindow}

\title{Fourier transform windows}

\description{
  Generates different Fourier Transform windows.
}

\usage{
ftwindow(wl, wn = "hamming")
}

\arguments{
  \item{wl}{window length}
  \item{wn}{window name: \code{bartlett}, \code{blackman},
    \code{flattop}, \code{hamming}, \code{hanning}, or \code{rectangle}
    (by default \code{hamming}).}
}


\value{
  A vector of length \code{wl}.
}

\references{Harris, F.J., 1978. On the use of windows for harmonic analysis
  with the discrete Fourier Transform. \emph{Proceedings of the IEEE}, 66(1): 51-83.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{Try the example to see windows shape.}

\seealso{\code{\link{covspectro}}, \code{\link{dfreq}}, \code{\link{meanspec}},
  \code{\link{spec}}, \code{\link{spectro}}, \code{\link{spectro3D}}}

\examples{
a<-ftwindow(512)
b<-ftwindow(512,wn="bartlett")
c<-ftwindow(512,wn="blackman")
d<-ftwindow(512,wn="flattop")
e<-ftwindow(512,wn="hanning")
f<-ftwindow(512,wn="rectangle")
all<-cbind(a,b,c,d,e,f)
matplot(all,type="l",col=1:6,lty=1:6)
legend(legend=c("hamming","bartlett","blackman","flattop","hanning","rectangle"),
x=380,y=0.95,col=1:6,lty=1:6,cex=0.75)
}

\keyword{ts}
\keyword{IO}