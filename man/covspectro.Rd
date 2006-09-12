\name{covspectro}

\alias{covspectro}

\title{Covariance between two spectrograms}

\description{
This function tests the similarity between two spectrograms by returning
their maximal covariance and the time shift related to it. 
}

\usage{
covspectro(wave1, wave2, f, wl = 512, wn = "hanning", n, plot = TRUE, plotval = TRUE,
method = "spearman", col = "black", colval = "red", cexval = 1,
fontval = 1, xlab = "Time (s)",
ylab = "Normalised covariance (cov)", ...)
}

\arguments{
  \item{wave1}{data or a \code{\link[Sound]{Sample}} object created loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}) describing a first time wave.}
  \item{wave2}{data or a \code{\link[Sound]{Sample}} object created loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}) describing a second time wave.}
  \item{f}{sampling frequency of \code{wave1} and \code{wave2} (in Hz).}
  \item{wl}{length of the window for the analysis
      (even number of points).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{n}{number of covariances computed between \code{wave1} and \code{wave2}
      when sliding \code{wave2} along \code{wave1}.}
  \item{plot}{logical, if \code{TRUE} plots r values against frequency shift
      (by default \code{TRUE}).}
  \item{plotval}{logical, if \code{TRUE} adds to the plot maximum R value
  and frequency offset (by default \code{TRUE}).}  
  \item{method}{a character string indicating which correlation coefficient is
      to be computed ("pearson", "spearman", or "kendall")
      (see \code{\link{cor}}).}
  \item{col}{colour of r values.}
  \item{colval}{colour of r max and frequency offset values.}
  \item{cexval}{character size of r max and frequency offset values.}
  \item{fontval}{font of r max and frequency offset values.}
  \item{xlab}{title of the frequency axis.}
  \item{ylab}{title of the r axis.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
Successive covariances between the spectrogram of \code{wave1} and
the spectrogram of \code{wave2} are computed when regularly sliding
forward and backward \code{wave2} along \code{wave1}.\cr
The maximal covariance is obtained at a particular shift (time offset). 
This shift may be positive or negative.\cr 
\code{n} sets in how many steps \code{wave2}will be slided along \code{wave1}.
Time process can be then decreased by setting low \code{n} value.\cr
Inverting \code{wave1} and \code{wave2} may give slight different results.
}

\value{
If \code{plot} is \code{FALSE}, \code{covspectro} returns a list containing
three components:
  \item{cov}{the successive covariance values between \code{wave1} and \code{wave2}.}
  \item{covmax}{the maximum covariance between \code{wave1} and \code{wave2}.}
  \item{t}{the time offset corresponding to \code{cov}.}
}

\references{
Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg.}

\author{Jérôme Sueur \email{jerome.sueur@univ-tours.fr}}

\seealso{\code{\link{corspec}}, \code{\link{corenv}}, \code{\link{spectro}},
        \code{\link{cor}},}

\examples{
# covariance between two echemes of a cicada calling song
data(orni)
orni1<-cutw(orni,f=22050,from=0,to=0.1,plot=FALSE)
orni2<-cutw(orni,f=22050,from=0.15,to=0.25,plot=FALSE)
covspectro(orni1,orni2,f=22050,wl=512,n=37,type="l")
}

\keyword{dplot}
\keyword{ts}
