\name{autoc}

\alias{autoc}

\title{Short-term autocorrelation of a time wave}

\description{
This function returns the fundamental frequency of a harmonic time wave. This is
achieved by computing a correlation of the signal with itself
after a time delay. 
}

\usage{
autoc(wave, f, wl = 512, fmin, fmax, threshold = NULL, plot = TRUE,
xlab = "Time (s)", ylab = "Frequency (kHz)", ylim = c(0, f/2000), pb =
FALSE, ...)
}

\arguments{
\item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{length of the window for the analysis (even number of points,
  by default = 512).}
  \item{fmin}{the minimum frequency to detect (in Hz). See details.}
  \item{fmax}{the maximum frequency to detect (in Hz). See details}
  \item{threshold}{amplitude threshold for signal detection (in \%).}
  \item{plot}{logical, if \code{TRUE} plots the fundamental frequency 
  against time (by default \code{TRUE}).}
  \item{xlab}{title of the x-axis.}
  \item{ylab}{title of the y-axis.}
  \item{ylim}{the range of y values.}
  \item{pb}{if \code{TRUE} returns a text progress bar in the console.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{'fmin' and 'fmax' can help by reducing computing time but can
  also produce less accurate results.}

\value{
When \code{plot} is \code{FALSE}, \code{autoc} returns a two-column matrix, the first column corresponding to time in seconds (\emph{x}-axis) and the second column corresponding to 
to fundamental frequency in kHz (\emph{y}-axis).\cr
NA corresponds to pause sections in \code{wave} (see \code{threshold}). 
}

\references{
Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg.}

\author{Jerome Sueur \email{sueur@mnhn.fr} and 
Thierry Aubin \email{thierry.aubin@u-psud.fr}}

\seealso{\code{\link{ceps}}, \code{\link{acf}}}

\examples{
data(sheep)
# fundamental frequency of a sheep
res <- autoc(sheep, f=8000, threshold=5, fmin=100, fmax=700, plot=FALSE)
spectro(sheep, f=8000, ovlp=75, scale=FALSE)
points(res, pch=20)
legend(0.5, 3.6, "Fundamental frequency", pch=20, bty=0, cex=0.7)
}

\keyword{dplot}
\keyword{ts}
