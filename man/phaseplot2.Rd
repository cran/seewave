\name{phaseplot2}

\alias{phaseplot2}

\title{Phase-phase 2D plot of a time wave}

\description{This functions returns a 2D representation of a time wave
  against a delayed version of itself.}

\usage{
phaseplot2(wave, f, tau = 1, type = "l",
xlab = "x(t)",
ylab = paste("x(t+", tau, ")", sep = ""), ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{tau}{the time delay to apply in number of samples.}
  \item{type}{type of plot that should be drawn. See \code{\link{plot}} for details (by default "l" for lines).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{\dots}{other \code{\link{plot}} parameters.}
}

\details{The principle consists in displaying in a single x-y graph the
  original time wave with a delayed version of itself. The delay is
  controlled with the argument \code{tau} that needs to be specified in
  number of samples. The conversion of \code{tau} in second is obtained by calculating
  \code{tau/f}, with \code{f} the sampling frequency.}

\value{Nothing is returned except an x-y plot.}

\references{Kantz H, Schreiber T (2003) \emph{Non linear time series analysis}. Cambridge University Press.}

\author{Jerome Sueur}

\note{Phase-phase plot can be used to test non-linearity.}

\seealso{\code{\link{phaseplot}}}

\examples{
s <- synth(d=0.05, f=44100, cf=440, out="Wave")
n <- noisew(d=0.05, f=44100, out="Wave")
par(mfrow=c(2,1))
phaseplot2(s)
phaseplot2(n)
}

\keyword{dplot}
\keyword{ts}
