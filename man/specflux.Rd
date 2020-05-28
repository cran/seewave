\name{specflux}

\alias{specflux}

\title{Spectral flux}

\description{Compute spectral flux}

\usage{specflux(wave, f, channel = 1,
         wl = 512, ovlp = 0, wn = "rectangle", flim = NULL,
         norm = FALSE,  p = 2,
         plot = TRUE, xlab = "Times (s)", ylab = "Flux", type = "l", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{channel}{channel of the R object, by default left channel (1).}
  \item{wl}{window length for the analysis
    (even number of points) (by default = 512).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"rectangle"}).}
  \item{flim}{a numeric vector of length 2 to select a frequency band
    (in kHz).}
  \item{norm}{if is \code{TRUE} then the normalised spectra are
    used. The spectra are normalised by their sum.} 
  \item{p}{the norm type, by default = 2.}
  \item{plot}{logical, if \code{TRUE} the spectral flux is displayed
  against time (s) (by default \code{TRUE}).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} parameters.}
}

\details{The spectral flux (\eqn{F}{F}) is the sum of the time (\eqn{t}{t}) derivative of
  the columns -- that is the successive spectra --
  (\eqn{s}{s}) of the normalized short-term Fourier transform (\eqn{z}{z}).\cr
  F is then computed according to:\cr
  \deqn{F = (\sum |s(t+1) - s(t)|^{p})^{\frac{1}{p}}}{%
  F = rowSums(abs(diff(t(z)))^p)^(1/p)}  
}

\value{When \code{plot} is \code{FALSE}, \code{specflux} returns a two-column matrix, the
  first column being time in seconds (\emph{x}-axis) and the second column being the spectral flux (\emph{y}-axis) computed along time.}

\references{
Scheirer E, Slaney M (1997). Construction and evaluation of a robust multifeature speech/music discriminator.
\emph{IEEE International Conference on Acoustics, Speech, and Signal Processing}, 2, 1221-1224.
}

\author{Jérôme Sueur}

\note{The sum of the successive spectral flux values could be used as an
ecoacoustic index, quite close to the acoustic complexity index
(\code{\link{ACI}}). See examples.}

\seealso{\code{\link{spectro}}, \code{\link{ACI}}}

\examples{
## default use
data(tico)
specflux(tico)
## norm 1
specflux(tico, p = 1)
## frequency limit between 2 and 4 kHz
specflux(tico, flim = c(2,4))
## index computation
sum(specflux(tico, plot=FALSE)[,2])
}

\keyword{ts}
