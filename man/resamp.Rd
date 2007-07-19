\name{resamp}

\alias{resamp}

\title{Resample a time wave}

\description{This function resamples (down- or over-samples) a time wave. This
corresponds to a sampling frequency change.}

\usage{resamp(wave, f, g)}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[Sound]{Sample}} object generated loading a wav file
  with \code{\link[Sound]{loadSample}} (package \pkg{Sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[Sound]{Sample}} object.}
  \item{g}{new sampling frequency of \code{wave} (in Hz).}}

\value{Data are returned as a one-column matrix.}

\author{Jérôme Sueur \email{sueur.jerome@neuf.fr}}

\note{Resampling might change frequency properties of the time wave.}

\examples{
data(peewit)
# downsampling
a<-resamp(peewit,f=22050,g=11025)
# oversampling
b<-resamp(peewit,f=22050,g=44100)
}

\keyword{dplot}
\keyword{ts}
