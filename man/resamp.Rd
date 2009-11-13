\name{resamp}

\alias{resamp}

\title{Resample a time wave}

\description{This function resamples (down- or over-samples) a time wave. This
  corresponds to a sampling frequency change.}

\usage{resamp(wave, f, g, output="matrix")}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{g}{new sampling frequency of \code{wave} (in Hz).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

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
