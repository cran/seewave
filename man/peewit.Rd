\name{peewit}
\alias{peewit}
\docType{data}
\title{Song of the bird Vanellus vanellus}
\description{
Recording of a song emitted by a peewit (lapwing) male
\emph{Vanellus vanellus}
}
\usage{data(peewit)}
\format{
  A data frame with 15561 observations on the following variable.
  \describe{
    \item{\code{V1}}{a numeric vector}
  }
}
\details{
Duration = 0.706 s. Sampling frequency = 22050 hz.
}
\source{
Recording by Thierry Aubin.
}
\examples{
data(peewit)
oscillo(peewit,f=22050)
}
\keyword{datasets}
