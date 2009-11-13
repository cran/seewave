\name{tico}

\alias{tico}

\docType{data}

\title{Song of the bird Zonotrichia capensis}

\description{
  Recording of a song emitted by a male of the neotropical sparrow
  \emph{Zonotrichia capensis}.
}

\usage{data(tico)}

\format{
  A data frame with 39578 observations on the following variable.
  \describe{
    \item{\code{V1}}{a numeric vector}
  }
}

\details{Duration = 1.795 s. Sampling frequency = 22050 hz.}

\source{Recording by Thierry Aubin.}

\examples{
data(tico)
oscillo(tico,f=22050)
}

\keyword{datasets}
