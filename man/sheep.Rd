\name{sheep}

\alias{sheep}

\docType{data}

\title{Sheep bleat}

\description{Recording of a sheep bleat.}

\usage{data(sheep)}

\format{
  A data frame with 19764 observations on the following variable.
  \describe{
    \item{\code{V1}}{a numeric vector}
  }
}

\details{
  Duration = 2.47 s. Sampling frequency = 8000 hz.
}

\source{Recording by Frederic Sebe.}

\examples{
data(sheep)
oscillo(sheep,f=8000)
}

\keyword{datasets}
