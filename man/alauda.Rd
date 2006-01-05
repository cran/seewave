\name{alauda}

\alias{alauda}

\docType{data}

\title{Song of the bird Alauda arvensis}

\description{
Recording of a song section emitted by a skylark 
\emph{Alauda arvensis}.
}

\usage{data(alauda)}
\format{
  A data frame with 50768 observations on the following variable.
  \describe{
    \item{\code{V1}}{a numeric vector}
  }
}

\details{
Duration = 2.302 s. Sampling frequency = 22050 hz.
}

\source{
Recording by Thierry Aubin.
}

\examples{
data(alauda)
oscillo(alauda,f=22050)
}

\keyword{datasets}
