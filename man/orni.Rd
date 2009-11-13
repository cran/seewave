\name{orni}

\alias{orni}

\docType{data}

\title{Song of the cicada Cicada orni}

\description{
  Recording of a calling song section of the Mediterranean cicada 
  \emph{Cicada orni}.
}

\usage{data(orni)}

\format{
  A data frame with 15842 observations on the following variable.
  \describe{
    \item{\code{V1}}{a numeric vector}
  }
}

\details{
  Duration = 0.719 s. Sampling frequency = 22050 Hz.
}

\source{
  Recording by Jerome Sueur.
}

\examples{
data(orni)
oscillo(orni,f=22050)
}

\keyword{datasets}
