\name{pellucens}

\alias{pellucens}

\docType{data}

\title{Calling song of the tree cricket Oecanthus pellucens}

\description{
  Recording of a calling song section emitted by the European tree cricket \emph{Oecanthus pellucens}.}

\usage{data(pellucens)}

\format{
  A data frame with 36476 observations on the following variable.
  \describe{
    \item{\code{V1}}{a numeric vector}
  }
}

\details{
  Duration = 3.309 s. Sampling frequency = 11025 hz.
}

\source{
  Recording by Jerome Sueur.
}

\examples{
data(pellucens)
oscillo(pellucens,f=11025)
}

\keyword{datasets}
