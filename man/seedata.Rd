\name{seedata}

\Rdversion{1.1}

\alias{seedata}

\title{A quick look at quantitative data}

\description{See quantitative data at a glance}

\usage{seedata(data, na.rm = FALSE, col = "grey")}

\arguments{
  \item{data}{a numeric vector describing quantitative data.}
  \item{na.rm}{logical, if \code{TRUE} removes \code{NA}.}
  \item{col}{main color.}
}

\details{The red curves depict the corresponding Normal law (same mean and sd as \code{data}).}

\value{A multi-plot graphic is returned.}

\author{Caroline Simonis \email{csimonis@mnhn.fr} and Jerome Sueur \email{sueur@mnhn.fr}.}
  
\examples{seedata(rnorm(1000))}

\keyword{dplot}
