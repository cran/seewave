\name{discrets}

\alias{discrets}

\title{Time series discretisation}

\description{This function transforms a numeric (time) series into a sequence of symbols}

\usage{discrets(x, symb = 5, collapse = TRUE, plateau=1)}

\arguments{
  \item{x}{a \code{vector}, a \code{matrix} (first column),
    an object of class \code{ts}, \code{Sample} (left channel),
    or \code{\link[tuneR]{Wave}} (left channel).}
  \item{symb}{the number of symbols used for the discretisation,
    can be set to 3 or 5 only.}
  \item{collapse}{logical, if \code{TRUE}, the symbols are pasted in a
    character string of length 1.}
  \item{plateau}{a numeric vector of length 1 taking the values \code{1} or \code{2}
    only. See details.}
  }

\details{The function partitions the numeric (time) series into a sequence of finite number of symbols.
  These symbols result of the comparaison of each series value with its temporal neighbours.\cr
  They are two discretisations available:\cr
  when \code{symb} is set to 3, each value will be replaced by either:\cr
  - \emph{I} if the series is \emph{I}ncreasing,\cr
  - \emph{D} if the series is \emph{D}ecreasing,\cr
  - \emph{F} if the series remains \emph{F}lat,\cr
  when \code{symb} is set to 5, each value will be replaced by either:\cr
  - \emph{I} if the series is \emph{I}ncreasing,\cr
  - \emph{D} if the series is \emph{D}ecreasing,\cr
  - \emph{F} if the series remains \emph{F}lat,\cr
  - \emph{P} if the series shows a \emph{P}eak,\cr
  - \emph{T} if the series shows  a \emph{T}rough.\cr
  \cr
  The argument \code{plateau} can be used to control the way a plateau
  is encoded. A plateau is an elevated flat region that can be either
  considered a 'flat peak' encoded as \emph{PF...FP} (\code{plateau
  = 1}) or as an increase, a flat region and a decrease encoded as
  \emph{IF...FD} (\code{plateau = 1}. The default value (\code{plateau
  = 1}) refers to Cazelles \emph{et al}. (2004).
  
}

\value{A character string of length 1 if \code{collapse} is \code{TRUE}.
  Otherwise, a character string of length \emph{n-2} if \code{symbol=5}
  (the first and last values cannot be replaced with a symbol)
  or \emph{n-1} if \code{symbol=3} (the first value cannot be replaced with a symbol.)
}

\references{Cazelles, B. 2004
  Symbolic dynamics for identifying similarity between rhythms of ecological time series.
  \emph{Ecology Letters}, \bold{7}: 755-763.}

\author{Jerome Sueur, improved by Laurent Lellouch}

\seealso{\code{\link{symba}}}

\examples{
# a random variable
discrets(rnorm(30))
discrets(rnorm(30),symb=3)
# a frequency spectrum
data(tico)
spec1<-spec(tico,f=22050,at=0.2,plot=FALSE)
discrets(spec1[,2])
}
\keyword{ts}
