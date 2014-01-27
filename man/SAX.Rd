\name{SAX}

\alias{SAX}

\title{Symbolic Aggregate approXimation}

\description{This function converts a numeric times seris into a series
  of letters with a specific length and alphabet.}

\usage{
SAX(x, alphabet_size, PAA_number,
breakpoints = "gaussian", collapse = NULL)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{alphabet_size}{a numeric vector of length 1 setting the size of
  the alphabet.}
  \item{PAA_number}{a numeric vector of length 1 setting the number of
  elements (subsequences) of the Piecewise Aggregate Approximation (PAA).}
  \item{breakpoints}{either a character vector ("gaussian", "quantiles")
  or a numeric vector specifying the sorted values of the breakpoints
  along the distribution of \code{x}. See details and examples.}
  \item{collapse}{a character vector of length 1, specifying the way to
  collapse the output letters, see \code{\link{paste}}. By default letters are
  returned separated.}
}

\details{The SAX method has been developed to reduce the dimensionality of a numerical series into a short chain of characters. SAX follows a two-step process: (1) Piecewise Aggregate Approximation (PAA) and (2) conversion a PAA sequence into a series of letters. \cr

PAA consists in a Z-normalisation, a segmentation of the series of
length \emph{n} into \emph{w} segments, and the computation of each segment average. \cr

The conversion of the PAA into a series of letters is achieved by attributing with
equiprobability each value of the PAA to a letter in reference to a
Gaussian distribution. This process therefore assumes that the
distribution of the numeric series \code{x} follows a Gaussian
distribution. To relax the constraints of normality we here added the possibility to directly work
on the quantiles of the original data distribution or to specify particular breakpoints along the
distribution of \code{x}. See the examples.
}

\value{A character vector of length (when \code{collapse} is
  \code{NULL}) or number of character (when \code{collapse} is
  not \code{NULL}) corresponding to \code{PAA_number} argument. }

\references{
Kasten, E.P., Gage, S.H., Fox, J. & Joo, W. (2012). The remote environmental assessment laboratory's acoustic library: an archive for studying soundscape ecology. \emph{Ecological Informatics}, 12, 50 - 67.\cr

Lin, J., Keogh, E., Lonardi, S., Chiu, B., June (2003). A symbolic representation of time series
with implications for streaming algorithms. Proceedings of the 8th ACM SIGMOD Workshop on Research Issues in Data Mining and Knowledge Discovery. San Diego, California, USA.
}

\author{Laurent Lellouch}

\note{SAX has been used recently to search similar times series in a
  soundcape data base (Kasten et al., 2012).}

\seealso{\code{\link{discrets}}, \code{\link{symba}}, \code{\link{soundscapespec}}}

\examples{
data(tico)
spec <- soundscapespec(tico, plot=FALSE)[,2]
SAX(spec, alphabet = 5, PAA = 10)

# change breakpoints
SAX(spec,  alphabet = 5, PAA = 10, breakpoints="quantiles")
SAX(spec,  alphabet = 5, PAA = 10, breakpoints=c(0, 0.5, 0.75, 1))
SAX(spec,  alphabet = 5, PAA = 10, breakpoints=c(0, 0.33, 0.66, 1))

# different output formats
SAX(spec,  alphabet = 5, PAA = 10, collapse="")
SAX(spec,  alphabet = 5, PAA = 10, collapse="-")
}

\keyword{ts}
