\name{symba}

\alias{symba}

\title{Symbol analysis of a numeric (time) series}

\description{This function analyses one or two sequences of symbols from numeric (time) series.}

\usage{
symba(x, y = NULL, symb = 5, collapse = TRUE, entropy = "abs",
plot = FALSE, type = "l", lty1 = 1, lty2 = 2, col1 = 2, col2 = 4,
cex1 = 0.75, cex2= 0.75, xlab = "index", ylab = "Amplitude", legend=TRUE, ...)
}

\arguments{
  \item{x}{a first R object.}
  \item{y}{a second R object}
  \item{symb}{the number of symbols used for the discretisation, can be set to 3 or 5 only.}
  \item{collapse}{logical, if \code{TRUE}, the symbols are pasted in a character string of length 1.}
  \item{entropy}{either "abs" for an absolute value or "rel" for a relative value, i. e. between 0 and 1.}
  \item{plot}{logical, if \code{TRUE} plots the series \code{x} (and \code{y}) and the respective symbols.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{lty1}{line type of the object \code{x} if \code{type}="l".}
  \item{lty2}{line type of the object \code{y} if \code{type}="l".}
  \item{col1}{colour of the object \code{x}.}
  \item{col2}{colour of the object \code{y}.}
  \item{cex1}{character size of \code{x} symbols.}
  \item{cex2}{character size of \code{y} symbols.}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{legend}{logical, if \code{TRUE} and if \code{y} is not \code{NULL} adds a legend to the plot.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{The analysis consists in transforming the series into a sequence of symbols (see the function
  \code{\link{discrets}}) and in computing the absolute frequency of each symbol within the sequence.\cr
  The entropy (\emph{H}) is then calculated using the symbol frequencies.
  Using the argument \code{entropy}, the entropy can be expressed along an absolute scale or as a relative
  value varying between 0 and 1.\cr
  If two numeric (time) series are provided (\code{x} and \code{y}) the absolute symbol
  frequencies and entropy of each series is returned. Besides the mutual information (\emph{I})
  is estimated according to:\cr
  \deqn{I = H_{x} + H_{y} - H{xy}}{% 
    I = Hx + Hy - Hxy}
  with \emph{Hx} the entropy of \code{x} symbol series,
    \emph{Hy} the entropy of \code{y} symbol series, and \`emph{Hxy}$ the joint entropy
    of \code{x} and \code{y} symbol series.
}

\value{
  If \code{y} is \code{NULL} a list of three items is returned (s1, freq1, h1).\cr
  If \code{y} is not \code{NULL}, a list of 6 items is returned (s1, freq1, h1, s2, freq2, h2, I):
  \item{s1}{the sequence of symbols of \code{x},}
  \item{freq1}{the absolute frequency of each \code{x} symbol,}
  \item{h1}{the entropy of \code{x} symbol sequence,}
  \item{s2}{the sequence of symbols of \code{y},}
  \item{freq2}{the absolute frequency of each \code{y} symbol,}
  \item{h2}{the entropy of \code{y} symbol sequence,}
  \item{I}{the mutual information between \code{x} and \code{y}.}
}

\references{Cazelles, B. 2004
  Symbolic dynamics for identifying similarity between rhythms of ecological time series.
  \emph{Ecology Letters}, 7: 755-763.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{
  It might be useful to round the values of the input series (see examples).\cr
  The mutual information (\emph{I}) should increase with the similarity
  between the series to compare (\code{x} and \code{y}).
}

\seealso{\code{\link{discrets}}, \code{\link{SAX}}}

\examples{
# analysis of a frequency spectrum
data(tico)
spec1<-spec(tico,f=22050,at=0.2,plot=FALSE)
symba(spec1[,2],plot=TRUE)
# it might be better to round the values
symba(round(spec1[,2],2),plot=TRUE)
# in that case the symbol entropy is almost similar to the spectral entropy
symba(round(spec1[,2],2),entrop="rel")$h1
sh(spec1)
# to compare two frequency spectra
spec2<-spec(tico,f=22050,wl=512,at=1.1,plot=FALSE)
symba(round(spec1[,2],2),round(spec2[,2],2),plot=TRUE)
}
\keyword{ts}
