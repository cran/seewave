\name{crest}

\Rdversion{1.1}

\alias{crest}

\title{Crest factor and visualization}

\description{This function returns the crest factor and localizes the
  different crest(s).}

\usage{crest(wave, f, plot = FALSE, col = 2, cex = 3, symbol = "*", ...)}

\arguments{
\item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{plot}{if \code{TRUE} plots the oscillogram of \code{wave} and
  indicates the location of the crest(s)}
  \item{col}{color of the symbol indicating the localisation of the crest(s)}
  \item{cex}{symbol magnification}
  \item{symbol}{symbol indicating the localisation of the crest(s)}
  \item{\dots}{other}
}
\details{
  The crest factor of a time series s is calculated according to:\cr
  \deqn{C = \frac{max(s)}{rms(s)}}{% 
      C = max(abs(s))/rms(s)}
    with rms the root-mean-square (see \code{\link{rms}}).
}

\value{
The function returns a list of three items
\item{C}{crest factor}
\item{val}{value of the crest(s)}
\item{loc}{location of the crest(s)}
}

\references{Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.
}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{There might be several crests (maxima) along the time wave but
  there is a single crest factor.}

\seealso{\code{\link{oscillo}}, \code{\link{rms}}}

\examples{
data(tico)
crest(tico, f=22050)
# see the crest location and change the default graphical parameters
crest(tico, f=22050, plot=TRUE, sym="-")
}

\keyword{ts}
