\name{NDSI}

\alias{NDSI}

\title{Normalized Difference Soundscape Index}

\description{This function computes the Normalized Difference Soundscape
  Index as described by Kasten et al. (2012).}

\usage{NDSI(x, anthropophony = 1, biophony = 2:8, max = FALSE)}

\arguments{
  \item{x}{a two-column numeric matrix computed with \code{\link{soundscapespec}}.}
  \item{anthropophony}{a numeric vector defining the frequency band(s) of
  the anthropophony (in kHz).}
  \item{biophony}{a numeric vector defining the frequency band(s) of the
    biophony (in kHz).}
  \item{max}{a logical, if \code{TRUE} then defines the biophony as the
  maximum - not the sum - of the 2 and 8 kHz frequency bands}
}

\details{
NDSI aims at estimating the level of anthropogenic disturbance on the soundscape by
computing the ratio of human-generated (anthropophony) to biological
(biophony) acoustic components found in field collected sound
samples. In terms of frequency, the anthropophony is defined as the [1-2[
kHz frequency bin and the biophony as the [2-8[ kHz frequency bins of a
soundscape frequency spectrum (see \code{\link{soundscapespec}}). \cr

NDSI is computed according to:
\deqn{NDSI = \frac{(biophony - anthropophony)}{(biophony + anthropophony)}}{%
  NDSI = (biophony - anthropophony) / (biophony + anthropophony)
   }
\cr
NDSI varies between -1 and +1, where +1 indicates a signal containing
no anthropophony.
}

\value{A numeric vector of length 1 giving the NDSI value.}

\references{
Kasten, E.P., Gage, S.H., Fox, J. & Joo, W. (2012). The remote
environmental assessment laboratory's acoustic library: an archive for
studying soundscape ecology. \emph{Ecological Informatics}, 12, 50-67.\cr
}

\author{Jerome Sueur}

\seealso{\code{\link{soundscapespec}}, \code{\link{SAX}}, \code{\link{NDSI}}}

\examples{
## Note that 'tico' is not a soundscape recording...
data(tico)
spec <- soundscapespec(tico, plot=FALSE)
NDSI(spec)
NDSI(spec, max=TRUE)
}

\keyword{ts}
