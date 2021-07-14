\name{TFSD}

\alias{TFSD}

\title{normalized Time and Frequency Second Derivative}

\description{This function computes the normalized Time and Frequency Second Derivative as described by Aumond et al. (2017).}

\usage{TFSD(wave, f, channel = 1, ovlp = 0,  wn = "hamming", flim = c(2,6), nbwindows = 1)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{channel}{channel of the R object, by default left channel (1).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{flim}{a numeric vector of length 2 to select a frequency band
  (in kHz). Cannot be NULL.}
  \item{nbwindows}{a numeric vector of length 1 specifying the number of
  windows (by default 1, ie a single window including the complete \code{wave} object.}
}

\details{
The TFSD aims at estimating the time of presence of avian or human vocalizations within a sound environment. It calculates the variation in time and frequency of a signal around frequencies of interest, normalized by the spectral time variation of a signal as a whole.\cr

Warning, this index was initially developed to work from a third octave spectrogram with a time sampling of 125 ms.\cr

TFSD is computed according to formulation in reference.\cr

The higher the TFSD varies between 0 and 1, the greater the temporal presence of avian or human vocalizations.  With the default configuration, a TFSD > 0.3 indicates a very important presence time of the vocalizations in the signal. The TFSD is always greater than 0.
}

\value{A numeric vector of length nbwindows giving the TFSD values.}

\references{
Aumond, P., Can, A., De Coensel, B., Botteldooren, D., Ribeiro, C., & Lavandier, C. (2017). Modeling soundscape pleasantness using perceptual assessments and acoustic measurements along paths in urban context. \emph{Acta Acustica united with Acustica}, 12, 50-67.\cr
\cr
Gontier, F., Lavandier, C., Aumond, P., Lagrange, M., & Petiot, J. F. (2019). Estimation of the perceived time of presence of sources in urban acoustic environments using deep learning techniques. \emph{Acta Acustica united with Acustica}, 105(6), 1053-1066.\cr
}

\author{Pierre Aumond, Guillaume Corbeau}

\seealso{\code{\link{ACI}}, \code{\link{NDSI}}}

\examples{
## Note that 'tico' is not a soundscape recording...
data(tico)
TFSD(tico)
## dividing the sound sample into 4 windows of equal duration
TFSD(tico, nbwindows=4)
## selection of a frequency band
TFSD(tico, flim=c(2,6))
}

\keyword{ts}