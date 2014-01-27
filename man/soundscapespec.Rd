\name{soundscapespec}

\alias{soundscapespec}

\title{Soundscape frequency spectrum of a time wave}

\description{This function returns a kHz binned spectrum as described by Kasten et
  al. (2012) for the description of a soundscape.}

\usage{
soundscapespec(wave, f, wl = 1024, wn = "hamming", ovlp = 50,
plot = TRUE, xlab = "Frequency (kHz)", ylim = c(0, 1), ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{length of the window for the analysis (even number of points, by default = 1024).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hamming"}).}
  \item{ovlp}{overlap between two successive analysis windows (in \%),
  by default = 50\%.}
  \item{plot}{if \code{TRUE} returns a barplot.}
  \item{xlab}{title of the barplot x axis.}
  \item{ylim}{range of the barplot y axis.}
  \item{\dots}{other \code{\link{barplot}} graphical parameters.}
}

\details{The soundscape frequency spectrum is based on the computation of a spectrogram power
  spectral density using Welch'smethod (Welch & June, 1967). Parameters used in
  Kasten et al. (2012) were a Hamming window of 1024 samples with 50\% of overlap and are used here as default values.}

\value{A two-column numeric matrix, the first column returning the
  frequency (kHz) bands and the second column returning the power
  value within each frequency band.\cr
  A barplot is returned when \code{plot} is \code{TRUE}.}

\references{
Kasten, E.P., Gage, S.H., Fox, J. & Joo, W. (2012). The remote
environmental assessment laboratory's acoustic library: an archive for
studying soundscape ecology. \emph{Ecological Informatics}, 12, 50-67.\cr
Welch, P.D., June (1967). The use of the fast Fourier transform for the estimation
of power spectra: a method based on time-averaging over short, modified
periodograms. \emph{IEEE Transactions on Audio and Electroacoustics}, 15: 70-73.
}

\author{Jerome Sueur and Eric Kasten}

\seealso{\code{\link{spec}}, \code{\link{meanspec}}, \code{\link{SAX}}, \code{\link{NDSI}}}


\examples{
## Note that 'tico' is not a soundscape recording...
data(tico)
soundscapespec(tico, plot=TRUE, col="darkgreen")
}

\keyword{dplot}
\keyword{ts}
