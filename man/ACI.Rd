\name{ACI}

\alias{ACI}

\title{Acoustic Complexity Index}

\description{This function computes the Acoustic Complexity Index (ACI)
  as described in Pieretti \emph{et al}. (2011)}

\usage{ACI(wave, f, wl = 512, ovlp = 0,  wn = "hamming", flim = NULL, nbwindows = 1)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{window length for the analysis
    (even number of points) (by default = 512)}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{flim}{a numeric vector of legnth 2 to select a frequency band (in kHz).}
  \item{nbwindows}{a numeric vector of length 1 specifying the number of
  windows (by default 1, ie a single window including the complete \code{wave} object.}
}

\details{
  The function computes first a short-term Fourier transform and then
  the ACI index. \cr
  The function returns only the ACI total, ACI tot in Pieretti \emph{et
    al}. (2010). \cr
  See the references for details on computation. }

\note{Values returned were checked with the results provided by the
  add-on Soundscapemeter for the software Wavesurfer \url{http://www.speech.kth.se/wavesurfer/}.
}

  
\value{A vector of length 1 returning the ACI total.}

\references{
Pieretti N, Farina A, Morri FD (2011) A new methodology to infer the
singing activity of an avian community: the Acoustic Complexity Index
(ACI). \emph{Ecological Indicators}, 11, 868-873.\cr
Farina A, Pieretti N, Piccioli L (2011) The soundscape methodology for long-term bird monitoring: a Mediterranean Europe case-study. \emph{Ecological Informatics}, 6, 354-363.
}

\author{Laurent Lellouch}

\seealso{\code{\link{spectro}}}

\examples{
data(tico)
ACI(tico)
## dividing the sound sample into 4 windows of equal duration
ACI(tico, nbwindows=4)
## selection of a frequency band
ACI(tico, flim=c(2,6))
}

\keyword{ts}
