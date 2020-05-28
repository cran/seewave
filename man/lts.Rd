\name{lts}

\alias{lts}

\title{Long-term spectrogram}

\description{A spectrogram computed over several survey files obtained
  with a Wildlife Acoustics SongMeter recorder}

\usage{
lts(dir, f, wl = 512,
wn = "hanning", ovlp = 0, rmoffset = TRUE, FUN = mean, col = spectro.colors(30),
fftw = FALSE, norm = FALSE, verbose = TRUE,
tlab = "Time", ntann = NULL, flab = "Frequency (kHz)",
recorder = c("songmeter", "audiomoth"), plot = TRUE, ...)
}

\arguments{
  \item{dir}{a character vector, the path to the directory where the
    .wav files are stored or directly the names of the .wav files to be processed.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in the .wav files contained in \code{dir}.}  
  \item{wl}{window length for the analysis
    (even number of points) (by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default
    \code{"hanning"}).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{rmoffset}{a logical to sepcify whether DC offset should be
  removed. By default \code{TRUE}.}
  \item{FUN}{the function to apply to compute the successive frequency spectra, by
    default \code{mean} for a mean spectrum but could be other as
    \code{median} or \code{var} for a median spectrum or variance spectrum.}
  \item{col}{a list of colors or the color palette with a number of colors}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
    library \code{fftw}. See Notes.}
  \item{norm}{a logical, to specify if each mean
    spectrum should be normalised between 0 and 1 (default \code{FALSE})
    before to concatenate the image.}
  \item{verbose}{a logical, if \code{TRUE} (default) the file number and name
    processed are displayed in the console.}
  \item{tlab}{label of the time axis.}
  \item{ntann}{a numeric of length 1, the number of axis annotations (all
    annotations by default).}
  \item{flab}{label of the frequency axis.}
  \item{recorder}{the type of automatic recorder used, either a Wildlife
  SongMeter or a Open Audio deveices Audiomoth.}
  \item{plot}{logical, if \code{TRUE} plots the spectrogram
    (by default \code{TRUE}).}
  \item{\dots}{other \code{\link{image}} graphical parameters.}
}

\details{
  The function reads each .wav file and computes its mean spectrum with
  \code{\link{meanspec}}. The successive mean spectra are then
    concatenated into a single image with the function
    \code{\link{image}}.
  The parameters \code{wl}, \code{ovlp}, and \code{wn} are those of the
  function \code{\link{meanspec}}.
}

\value{
  This function returns a list of three items:
  \item{time}{a numeric vector corresponding to the time axis.}
  \item{freq}{a numeric vector corresponding to the frequency axis.}
  \item{amp}{a numeric or a complex matrix corresponding to the amplitude values.
    Each column is a Fourier transform of length \code{wl/2}.}
}

\author{Jerome Sueur}

\seealso{
  \code{\link{spectro}}, \code{\link{meanspec}},
  \code{\link{image}},
  \code{\link{spectro3D}}, \code{\link{ggspectro}},
  \code{\link{songmeter}}, \code{\link{audiomoth}}
}

\examples{
\dontrun{
## if 'dir' contains a set of files recorded with a Wildlife Acoustics
# songmeter recorder then a direct way to obtain
# the spectrogram of all .wav files is
dir <- "pathway-to-directory-containing-wav-files"
lts(dir)
# to normalise each mean spectrum
lts(dir, norm=TRUE)
# to change the STFT parameters used to obtain each mean spectrum
lts(dir, wl=1024, wn="hamming", ovlp=50)
# to change the colors and the number of time labels and to make it quiet
lts(dir, col=cm.colors(20), ntann=10, verbose=FALSE)
## direct use of files names stored in the working directory
files <- c("S4A09154_20190213_150000.wav", "S4A09154_20190213_153000.wav",
"S4A09154_20190213_160000.wav", "S4A09154_20190213_163000.wav",
"S4A09154_20190213_170000.wav", "S4A09154_20190213_173000.wav",
"S4A09154_20190213_180000.wav", "S4A09154_20190213_183000.wav",
"S4A09154_20190213_190000.wav", "S4A09154_20190213_193000.wav")
lts(files)
}
}

\keyword{dplot}
\keyword{ts}
