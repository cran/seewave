\name{scd}

\alias{scd}

\title{Soundcape chord diagram}

\description{Computes and displays a chord diagram of a set of audio
  files or of a set segments extracted from a single audio file.}

\usage{
scd(input, f,  sl, wl = 512, wn = "hanning", ovlp = 0, flim = NULL,
rmoffset = TRUE, threshold = NULL, HCA = TRUE,
grid.col = terrain.colors, names, plot = TRUE, verbose = TRUE, ...)
}

\arguments{
  \item{input}{a character vector, either (i) the path to the directory where
    .wav files are stored, (ii) directly the names of the .wav files to be
    processed, or (iii) a single .wav file to be segmented by the
    duration set with the argument \code{sl}.}
  
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in the .wav files contained in \code{input}.}  

  \item{sl}{segment length in duration if \code{input} is a single file.}
  
  \item{wl}{window length for the analysis spectral
    (even number of points) (by default = 512).}

  \item{wn}{window name for the spectral analysis, see
  \code{\link{ftwindow}} (by default \code{"hanning"}).}

  \item{ovlp}{overlap between two successive windows (in \%) for the
    spectral analysis.}

  \item{flim}{a numeric vector of length 2 to select a frequency band
    (in kHz).}
  
  \item{rmoffset}{a logical to sepcify whether DC offset should be
    removed. By default \code{TRUE}.}

  \item{threshold}{a numeric value in ]0,1[ to be applied to the
    similarity distance. All similairty distances below this threshold
    will not depicted.}
  
  \item{HCA}{logical, if \code{TRUE} the sectors are colored according
    to clusters automatically obtained with a hierarchical cluster analysis (HCA).}
  
  \item{grid.col}{name of color palette to color the sectors and
    the links). By default \code{terrain.colors}.}

  \item{names}{names of the sectors, if empty then the names of the .wav
    files or the time stamps of the segments.}

  \item{plot}{logical, if \code{TRUE} plots the chord diagram.}

  \item{verbose}{logical, if \code{TRUE} prints the computation progress
    of the mean spectra.}
  
  \item{\dots}{other \code{chordDiagram} parameters.}
}

\details{
  The soundscape chord diagram (SCD) aims at representing similarities
  between audio files or audio segments extracted from a single audio file. The
  mean frequency spectrum of each file/segment is computed using a
  STFT. These frequency spectra are then (1) pairwised compared using a
  similarity distance (see function \code{\link{diffcumspec}}, and (2)
  automatically clustered with a hierarchical cluster analysis (HCA) (see function
  \code{HCPC} of \code{FactoMiner}). The
  resulting similarity matrix is then given as an input to the function
  \code{chordDiagram}. The width of the sectors and the links are
  based on the spectral similarity matrix. The color of the sectors and the links
  follow the HCA classification.
}

\value{
  The function returns a list of two items:
  \item{m}{spectral similarity matrix}
  \item{resHCA}{the classification result of the HCA, if \code{HCA} is \code{TRUE}}
}


\references{de Baudouin, A, Couprie P, Michaud F, Haupert S, Sueur J --
  Similarity visualization of nature and music soundscapes, in prep.}

\author{Adèle de Baudouin, Jérôme Sueur}

\note{The function call the function \code{HCPC} of the package
  \code{FactoMineR} and the function \code{chordDiagram} of the
  package \code{circlize}.
}

\seealso{\code{\link{diffcumspec}}}

\examples{
\dontrun{
## 1 ##
# if 'dir' contains a set of files recorded with a Wildlife Acoustics
# songmeter recorder or an Audiomoth  then a direct way to obtain
# the soundscape chord diagram (SCD)  of all .wav files is
dir <- "pathway-to-directory-containing-wav-files"
scd(dir)
# to change the STFT parameters used to obtain each mean spectrum
lts(dir, wl=1024, wn="hamming", ovlp=50)
# to select only high similarities, here above 0.6
scd(dir, threshold=0.6)
# to change the colors
scd(dir, grid.col=colorRampPalette(c("darkblue", "yellow", "grey")))
# to name manually the sectors
scd(dir, names=as.character(0:23))
# to name automatically the sectors from the name of songmeter files
# here according to hour of recording 
scd(dir, names=as.character(songmeter(files)$hour))  

## 2 ##
# to directly use files names stored in the working directory
files <- c("S4A09154_20190213_150000.wav", "S4A09154_20190213_153000.wav",
"S4A09154_20190213_160000.wav", "S4A09154_20190213_163000.wav",
"S4A09154_20190213_170000.wav", "S4A09154_20190213_173000.wav",
"S4A09154_20190213_180000.wav", "S4A09154_20190213_183000.wav",
"S4A09154_20190213_190000.wav", "S4A09154_20190213_193000.wav")
scd(files)

## 3 ##
# to  use of single files which is segmented in successive time segments
# lasting each 60 s
file <- "a-very-nice-soundscape.wav")
scd(file, sl = 60)
}

}

\keyword{plot}
\keyword{ts}
