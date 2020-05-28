\name{timelapse}

\alias{timelapse}

\title{Time lapse}

\description{Append successive input sounds into a single output sound}

\usage{
timelapse(dir, from = 1, to = Inf,
units = c("samples", "seconds", "minutes", "hours"), verbose = TRUE)
}

\arguments{
  \item{dir}{a character vector, the path to the directory where the
    .wav files are stored or directly the names of the .wav files to be appended.}
  \item{from}{where to start reading the input files, in \code{units}. See
  \code{readWave} of the package \code{tuneR}.}
  \item{to}{where to stop reading, in \code{units}. See
  \code{readWave} of the package \code{tuneR}.}
  \item{units}{time units in which \code{from} and \code{to} is given, the default is
          \code{"samples"}, but can be set to time intervals such as \code{"seconds"}. See
  \code{readWave} of the package \code{tuneR}.}
  \item{verbose}{a logical, if \code{TRUE} (default) the file number and name
    processed are displayed in the console.}
}

\details{The function takes the \code{.wav} files which names are provided in the
  argument \code{dir} and append (paste) them successively so that a single object
  is obtained. This can be used to produce sound time lapse based on a series
  of ordered files as those produced by an automatic recorder
  (e.g. SongMeter of the society 'Wildlife Acoustics').\cr
  Only a section of each file can be extracted by using the arguments
  \code{from} and \code{to}. 
  The function is based on \code{readWave} and \code{bind} of the package \code{tuneR}.
}

\value{A \code{Wave} object, a class defined in the package \code{tuneR}.}

\author{Jérôme Sueur}

\note{The characteristics (sampling rate, number of bits, stereo/mono)
  of the output object are those of the .wav file.\cr
  The files should be alphabatically ordered according to time to
  ensure a proper time lapse.\cr
  You should use either \code{savewav} or \code{writeWave} to save the results as a \code{.wav} file. 
}


\seealso{\code{\link{pastew}}}

\examples{
\dontrun{
## if 'dir' contains a set of files recorded with a Wildlife Acoustics
# songmeter recorder then a direct way to obtain
# the spectrogram of all .wav files is
dir <- "pathway-to-directory-containing-wav-files"
res <- timelapse(dir)
# to extract a selection of each file (here a section starting
# at 10 s and ending at 12 s)
res <- timelapse(dir, from=10, to=12, unit="seconds")
}
}

\keyword{ts}