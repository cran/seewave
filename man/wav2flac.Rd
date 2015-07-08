\name{wav2flac}

\alias{wav2flac}

\title{wav-flac file conversion}

\description{This function converts .wav files into .flac files and reversely}

\usage{wav2flac(file, reverse = FALSE, overwrite = FALSE,
exename = NULL, path2exe = NULL)}

\arguments{
  \item{file}{the .wav or .flac file to convert.}
  \item{reverse}{logical, if \code{TRUE} converts a .flac file into a .wav file.}
  \item{overwrite}{logical, if \code{TRUE} overwrites the file to convert.}
  \item{exename}{a character string specifying the name of the FLAC binary file. If \code{NULL}, the dedault name "flac" will be used for Linux OS and "flac.exe" for Windows OS.}
  \item{path2exe}{a character string giving the path to the FLAC binary file. If \code{NULL}, the dedault path "c:/Program Files/FLAC/" will be used for Windows OS.}
}

\details{The function runs FLAC. FLAC has then to be installed first: \url{http://flac.sourceforge.net/}, if not the function will not work.}

\value{A new file is created.}

\references{FLAC website: \url{http://flac.sourceforge.net/}}

\note{FLAC must be installed to use this function but not to install the
  package \code{seewave}. Free Lossless Audio Codec (FLAC) is a file format by Josh Coalson for lossless audio data compression. FLAC reduces bandwidth and storage requirements without sacrificing the integrity of the audio source. Audio sources encoded to FLAC are typically reduced in size 40 to 50 percent.}

\author{Luis J. Villanueva-Rivera}

\seealso{\code{\link{savewav}}}

\examples{
\dontrun{
# synthesis of a 1kHz sound
a<-synth(d=10,f=8000,cf=1000)
# save it as a .wav file in the default working directory
savewav(a,f=8000)
# compress it to FLAC format and overwrite on the file a.wav
wav2flac("a.wav", overwrite=TRUE)
# back to .wav format
wav2flac("a.flac", reverse=TRUE)
# remove the files
unlink(c("a.wav","a.flac"))
}
}

\keyword{IO}
