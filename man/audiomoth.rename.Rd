\name{audiomoth.rename}

\alias{audiomoth.rename}

\title{Rename audiomoth files in a readable format}

\description{This function renames or copies files created with an
  Audiomoth device in a readable format including the data and time of recording.}

\usage{audiomoth.rename(dir, overwrite = FALSE, tz = "", prefix = "")}

\arguments{
  \item{dir}{a character vector, path to directory where the .WAV files are
  stored.}
  \item{overwrite}{a logical, to specify if the files should be renamed
  or copied, if \code{TRUE} the files are copied, if \code{FALSE} the files are renamed.}
  \item{tz}{a character vector defining a time zone specification. See \code{as.POSIXct}}
  \item{prefix}{a charcter vector for a prefix name to be added at the
  beginning of the file name.}
}

\details{The format of the new file names follows the format of the
  SongMeter SM2/SM4 deveices: \code{PREFIX_YYYYMMDD_HHMMSS.wav}.}

\value{1 logical vector indicating which operation succeeded for each of the files attempted.}

\author{Jerome Sueur}

\note{Be careful if you overwrite the files.}

\seealso{\code{\link{audiomoth}}, \code{\link{songmeter}}}

\keyword{input}
\keyword{data}
