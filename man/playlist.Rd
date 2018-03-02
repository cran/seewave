\name{playlist}

\alias{playlist}

\title{Play a list of sound files}

\description{This function works as a playlist, ie it plays back a
  list of sound files.}

\usage{playlist(directory, sample = FALSE, loop = 1)}

\arguments{
  \item{directory}{a character vector indicating the path to the
  directory where sound files to played are saved.}
  \item{sample}{a logical, if \code{TRUE} the order of sounds files to
  be played back is shuffled.}
  \item{loop}{a numeric vector of length 1, number of loops.}
}

\details{The success of using this function depends on the wave player in
use. This works particularly well with SoX under Linux. The type of
files (.mp3, .wav, .ogg etc) depends on the wave player as well)}

\value{None. Listen and enjoy!}

\author{Jérôme Sueur}

\note{The function is mainly based on \code{\link[tuneR]{play}}}

\seealso{\code{\link[tuneR]{play}}, \code{\link{listen}}}

\examples{
\dontrun{
playlist("MyMusic", sample = TRUE, loop=2)
}
}

\keyword{ts}
