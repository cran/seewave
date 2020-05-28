\name{beep}

\alias{beep}

\title{Beep sound}

\description{Generate a simple beep to be used as an alert, for instance
  at the end of a loop of when ending up a long script.}

\usage{beep(d = 0.5, f = 8000, cf = 1000)}

\arguments{
  \item{d}{duration (in s)}
  \item{f}{sampling frequency (in Hz)}
  \item{cf}{carrier frequency (in Hz)}
}

\value{Nothing returned, a pure tone sound is played back. The default
  duration is 0.5 s and the default frequency is 1000 Hz}

\author{Jerome Sueur}

\note{The function uses \code{\link{listen}} of \code{seewave} which
  calls \code{play} of \code{tuneR}. You might need to set up your sound
  player with \code{setWavPlayer} of \code{tuneR}.}

\examples{\dontrun{
# default settings
beep()
# change the duration and the frequency
beep(d=1, cf=880)
}
}

\keyword{ts}
