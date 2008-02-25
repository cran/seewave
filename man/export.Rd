\name{export}

\alias{export}

\title{Export sound data}

\description{
Export sound data as a text file that can be read by a sound player like
'Goldwave')
}

\usage{
export(wave, f,
filename = NULL,
...)
}

\arguments{
  \item{wave}{data describing the time wave to be exported.}
  \item{f}{sampling frequency of \code{wave} (in Hz).}
  \item{filename}{name of the new file. (by default the name of \code{wave}).}
  \item{\dots}{other \code{\link{write.table}} parameters.}
}

\details{
Creates a new text file with a header describing the main features of the sound
(wave). For instance, for a 2 s sound with a sampling frequency of 8000 Hz,
the header will be: [ASCII 8000Hz, Channels: 1, Samples: 160000, Flags: 0].
This type of file can be read by sound players like
Goldwave (\url{http://www.goldwave.com/}).
}

\author{Jérôme Sueur \email{sueur@mnhn.fr}}

\examples{
a<-synth(f=8000,d=2,cf=2000,plot=FALSE)
export(a,f=8000)
unlink("a.txt")
}

\keyword{IO}
