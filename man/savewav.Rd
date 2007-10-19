\name{savewav}
\alias{savewav}


\title{Save .wav file}

\description{
Save sound data as .wav file
}

\usage{
savewav(wave, f,
filename = NULL)
}


\arguments{
  \item{wave}{data describing the time wave to be exported.}
  \item{f}{sampling frequency of \code{wave} (in Hz).}
  \item{filename}{name of the new file. (by default the name of \code{wave}).}
}

\details{
This functions uses two functions from the package \pkg{sound}:
\code{\link[sound]{Sample}} and \code{\link[sound]{saveSample}}
}

\author{Jérôme Sueur \email{sueur@mnhn.fr}}

\note{The  file automatically owerwrites an existing file
with the same name.
}

\seealso{\code{\link[sound]{as.Sample}},
\code{\link[sound]{saveSample}}, \code{\link{export}}.}

\examples{
a<-synth(f=8000,d=2,cf=2000,plot=FALSE)
# the name of the file is automatically the name of the object
# here: "a.wav"
savewav(a,f=22050)
# if you wish to to change the name, use 'file' argument
savewav(a,f=22050,file="b.wav")
}

\keyword{IO}