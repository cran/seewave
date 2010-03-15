\name{savewav}

\alias{savewav}

\title{Save a .wav file}

\description{Save sound data as .wav file}

\usage{savewav(wave, f, filename = NULL)}


\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
\item{filename}{name of the new file. (by default the name of \code{wave}).}
}

\details{
This function uses three functions from the package \pkg{tuneR}:
\code{\link[tuneR]{Wave}}, \code{\link[tuneR]{normalize}} and \code{\link[tuneR]{writeWave}}.
}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{The  file automatically owerwrites an existing file
with the same name.}

\seealso{\code{\link{export}}.}

\examples{
require(tuneR)
a<-synth(f=8000,d=2,cf=2000,plot=FALSE)
# the name of the file is automatically the name of the object
# here: "a.wav"
savewav(a,f=22050)
unlink("a.wav")
# if you wish to change the name, use 'file' argument
savewav(a,f=22050,file="b.wav")
unlink("b.wav")
}

\keyword{IO}