\name{savewav}

\alias{savewav}

\title{Save a .wav file}

\description{Save sound data as .wav file}

\usage{savewav(wave, f = NULL, filename = NULL, rescale = NULL)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{filename}{name of the new file. (by default the name of
    \code{wave}).}
  \item{rescale}{a numeric vector of length 2 giving the lower (negative
  value) and upper (positive value) amplitude limits of the .wav file to be exported.}
}

\details{
This function uses three functions from the package \pkg{tuneR}:
\code{\link[tuneR]{Wave}}, \code{\link[tuneR]{normalize}} and \code{\link[tuneR]{writeWave}}.
}

\author{Jerome Sueur \email{sueur@mnhn.fr}, Ethan C. Brown for the
  argument 'rescale'}

\note{The  file automatically owerwrites an existing file
  with the same name.\cr
The amplitude (volume) of the .wav file is normalized by defaults but can be changed with the argument
  \code{rescale}. See examples}

\seealso{\code{\link{export}}.}

\examples{
require(tuneR)
a<-synth(f=8000,d=2,cf=2000,plot=FALSE)
# the name of the file is automatically the name of the object
# here: "a.wav"
savewav(a,f=22050)
unlink("a.wav")
# if you wish to change the name, use the 'file' argument
savewav(a,f=22050,file="b.wav")
unlink("b.wav")
# if you wish to change the amplitude of the file, use the argument 'rescale'
# this will turn down the volume of a 16 bit sound
# which amplitude was originally ranging between -2^15 and +2^15
savewav(a, f=22050, file="c.wav", rescale=c(-1500,1500))
unlink("c.wav")
}

\keyword{IO}