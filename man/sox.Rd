\name{sox}

\alias{sox}

\title{Calls SoX}

\description{This function calls SoX, the Swiss Army knife of sound processing programs.}

\usage{sox(command, exename = NULL, path2exe = NULL)}

\arguments{
  \item{command}{the SoX command to invoke.}
  \item{exename}{a character string specifying the name of the SoX binary file. If \code{NULL}, the default name "sox" will be used for Linux OS and "sox.exe" for Windows OS.}
  \item{path2exe}{a character string giving the path to the SoX binary file. If \code{NULL}, the default path "c:\\sox-14-4-0\\" will be used for Windows OS.}
}

\details{See the documentation of SoX for proper use: \url{http://sox.sourceforge.net/Docs/Documentation}}

\references{SoX website: \url{http://sox.sourceforge.net/}}

\author{Jerome Sueur and Stefanie LaZerte}

\note{Sox must be installed to use this function but not to install the
  package \code{seewave}. As mentioned on the SoX webpage, the primary development
  platform is Linux. Using SoX with Windows from \R might not be
  straightforward. In particular, it is advisable to pay attention to file path and exe name.}

\examples{
\dontrun{
## Generate a simple sound file at 440 Hz
s <- synth(cf=440, f= 8000, d=1, output="Wave")
savewav(s, file="mysound.wav")
## Plays the file
sox("mysound.wav", exe="play")
## Slows down  the  audio  tempo  (but  not  its  pitch)
sox("mysound.wav myslowsound.wav tempo 0.5")
## Cuts the file
sox("myslowsound.wav myslowcuttedsound.wav trim 0.25, 0.75")
## Deletes example files
file.remove("mysound.wav", "myslowsound.wav", "myslowcuttedsound.wav")
}
}

\keyword{IO}
