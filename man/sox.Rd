\name{sox}

\alias{sox}

\title{Calls SoX}

\description{This function calls SoX, the Swiss Army knife of sound processing programs.}

\usage{sox(command, exename = NULL, path2exe = NULL, option = NULL,
shQuote_type = NULL)}

\arguments{
  \item{command}{the SoX command to invoke.}
  \item{exename}{a character string specifying the name of the SoX binary file. If \code{NULL}, the default name "sox" will be used for Linux OS.}
  \item{path2exe}{a character string giving the path to the SoX binary
    file}
g  \item{option}{option to be passed to the SoX command}
  \item{shQuote_type}{type of shell quotes ("cmd" or "cmd2", for Windows
  OS; "sh" or "csh" Unix OS)}
}

\details{See the documentation of SoX for proper use.}

\references{\url{https://en.wikipedia.org/wiki/SoX}}

\author{Jerome Sueur, Stefanie LaZerte, Andre Mikulec}

\note{Sox must be installed to use this function but not to install the
  package \code{seewave}. As mentioned on the SoX webpage, the primary development
  platform is Linux. Using SoX with Windows from \R might not be
  straightforward. In particular, it is advisable to pay attention to file path and exe name.}

\examples{
\dontrun{
##########
## data ##
##########
## Generate a simple sound file at 440 Hz
s <- synth(cf=440, f=8000, d=1, fm=c(0,0,1000,0,0), output="Wave")
savewav(s, file="mysound.wav")
##############
## Linux OS ##
##############
## Play the file
sox("mysound.wav", exename="play")
## Slow down  the  audio  tempo  (but  not  its  pitch)
sox("mysound.wav myslowsound.wav tempo 0.5")
## Cut the file
sox("myslowsound.wav myslowcutsound.wav trim 0.25 0.75")
################
## Windows OS ##
################
## path with simple slash
path <- "C:/Program Files (x86)/sox-14-4-2"
## or path with double backslash
## path <- "C:\\Program Files (x86)\\sox-14-4-2"
sox("mysound.wav", path2exe=path, option="-t waveaudio")
## with the option directly passed to the command 
sox("mysound.wav -t waveaudio", path2exe=path)
## Slow down  the  audio  tempo  (but  not  its  pitch)
sox("mysound.wav myslowsound.wav tempo 0.5", path2exe=path)
## Cut the file
sox("myslowsound.wav myslowcutsound.wav trim 0.25 0.75", path2exe=path)
###########
## clean ##
###########
file.remove("mysound.wav", "myslowsound.wav", "myslowcutsound.wav")
}
}

\keyword{IO}
