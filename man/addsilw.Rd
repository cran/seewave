\name{addsilw}

\alias{addsilw}

\title{Add or insert a silence section}

\description{Add or insert a silence section to a time wave.}

\usage{addsilw(wave, f, at = "middle", d = FALSE, plot = FALSE, Sample = FALSE,...)}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{at}{position where to add the silence section (in s).
  Can be also specified as "start", "middle" or "end".}
  \item{d}{duration of the silence section to add (in s).}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot
  of \code{wave} with the new silence section (by default \code{TRUE}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jérôme Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{cutw}},
\code{\link{pastew}}, \code{\link{deletew}}, \code{\link{revw}},
\code{\link{mute}}}

\examples{
data(tico)
addsilw(tico,f=22050,d=0.2)
addsilw(tico,f=22050,at="end",d=0.5)
addsilw(tico,f=22050,at=0.33,d=0.46)
}

\keyword{dplot}

\keyword{ts}

