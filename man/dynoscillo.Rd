\name{dynoscillo}

\alias{dynoscillo}

\title{Dynamic oscillogram}

\description{This graphical function displays a time wave as an windowed oscillogram.}

\usage{dynoscillo(wave, f, wd = NULL, wl = NULL, wnb = NULL, title = TRUE, ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wd}{a numerical vector, duration of the window (in seconds)}
  \item{wl}{a numerical vector, length of the window (in number of points).}
  \item{wnb}{a numerical vector, number of windows (no unit).}
  \item{title}{a logical, if \code{TRUE} displays a title with information regarding window size and number.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}
\details{
The arguments \code{wd}, \code{wl} and \code{wn} have to be used isolated, not in conjunction. They basically do the same, ie they set the duration of the zooming window that is slided along the signal. For instance, for a 5 seconds sound with a sampling rate (\code{f}) at 44.1 kHz, \code{wl = 4096} is equivalent to \code{wd = 4096 / 44100 = 0.093 s} and equivalent to \code{wnb = 5*4096 / 44100 = 53}.} 

\author{Jerome Sueur}

\note{This function requires the package \code{rpanel}.}

\seealso{\code{\link{oscillo}}, \code{\link{oscilloST}}, \code{\link{dynspec}}.}

\examples{\dontrun{
require(rpanel)
data(tico)
dynoscillo(tico, wn=4)}
}

\keyword{dplot}
\keyword{ts}
