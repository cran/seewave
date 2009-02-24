\name{rmnoise}

\alias{rmnoise}

\title{Remove noise}

\description{
  This function removes background noise by smoothing 
}

\usage{
rmnoise(wave, f, Sample = FALSE, ...)
}

\arguments{
  \item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{Sample}{a logical, if \code{TRUE} returns an object
    of class \code{\link[sound]{Sample}}}.      
  \item{\dots}{other \code{\link[stats]{smooth.spline}} arguments.}
}

\details{
 This function is based on \code{\link[stats]{smooth.spline}}. You can
 use the arguments of the later to modify the smoothing.
}


\value{A new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{Low frequency noise might not be removed out properly.}

\seealso{\code{\link{afilter}},  \code{\link{noise}}}

\examples{
# synthesis of a 440 Hz sound with background noise
n <- noise(d=1,f=8000)
s <- synth(d=1,f=8000,cf=440)
ns <- n+s
# remove noise (but low frequency content still there)
a <- rmnoise(ns,f=8000)
}

\keyword{ts}

