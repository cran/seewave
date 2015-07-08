\name{cutspec}

\alias{cutspec}

\title{Cut a frequency spectrum}

\description{This function can be used to select (cut)
a specific part of a frequency spectrum.}

\usage{cutspec(spec, f = NULL, flim, mel = FALSE, norm = FALSE, PMF = FALSE)}

\arguments{
  \item{spec}{a vector or a two-column matrix set resulting of a spectral analysis (not in dB)
This can be the value obtained with \code{\link{spec}} or \code{\link{meanspec}}.}
  \item{f}{sampling frequency of \code{spec} (in Hz).}
  \item{flim}{a vector of length 2 to specify the new frequency range
    (in kHz).}
  \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}
  \item{norm}{a logical, if \code{TRUE} the spectrum returned is normalised between 0 and 1.}
  \item{PMF}{a logical, if \code{TRUE} the spectrum returned is a probability mass function.}
}

\value{A new spectrum is returned.
The class of the returned object is the one of the input object (\code{spec})}

\author{Jerome Sueur, improved by Laurent Lellouch}

\note{The sampling frequency \code{f} is not necessary if \code{spec} has been obtained with
either \code{spec} or \code{meanspec}.\cr
This function can be used before calling analysis function like \code{\link{sh}} or 
\code{\link{sfm}}. See examples.}

\seealso{\code{\link{spec}}, \code{\link{meanspec}}}

\examples{
data(orni)
a <- meanspec(orni,f=22050,plot=FALSE)
b <- cutspec(a,flim=c(4,8))
## quick check with a plot
plot(b,type="l")
## effects on spectral entropy
sfm(a)
sfm(b)
## mel scale
require(tuneR)
mel <- melfcc(orni, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
melspec.mean <- apply(mel$aspectrum, MARGIN=2, FUN=mean)
c <- cutspec(melspec.mean, f=22050, flim=c(4000,8000), mel=TRUE)
}

\keyword{ts}
