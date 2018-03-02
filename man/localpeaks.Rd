\name{localpeaks}

\alias{localpeaks}

\title{Local maximum frequency peak detection}

\description{This functions searches for local peaks of a frequency spectrum}

\usage{localpeaks(spec, f = NULL, bands = 10, mel = FALSE, plot = TRUE,
xlab = NULL, ylab = "Amplitude", labels = TRUE, ...)}

\arguments{
 \item{spec}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}}. Can be in dB.}
 \item{f}{sampling frequency of \code{spec} (in Hz). Not requested if
  the first column of \code{spec} contains the frequency axis.}
 \item{bands}{a numeric vector. If vector of length 1, then sets the
 number of bands dividing in equal parts the spectrum. If of length > 1,
 then takes the values as kHz limits of the bands dividing the spectrum. These
 bands can be of different size. See details and examples.}
 \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}
 \item{plot}{logical, if \code{TRUE}, a plot showing the peaks is
    returned.}
 \item{xlab}{label of the x-axis.}
 \item{ylab}{label of the y-axis.}
 \item{labels}{logical, if \code{TRUE} peak labels are plotted.}
 \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  The function proceed as follows
  \itemize{
    \item divides the spectrum in bands. The limits of the bands are set
    with the argument \code{bands}. There are two options:
    \itemize{
      \item you set a number of bands with equal size by giving a single
  value to \code{bands}. For instance, setting \code{bands} to a value
  of 10 will slice the spectrum in 10 equal parts and return 10 local
  peaks.
  \item you set the limits of the bands. This is achieve by giving a
  numeric vector to \code{bands}. The limits can follow a regular or
  irregular series. For instance attributing the vector c(0,2,4,8) will
  generate the following bands [0,2[, [2,4[, [4,8] kHz. Be aware that
  the last value should not exceed half the sampling frequency used to
  obtain the spectrum \code{spec}.
}
\item uses the function \code{\link{fpeaks}} with the argument
  \code{nmax} set to 1.
    }
}

\value{A two-column matrix, the first column corresponding to the
  frequency values (x-axis) and the second column corresponding to the
  amplitude values (y-axis) of the peaks.}

\author{Jerome Sueur}

\seealso{\code{\link{fpeaks}}, \code{\link{meanspec}}, \code{\link{spec}}}

\examples{
data(sheep)
spec <- meanspec(sheep, f=8000)
# a specific number of bands with all the same size
localpeaks(spec, bands=5)
# bands directly specified  with a regular sequence
localpeaks(spec, bands=seq(0,8/2,by=0.5))
# bands directly specified  with an irregular sequence
localpeaks(spec, bands=c(0,0.5,1,1.5,3,4))
# Amaj octave bands, note that there is no peak detection
# in the higher part of the spectrum as sequence stops at 3520 Hz
localpeaks(spec, bands=octaves(440, below=3, above=3)/1000)
# melscale
require(tuneR)
mel <- melfcc(sheep, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
melspec.mean <- apply(mel$aspectrum, MARGIN=2, FUN=mean)
melspec.mean <- melspec.mean/max(melspec.mean) # [0,1] scaling 
localpeaks(melspec.mean, f=8000, bands=8)
}

\keyword{dplot}
\keyword{ts}