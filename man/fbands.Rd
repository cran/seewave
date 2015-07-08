\name{fbands}

\alias{fbands}

\title{Frequency bands plot}

\description{This graphical function returns a frequency spectrum as a
  bar plot.}

\usage{fbands(spec, f = NULL, bands = 10, width = FALSE, mel = FALSE, plot = TRUE,
xlab = NULL, ylab = "Relative amplitude", ...)}

\arguments{
  \item{spec}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}}. Can be in dB.}
  \item{f}{sampling frequency of \code{spec} (in Hz). Not requested if
    the first column of \code{spec} contains the frequency axis.}
  \item{bands}{a numeric vector. If vector of length 1, then sets the
    number of bands dividing in equal parts the spectrum. If of length > 1,
    then takes the values as kHz limits of the bands dividing the spectrum. These
    bands can be of different size. See details and examples.}
  \item{width}{logical, if \code{TRUE} and that \code{bands} is an
    irregular series of values, then the width of the bands will be
    proportional to the frequency  limits defined in \code{bands}.}
  \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}
  \item{plot}{logical, if \code{TRUE}, a plot showing the peaks is
    returned.}
  \item{xlab}{label of the x-axis.}
  \item{ylab}{label of the y-axis.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  The function proceeds as follows
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
    \item uses the function \code{\link{barplot}}.
  }
}


\value{A two-column matrix, the first column corresponding to the
  frequency values (x-axis, mean of the bars limits) and
  the second column corresponding to height values (y-axis) of the bars.}

\author{Jerome Sueur, improved by Laurent Lellouch}

\note{The value below bars is the mean between the corresponding frequency limits.}

\seealso{\code{\link{meanspec}}, \code{\link{spec}}, \code{\link{barplot}}.}

\examples{
data(sheep)
spec <- meanspec(sheep, f=8000, plot=FALSE)
# default plot
fbands(spec)
# setting a specific number of bands
fbands(spec, bands=6)
#setting specific regular bands limits
fbands(spec, bands=seq(0,4,by=0.25))
# some plot tuning
op <- par(las=1)
fbands(spec, bands=seq(0,4,by=0.1),
       horiz=TRUE, col=heat.colors(41),
       xlab="", ylab="",
       cex.axis=0.75, cex.names = 0.75,
       axes=FALSE)
par(op)
# showing or not the width of the bands
oct <- octaves(440,3)/1000
op <- par(mfrow=c(2,1))
fbands(spec, bands=oct, col="blue")
fbands(spec, bands=oct, width = TRUE, col="red")
par(op)
# kind of horizontal zoom
op <- par(mfrow=c(2,1))
fbands(spec, bands=seq(0,4,by=0.2), col=c(rep(1,10),
   rep("orange",5),rep(1,5)), main="all frequency range")
fbands(spec, bands=seq(2,3,by=0.2),
   col="orange", main="a subset or zoom in")
par(op)
# kind of dynamic frequency bands
specs <- dynspec(sheep, f=8000, plot= FALSE)$amp
out <- apply(specs, f=8000, MARGIN=2,
     FUN = fbands, bands = seq(0,4,by=0.2),
     col = 1, ylim=c(0,max(specs)))
# mel scale
require(tuneR)
mel <- melfcc(sheep, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
melspec.mean <- apply(mel$aspectrum, MARGIN=2, FUN=mean)
melspec.mean <- melspec.mean/max(melspec.mean) # [0,1] scaling 
fbands(melspec.mean, f=8000, bands=8)
}

\keyword{dplot}
\keyword{ts}
