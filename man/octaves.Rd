\name{octaves}

\alias{octaves}

\title{Octave values}

\description{This functions returns the frequency values of the octaves
  below and above a specific frequency}

\usage{octaves(x, below = 3, above = 3)}

\arguments{
  \item{x}{a numeric vector, frequency of the note in Hz or Khz.}
  \item{below}{the number of octaves below \code{x}.}
  \item{above}{the number of octaves above \code{x}.}
}

\value{A numeric vector with the octave series in frequency (Hz or kHz
  depending on \code{x} unit).}

\author{Jerome Sueur}

\seealso{\code{\link{notefreq}}}

\examples{
names <- c("C","D","E","F","G","A","B")
values <- c(261.63, 293.66, 329.64, 349.23, 392, 440, 493.88)
res <- sapply(values, FUN=octaves)/1000
op <- par(las=1,mfrow=c(2,1))
par(mar=c(0,4,1,1))
matplot(x=1:7, y=res, t="o", pch=names, xlab="",
    ylab="Frequency (kHz) [linear scale]", col=rainbow(7), xaxt="n")
par(mar=c(4.5,4,0,1))
matplot(x=1:7, y=res, t="o", pch=names, xlab="Octave",
    ylab="Frequency (kHz) [log scale]", col=rainbow(7), ylog=TRUE, log="y")
par(op)
}

\keyword{maths}
