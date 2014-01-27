\name{TKEO}

\alias{TKEO}

\title{Teager-Kaiser energy tracking operator}

\description{This function computes the Teager-Kaiser energy operator.}
\usage{
TKEO(wave, f, m = 1, M = 1, plot = TRUE,
xlab = "Time (s)", ylab = "Energy",
type = "l", bty = "l", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{m}{a numeric vector of length 1 for the lag parameter. See details.}
  \item{M}{a numeric vector of length 1 for the exponent parameter. See details.}
  \item{plot}{logical, if \code{TRUE} returns a plot of the TK energy
  along time (by default \code{TRUE}).}
  \item{xlab}{Label of time x-axis.}
  \item{ylab}{Label of energy y-axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{bty}{the type of box to be drawn around the energy plot.} 
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{The discrete version of the Teager-Kaiser operator is computed
  according to: \cr
  \deqn{y_{n} = x_{n}^{2/m} - (x_{n-M} \times x_{n+M})^{1/m}}{%
    y[n] = x[n]^{2/m} - (x[n-M]*x[n+M])^{1/m}}, \cr
  with \emph{M} the lag and \emph{m} the exponent which usually equal to 1 for
  the conventional operator. \cr
  The Teaser-Kaiser operator can be used to filter out a signal or to
  track amplitude modulations (AM) and/or frequency modulations (FM).\cr
  See examples.
}

\value{This function returns a two-column matrix, the first column is
  time and the second column includes the successive energy values. \cr
  m/2 NA values are added at the start and end of the vector.}

\references{Kvedalen, E. (2003). \emph{Signal processing using the Teager Energy
  Operator and other nonlinear operators}. University of Oslo,
Department of Informatics, PhD Thesis, x + 100 p.}

\author{Jerome Sueur}

\seealso{\code{\link{env}}, \code{\link{ifreq}}.}

\examples{
op <- par(mfrow=c(2,1))

## sinusoid AM 
s1 <- synth(f=8000, d=0.1, cf=200, am=c(100,10), output="Wave") 
oscillo(s1)
TKEO(s1)

## linear AM decrease
s2 <- synth(f=8000, d=0.1, cf=200, shape="decr", output="Wave") 
oscillo(s2)
TKEO(s2)

## sinusoid FM
s3 <- synth(f=8000, d=0.1, cf=200, fm=c(150,50,0), output="Wave") 
oscillo(s3)
TKEO(s3)

## linear FM increase
s4 <- synth(f=8000, d=0.1, cf=200, fm=c(0,0,600), output="Wave") 
oscillo(s4)
TKEO(s4)

## AM and FM
s5 <- synth(f=8000, d=0.1, cf=200, am=c(100,10), fm=c(150,50,0), output="Wave")
oscillo(s5)
TKEO(s5)

par(op)
}

\keyword{ts}
