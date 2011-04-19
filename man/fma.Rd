\name{fma}

\alias{fma}

\title{Frequency modulation analysis}

\description{
  This function computes the Fourier analysis of the instantaneous
  frequency of a time wave. This allows to detect periodicity 
  in frequency modulation.}

\usage{fma(wave, f, threshold = NULL, plot = TRUE, ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{threshold}{amplitude threshold for signal detection (in \% ).}
  \item{plot}{logical, if \code{TRUE} the spectrum of the instantaneous frequency 
    (by default \code{TRUE}).}
  \item{\dots}{other \code{\link{spec}} parameters.}
}

\details{
  This function is based on \code{ifreq} and \code{spec}.\cr 
  The instantaneous frequency of \code{wave} is first computed and the spectrum of this
  frequency modulation is then processed. All \code{env} and \code{spec} arguments can be
  set up.
}

\value{
  If \code{plot} is \code{FALSE}, \code{fma} returns a numeric vector
  corresponding to the computed spectrum.
  If \code{peaks} is not \code{NULL}, \code{fma} returns a list with
  two elements:
  \item{spec}{the spectrum computed}
  \item{peaks }{the peaks values (in kHz).}
}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{ifreq}}, \code{\link{hilbert}}, \code{\link{spec}}, \code{\link{ama}}}

\examples{
# a sound with a 1 Khz sinusoid FM
a<-synth(d=1, f=8000, cf=1500, fm=c(1000,1000,0))
fma(a,f=8000)
}

\keyword{dplot}
\keyword{ts}