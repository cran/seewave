\name{combfilter}

\alias{combfilter}

\title{Comb filter}

\description{This function processes a feedforward comb filter and plots a spectrogram of the filtered wave asso-
ciated with the frequency response of the filter.}

\usage{combfilter(wave, f, alpha, K, units = c("samples", "seconds"),
plot = FALSE, output = "matrix", ...)}

\arguments{
  \item{wave}{an \code{R} object}
  \item{f}{sampling frequency (in Hz). Does not need to be specified if
  embedded in \code{wave}}
  \item{alpha}{a numeric vector of length 1 for the scaling factor. See Details.}
  \item{K}{a numeric vector of lenght 1 for the delay length, in
  \code{units}. See Details.}
  \item{units}{units in which \code{K} is given, the default is
  \code{'samples'} but can be set to \code{'seconds'}.}
  \item{plot}{a logical, if \code{TRUE} plots the spectrogram of the
  filtered \code{wave} and the frequency response of the comb filter.}
  \item{output}{character string, the class of the object to return,
  either \code{'matrix'}, \code{'Wave'}, \code{'Sample'},
  \code{'audioSample'} or \code{'ts'}.}
  \item{\dots}{other arguments to be passed to \code{\link{spectro}}
  except \code{scale} and \code{osc} that are set by default to \code{FALSE}.}
}

\details{A comb filter consists in adding a delayed version of a signal to itself resulting in constructive and destructive interference. The feedforward version of a comb filter can be written following:
  \deqn{y(n) = x(n) + \alpha \times x(n - K)}{%
    y(n) = x(n) + alpha * x(n - K)
  }
  where alpha is the scaling factor and K the dealy length.
  The frequency response of the filter is obtained with:
  \deqn{H(f) = \sqrt{(1+\alpha^2)+2 \times \cos(\omega K)}}{%
    H(f) = sqrt((1+alpha^2)+2*cos(omega*K))}
  The frequency response is periodic. The depth of the cycles is controlled with alpha and the number of cycles with K.
}

\value{A new wave is returned. The class of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur}

\note{Setting K to high values may generate unwanted results. \cr
  The feedback form of the combfilter is not implemented yet.}

\seealso{\code{\link{fir}}, \code{\link{squarefilter}}, \code{\link{drawfilter}}, \code{\link{ffilter}}}

\examples{
\dontrun{f <- 44100
## chirp
s1 <- synth(f=f, cf=1, d=2, fm=c(0,0,f/2), out="Wave")
combfilter(s1, alpha=1, K=50, plot=TRUE)
## harmonic sound
s2 <- synth(f=f, d=2, cf=600, harmonics=rep(1, 35), output="Wave")
combfilter(s2, alpha=1, K=10, plot=TRUE)
## noise, units in seconds
s3 <- noisew(d=2, f=44100, out="Wave")
combfilter(s3, alpha=0.5, K=1e-4, units="seconds", plot=TRUE)
}
}

\keyword{ts}
\keyword{filter}
