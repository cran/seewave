\name{preemphasis}

\alias{preemphasis}

\title{Pre-emphasis speech filter}

\description{A pre-emphasis frequency filter for speech}

\usage{preemphasis(wave, f, alpha = 0.9,
plot = FALSE, output = "matrix", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{alpha}{time constant, see \code{Details}.}
  \item{plot}{a logical, if \code{TRUE} plots the spectrogram of the
  filtered \code{wave} and the frequency response of the comb filter.}
  \item{output}{character string, the class of the object to return,
  either \code{'matrix'}, \code{'Wave'}, \code{'Sample'},
  \code{'audioSample'} or \code{'ts'}.}
  \item{\dots}{other arguments to be passed to \code{\link{spectro}}
  except \code{scale} and \code{osc} that are set by default to \code{FALSE}.}
}

\details{The function applies a pre-emphasis filter usually applied in
  speech analysis. The filter is a kind of high-pass frequency filter
  that amplifies the high-frequency content of the sample. The filter is
  defined with:
 \deqn{y(n) = x(n) - \alpha \times x(n - 1)}{%
    y(n) = x(n) - alpha * x(n - 1)
  }
  where alpha is a time constant usually set between 0.9 and 1.
 
  The frequency response of the filter is obtained with:
  \deqn{H(f)=1+a^{2}-2 \times \alpha \times \cos(2 \times \pi \times
  f/f_{s})}{%
    H(f) = 1+a^2 - 2*a*cos(2*pi*freq*1000/f)}

}

\value{A new wave is returned. The class of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur}

\seealso{\code{\link{bwfilter}}, \code{\link{combfilter}}, \code{\link{ffilter}}, \code{\link{fir}},\code{\link{lfs}}, \code{\link{afilter}}}

\examples{
data(sheep)
fc <- 150
f <- sheep@samp.rate
alpha <- exp(-2*pi*fc/f)
res <- preemphasis(sheep, alpha=alpha, output="Wave")
}

\keyword{ts}
\keyword{filter}
