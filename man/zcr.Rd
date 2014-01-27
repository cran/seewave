\name{zcr}

\alias{zcr}

\title{Zero-crossing rate}

\description{This functions computes the zero-crossing rate of a time
  function, i. e. the average number the sign of a time wave changes.}

\usage{zcr(wave, f, wl = 512, ovlp = 0, plot = TRUE, type = "o", xlab =
"Time (s)", ylab = "Zero crossing rate", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{length of the window for the analysis (even number of
    points, by default = 512). If \code{NULL} the zero-crossing rate is
    computed of the complete signal.}
  \item{ovlp}{overlap between two successive analysis windows (in \%) if
  \code{wl} is not \code{NULL}.}
  \item{plot}{a logical, if \code{TRUE} plots a the zero-crossing rate
    results along time.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{xlab}{if \code{plot} is \code{TRUE}, label of the x axis.}
  \item{ylab}{if \code{plot} is \code{TRUE}, label of the y axis.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{The zero-crossing rate is computed according to: \cr
  \deqn{zcr = \frac{1}{2 \times N} \sum_{t=0}^{N-1} |sgn[x(t+1)] -sgn[x(t)]|}{%
  zcr = 0.5 * mean(abs(sign(x(t+1)) - sign(x(t))))
  }
  with:\cr
  \emph{N} the length of the signal \emph{x} \cr
  and where:
  \deqn{sgn[x(t)] = 1}{%
  sgn(x(t)) = 1}
  if \deqn{x(t) \leq 0}{%
  x(t) >= 0
  }
  and \deqn{sgn[x(t)] = -1}{%
  sgn(x(t)) = -1}
  if \deqn{x(t) > 0}{%
  x(t) > 0.  
    }
}

\note{The are two possibilities:
  \enumerate{
    \item if \code{wl} is \code{NULL} then the zero-crossing rate is
    computed for the complete signal.
    \item if \code{wl} is not \code{NULL} the the zero-crossing rate is
    computed for for a window sliding along the time wave.}
  The ZCR is supposed to help in detection of voiced/unvoiced sound sections.
}


\value{The are two possibilities:
  \enumerate{
    \item a numeric vector of length 1 if \code{wl} is \code{NULL},
    \item a numeric two-column matrix is returned with the first columnn being
    time (s) and the second colum being the zero-crossing rate (no scale) if \code{wl} is not \code{NULL}. 
  }
}

\references{\url{http://en.wikipedia.org/wiki/Zero-crossing_rate}}

\author{Jerome Sueur}


\seealso{\code{\link{zc}}}

\examples{
data(tico)
## a single value for the complete signal, no plot
zcr(tico, wl=NULL)
## a series of values computed for a sliding window oof 512 samples, plot
zcr(tico)
}

\keyword{ts}
