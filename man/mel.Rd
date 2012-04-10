\name{mel}

\alias{mel}

\title{Hertz / Mel conversion}

\description{This function converts Hertz data in Mel data.}

\usage{mel(x, inverse = FALSE)}

\arguments{
  \item{x}{a value in Hertz (or in Mel if \code{inverse} is \code{TRUE})}
  \item{inverse}{logical, if \code{TRUE} converts the Mel data in Hertz data.}
}

\details{
  Hertz to mel conversion is computed according to:\cr
  \deqn{m = 1127.01048 \times {\log{(1+(\frac{f}{700}))}}}{% 
    m = 1127.01048*log(1+(f/700))}
  with \emph{m} in Mel and \emph{f} in Hertz.

  Mel to Hertz conversion (when \code{inverse} is \code{TRUE}) 
  is therefore computed according to:\cr
  \deqn{f = 700 \times{(e^{\frac{m}{1127.01048}}-1)}}{% 
    f = 700*(exp(m/1127.01048)-1)}
  with \emph{f} in Hertz and \emph{m} in Mel.
}

\value{A corresponding \bold{R} object is returned.}

\references{Stevens, S. S., Volkman, J. and Newman, E. B. 1937.
  A scale for the measurement of psychological magnitude pitch.
  \emph{Journal of the Acoustical Society of America}, 8: 185-190.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{
  The Mel scale is a perceptual scale of pitches judged by listeners
  to be equal in distance from one another. The name Mel comes from the word
  melody to indicate that the scale is based on pitch comparisons. The reference
  point between this scale and normal frequency measurement is defined by
  equating a 1000 Hz tone, 40 dB above the listener's threshold,
  with a pitch of 1000 mels. 
}

\examples{
x<-seq(0,10000,by=50)
y<-mel(x)
plot(x,y,type="l",xlab = "f (hertz)", ylab = "f (mel)",
  main = "Mel scale", col="red")
}

\keyword{math}

