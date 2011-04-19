\name{notefreq}

\alias{notefreq}

\title{Frequency of a muscical note}

\description{This function computes the frequency of a musical note
  (Equal temperament)}

\usage{notefreq(note, ref = 440, octave = 3)}

\arguments{
  \item{note}{a numerical or a character vector. See \code{Note}.}
  \item{ref}{a numerical vector of length 1 for the reference frequency.}
  \item{octave}{a numerical vector of length for the octave number.}
}

\details{The frequency is computed according to:\cr
  \deqn{f = ref \times 2 ^{octave - 3 + \frac{note-10}{12}}}{% 
      f <- ref*2^((octave-3) + ((note-10)/12))}

  with:\cr
  \emph{ref} = reference frequency, \cr
  \emph{octave} = octave number, and \cr
  \emph{note} = rank of the note along the scale.
}

\value{The frequency in Hz is returned.}

\author{Jerome Sueur}

\note{The note can be given in two ways. The first solution is to give
  the rank of the note along the scale (e.g. rank 10 for A) or to give its names in
  characters with the following notation: C, D, E, F, G, A, B.}

\seealso{\code{\link{octaves}}}

\examples{
# Some notes frequency (use apply-like functions when dealing with character strings)
sapply(c("C", "A", "Gb"), notefreq)

# C major scale plot
n <- 1:12
freq <- notefreq(n)
names <- c("C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B") 
plot(n, freq, pch=19, cex=1.5,
     xlab = "Note name",
     ylab = "Frequency (Hz)",
     xaxt="n", las=1, main="Third octave")
axis(side=1, at=n, labels=names)
abline(h=freq, col="lightgrey")

# C major scale sound
f <- 2000 # sampling rate
s <- NULL
for (i in 1:length(freq))
  {
    tmp <- synth(d=0.5, f=f, cf=freq[i])
    s <- pastew(s, tmp, at="start", f)
  }
spectro(s, f, ovlp=75)
}

\keyword{maths}