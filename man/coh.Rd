\name{coh}

\alias{coh}

\title{Coherence between two time waves}

\description{
This function returns the frequency coherence between two time
waves.}

\usage{
coh(wave1, wave2, f, plot =TRUE, xlab = "Frequency (kHz)",
ylab = "Coherence", xlim = c(0,f/2000), type = "l",...)
}

\arguments{
\item{wave1}{data or a \code{\link[sound]{Sample}} object generated loading a wav
  file with \code{\link[sound]{loadSample}} (package \pkg{sound}) describing a first time wave.}
\item{wave2}{data or a \code{\link[sound]{Sample}} object generated loading a wav
  file with \code{\link[sound]{loadSample}} (package \pkg{sound}) describing a second time wave.}
\item{f}{sampling frequency of \code{wave1} and \code{wave1} (in Hz).
  Does not need to be specified if \code{wave1} and/or \code{wave2} are/is (a)
  \code{\link[sound]{Sample}} object(s).}
\item{plot}{logical, if \code{TRUE} plots the continuous coherence function
        (by default \code{TRUE}).}
\item{xlab}{title of the frequency X-axis.}
\item{ylab}{title of the coherence Y-axis.}
\item{xlim}{range of frequency X-axis.}
\item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
  See \code{\link{plot}} for details (by default "l" for lines).}
\item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{Coherence is a frequency domain function computed to show the degree of a relationship
between two signals. The value of the coherence function ranges between zero and one,
where a value of zero indicates there is no causal relationship between the signals.
A value of one indicates the existence of linear frequency response between the
two signals. This can be used, for instance, to compare the input and output
signals of a system.
}

\value{
If \code{plot} is \code{FALSE}, this function returns a one-column matrix
corresponding to the coherence.}

\author{Jérôme Sueur \email{sueur@mnhn.fr} but this function
is based on \code{\link{spec.pgram}}
by Martyn Plummer, Adrian Trapletti and B.D. Ripley.}

\note{This function is based on \code{\link{spec.pgram}}.}

\seealso{\code{\link{ccoh}}, \code{\link{spectro}}, \code{\link{spec.pgram}}.}


\examples{
wave1<-synth(d=1,f=4000,cf=500)
wave2<-synth(d=1,f=4000,cf=800)
ccoh(wave1,wave2,f=4000)
}

\keyword{dplot}
\keyword{ts}