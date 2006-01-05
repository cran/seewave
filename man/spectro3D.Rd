\name{spectro3D}

\alias{spectro3D}

\title{3D-spectrogram of a time wave}

\description{
This functions returns a three-dimension spectrographic representation of a
time wave. The function corresponds to short-term Fourier transform.
}

\usage{
spectro3D(wave, f, wl, zp = 0, ovlp = 0, plot = TRUE,
magt = 10, magf = 10, maga = 2, palette = rev.terrain.colors)
}

\arguments{
  \item{wave}{data describing a time wave.}
  \item{f}{sampling frequency of \code{wave}.}
  \item{wl}{length of the window for the analysis
    (even number of points).}
  \item{zp}{zero-padding (even number of points), see \code{Details}.}
  \item{ovlp}{overlap between two successive windows (in \% ).}
  \item{plot}{logical, if \code{TRUE} plots the spectrogram
        (by default \code{TRUE}).}
  \item{magt}{magnification of the time axis.}
  \item{magf}{magnification of the frequency axis.}
  \item{maga}{magnification of the amplitude axis.}
  \item{palette}{a color palette function to be used to assign colors in
      the plot, see \code{Details}.}
}

\details{
A Hanning function is applied to the analysis window.\cr
\code{zp} adds 0 values on both sides of the analysis window. This increases
frequency resolution without altering time resolution.\cr
Any colour palette can be used. In particular, it is possible to use other palettes coming with
\pkg{seewave}: \code{rev.gray.colors.1}, \code{rev.gray.colors.2},
\code{spectro.colors}, \code{rev.heat.colors}, 
\code{rev.cm.colors}, \code{rev.topo.colors}, corresponding to the reverse
of \code{heat.colors},\code{topo.colors},
\code{cm.colors}.\cr
Use \code{magt}, \code{magf} and \code{maga} to resize the plot.   
}

\value{
If \code{plot} is \code{FALSE}, this function returns a matrix. Each column
corresponds to a Fourier transform of length \code{wl}.}

\author{Jérôme Sueur \email{jerome.sueur@univ-tours.fr} and 
Caroline Simonis-Sueur \email{csimonis@mnhn.fr}.}

\note{This function runs with \pkg{rgl} and is based on \code{\link{fft}}.
See examples of \code{\link{spectro}} for analysis arguments
(\code{wl},\code{zp}, \code{ovlp}).}

\seealso{\code{\link{spectro}}, \code{\link{fft}}.}

\examples{
data(tico)
spectro3D(tico,f=22050,wl=512,ovlp=75,zp=16,maga=4,palette=rev.terrain.colors)
}

\keyword{dplot}
\keyword{ts}
