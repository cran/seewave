\name{bstd}

\alias{bstd}

\alias{BSTD}

\title{Brigthness, Standard Deviation time plot}

\description{This function displays a plot inspired from original
timber BSTD plot}

\usage{
bstd(wave, f, channel = 1, wl = 512, ovlp = 0,
xlab = "Time (s)", ylab = "",
xlim = NULL, ylim = NULL,
palette = colorRamp(c("yellow", "red")), plot = TRUE, ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{channel}{channel of the R object, by default left channel (1).}
  \item{wl}{length of the window for the analysis (even number of points, by default = 512).}
  \item{ovlp}{overlap between two successive analysis windows (in \%).}
  \item{xlab}{title of the time X-axis.}
  \item{ylab}{title of the rms Y-axis.}
  \item{xlim}{range of time X-axis (in s).}
  \item{ylim}{range of rms scaled Y-axis (within [0, 1]).}
  \item{palette}{color palette using \code{\link{colorRamp}} function.}
  \item{plot}{logical, if \code{TRUE} plots BSTD grpahic (by default \code{TRUE}).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\value{The function returns a 4-column matrix:
    \item{time}{the time axis (X-axis)}
    \item{rms}{the rms value (Y-axis)}
    \item{sf}{the spectral flux (thickness)}
    \item{zcr}{the zero crossing rate (color)}
  }

\references{Malt M, Jourdan E (2015) Le “BSTD” – Une représentation
  graphique de la brillance et de l’écart type spectral, comme possible
  représentation de l’évolution du timbre sonore.
  L’analyse musicale aujourd’hui, Delatour, 2015, 9782752102430. hal-01265269}

\author{Jérôme Sueur}

\note{As originally proposed by Malt and Jourdan (2015), the BSTD was a
combination of the spectral centroid (height) to describe the timber brightness,
the spectral standard deviation (thickness) and the amplitude of the signal (color).
}


\seealso{\code{\link{oscillo}}}

\examples{
data(tico)
## default parameters
bstd(tico)
## change of the STFT parameters
bstd(tico, wl=256, ovlp=90)
## change of the color palette
bstd(tico, palette=colorRamp(c("yellow", "green", "red")), wl=256, ovlp=90) 
}

\keyword{dplot}
\keyword{ts}
