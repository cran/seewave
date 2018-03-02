\name{ggspectro}

\alias{ggspectro}

\title{2D-spectrogram of a time wave using ggplot2}

\description{This function returns a ggplot object to draw a spectrogram
  with the package ggplot2. This is an alternative to \code{\link{spectro}}.}

\usage{ggspectro(wave, f, tlab = "Time (s)",
flab = "Frequency (kHz)", alab = "Amplitude\n(dB)\n", ...)}

\arguments{
  \item{wave}{an R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{tlab}{label of the time axis.}
  \item{flab}{label of the frequency axis.}
  \item{alab}{label of the amplitude axis.}
  \item{\dots}{other non-graphical arguments to be passed to spectro
  (\code{wl}, \code{ovlp} etc).}
}

\details{This function return the fist layer (data and aesthetic
  mapping) of a ggplot2 plot. \cr
  See the example section to understand how to build a spectrogram and consult ggplot2 help to get what you
  exactly need.\cr
  There is no way to plot the oscillogram as \code{\link{spectro} does.}
}

\value{A ggpot layer.}

\references{Wickham H (2009) -- \emph{ggplot2: elegant graphics for data analysis}. UseR! Springer.}

\author{Jerome Sueur}

\note{This function requires \pkg{ggplot2} package.}

\seealso{\code{\link{spectro}}, \code{\link{spectro3D}}, \code{\link{dynspec}}}

\examples{
\dontrun{
require(ggplot2)
## first layer
v <- ggspectro(tico, ovlp=50)
summary(v)
## using geom_tile ##
v + geom_tile(aes(fill = amplitude)) + stat_contour()
## coordinates flip (interest?)
v + geom_tile(aes(fill = amplitude)) + stat_contour() + coord_flip()
## using stat_contour ##
# default (not nice at all)
v + stat_contour(geom="polygon", aes(fill=..level..))
# set up to 30 color levels with the argument bins
(vv <- v + stat_contour(geom="polygon", aes(fill=..level..), bins=30))
# change the limits of amplitude and NA values as transparent
vv + scale_fill_continuous(name="Amplitude\n(dB)\n", limits=c(-30,0), na.value="transparent")
# Black-and-white theme
(vv + scale_fill_continuous(name="Amplitude\n(dB)\n", limits=c(-30,0),
  na.value="transparent", low="white", high="black") + theme_bw())
# Other colour scale (close to spectro() default output)
v + stat_contour(geom="polygon", aes(fill=..level..), bins=30)
  + scale_fill_gradientn(name="Amplitude\n(dB)\n", limits=c(-30,0),
  na.value="transparent", colours = spectro.colors(30))
}}

\keyword{dplot}
\keyword{ts}
