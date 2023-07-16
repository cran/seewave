\name{oscilloEQ}
\alias{oscilloEQ}

\title{Oscillogram 'equalizer'}
\description{
  A multipanel plot of a time wave displaying the oscillogram of a bank
  of frequency filters like in an 'equalizer'.
}

\usage{
oscilloEQ(wave, f, channel = 1, flim = NULL, colwave = 1,
xlab = "Time (s)", ylab = "Frequency band (kHz)",
cexlab = 1, collab = 1, fontlab = 1,
savedir = ".", plot = TRUE, ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{channel}{channel of the R object, by default left channel (1).}
  \item{flim}{a numeric vector giving the ordered limites of the
  frequency filters to be applied. By default, 1 kHz frequency filters.
}
  \item{colwave}{colour of the oscillogram.}
  \item{xlab}{label of the x axis.}
  \item{ylab}{label of the y axis.}
  \item{cexlab}{character size for axes labels.}
  \item{collab}{color for axes labels.}
  \item{fontlab}{font for axes labels.}
  \item{savedir}{the path were the \code{.wav} files will be stored if
    \code{plot} is \code{FALSE}. By default, the working directory.}
  \item{plot}{a logical, if \code{TRUE} a plot is created, if
    \code{FALSE} then a series of \code{.wav} files are saved. Each file
    corresponds to a frequency band.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
  }
  
\details{
The function applies a bank of filters as delimited with the argument
\code{flim}.
If \code{plot} is \code{TRUE}, then the function displays the wave on a multiframe plot
so that the time*amplitude dynamics of each frequency filter can
be estimated. The filtered waves are generated using the function \code{\link{fir}}.
If \code{plot} is \code{FALSE}, then the corresponding waves are saved
as separated \code{.wav} file. Each file corresponds to a frequency
filter. 
}

\value{If \code{plot} is \code{FALSE} then a series of \code{.wav} files are saved. Each file
 corresponds to a frequency filter.
}

\author{Jerome Sueur}

\seealso{\code{\link{oscillo}}, \code{\link{oscilloST}}}

\examples{
data(peewit)
## default 1 kHz frequency filter
oscilloEQ(peewit)
## change de frequency filter limits
oscilloEQ(peewit, flim=c(0, 4, 8, 10))
oscilloEQ(peewit, flim=seq(2, 10, by=0.5))
## play with colors
oscilloEQ(peewit, colwave=c(1,2))
oscilloEQ(peewit, colwave=heat.colors)
blue.gray  <- colorRampPalette(c("darkblue", "lightgrey"))
oscilloEQ(peewit, colwave=blue.gray)
## save files instead of visualizing them
\dontrun{
oscilloEQ(peewit, plot=FALSE)
}
}

\keyword{dplot}
\keyword{ts}
