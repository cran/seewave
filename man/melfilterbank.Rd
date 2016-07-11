\name{melfilterbank}

\alias{melfilterbank}

\title{Mel-filter bank for MFCC computation}

\description{This functions returns graphically and numerically the
  Mel-filters used to compute MFCC.}

\usage{
melfilterbank(f = 44100, wl = 1024,
minfreq = 0, maxfreq = f/2, m = 20,
palette, alpha = 0.5, plot = FALSE)
}

\arguments{
  \item{f}{sammpling frequency (in Hz).}
  \item{wl}{the Fourier window length (in number of samples).}
  \item{minfreq}{the minimum (or lower) frequency of the filter bank (in
  Hz).}
  \item{maxfreq}{the maximum (or upper) frequency of the filter bank (in Hz).}
  \item{m}{the total number of filters.}
  \item{palette}{an optional colour palette if \code{plot} is \code{TRUE}.}
  \item{alpha}{alpha-transparency when a colour palette is used.}
  \item{plot}{if \code{TRUE} all filters are displayed in a single plot.}
}

\value{A list of 3 items:
 \item{central.freq}{the kHz central frequencies of the filters,}
 \item{freq}{the kHz frequency scale,}
 \item{amp}{the amplitude of the filters, scaled between 0 and 1.}
}

\references{Sharan RV & Moir TJ (2016) Applications and advancements in
  automatic sound recognition. \emph{Neurocomputing}.}

\author{Jerome Sueur}

\note{These triangular filters are used for computing MFCCs.}

\seealso{\code{\link{mel}}}

\examples{
## default values
melfilterbank(plot=TRUE)
## with color surfaces
melfilterbank(palette=cm.colors, plot=TRUE)
## values changed
res <- melfilterbank(f=16000, wl=512, minfreq=300, plot=TRUE)
## plot the 1st filter only
plot(res$freq, res$amp[,1], type="l", xlab="Frequency (kHz)", ylab="Amplitude")
## plot the last filter only
plot(res$freq, res$amp[,ncol(res$amp)], type="l", xlab="Frequency (kHz)", ylab="Amplitude")
## get the kHz central frequencies of the succesive filters
res$central.freq
}

\keyword{ts}
\keyword{math}
