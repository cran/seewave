\name{fpeaks}

\alias{fpeaks}

\title{Frequency peak detection}

\description{This function searches for peaks of a frequency spectrum.}

\usage{
fpeaks(spec, f = NULL,
nmax = NULL, amp = NULL, freq = NULL, threshold = NULL,
mel =FALSE,
plot = TRUE, title = TRUE,
xlab = NULL, ylab = "Amplitude",
labels = TRUE, legend = TRUE, collab = "red", ...)
}

\arguments{
  \item{spec}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}}. Can be in dB.}
  \item{f}{sampling frequency of \code{spec} (in Hz). Not requested if
  the first column of \code{spec} contains the frequency axis.}
  \item{nmax}{maximal number of peaks detected. Overrides \code{amp},
    \code{freq} and \code{threshold}. See details.}
  \item{amp}{amplitude slope parameter, a numeric vector of length
    2. Refers to the amplitude slopes of the peak. The first value is the
    left slope and the second value is the right slope. Only peaks with
    higher slopes than threshold values will be kept. See details.}
  \item{freq}{frequency threshold parameter (in Hz). If the frequency
    difference of two successive peaks is less than this threshold, then 
    the peak of highest amplitude will be kept only. See details.}
  \item{threshold}{amplitude threshold parameter. Only peaks above this
    threshold will be considered. See details.}
  \item{mel}{a logical, if \code{TRUE} the (htk-)mel scale is used.}
  \item{plot}{logical, if \code{TRUE}, a plot showing the peaks is returned.}
  \item{title}{logical, if \code{TRUE} add the number of peaks detected
    as a plot title.}
  \item{xlab}{label of the x-axis.}
  \item{ylab}{label of the y-axis.}
  \item{labels}{logical, if \code{TRUE} peak labels are plotted.}
  \item{legend}{logical, if \code{TRUE} a legend returning the different
  selection parameters (\code{nmax}, \code{amp}, \code{freq},
  \code{threshold}, \code{threshold}) is added to the plot.}
  \item{collab}{labels color.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{Here are some details regarding the different selection
  parameters:\cr
  \itemize{
  \item \code{nmax}: this parameter is to be used if you wish to get a
  specific number of peaks. The peaks selected are those with the
  highest slopes. It then does not work in conjunction with the other parameters.
  \item \code{freq}: this parameter allows to remove from the selection
  successive peaks with a small frequency difference. Imagine you have two
  successive peaks at 1200 Hz and 1210 Hz and at 0.5 
  and 0.25 in amplitude. If you set \code{freq} to 50 Hz, then only the first
  peak will be kept.
  \item \code{amp}: this parameter allows to remove from the selection
  peaks with low slopes. You can make the selection on both slopes or on
  a single one. Imagine you have an asymetric peak with a 0.01 left slope and a
  0.02 right slope. The peak will be discarded for the following
  settings: both values higher than 0.02 (\emph{e.g.} \code{amp =
    c(0.03,0.04)}), the first value higher than 0.01 (\emph{e.g.} \code{amp =
    c(0.02,0.001)}), the second value higher than 0.02 (\emph{e.g.} \code{amp =
    c(0.001,0.03)}). If you do not want apply the selection
  on one of the slope use 0. For instance, a selection on the left slope
  only will be achieved with: \code{amp = c(0.02,0)}.
  \item \code{threshold}: this parameter can be used to do a rough
  selection on the spectrum. Peaks with an amplitude value (not a slope)
  lower than this threshold will be automatically discarded. This can be
  useful when you want to remove peaks of a low-amplitude background noise.
  }
}

\value{A two-column matrix, the first column corresponding to the
  frequency values (x-axis) and the second column corresponding to the
  amplitude values (y-axis) of the peaks.}

\note{You can also use \code{fpeaks} with other kind of spectrum, for
  instance a cepstral spectrum. See examples.}

\author{Jerome Sueur and Amandine Gasc}

\seealso{\code{\link{localpeaks}}, \code{\link{meanspec}}, \code{\link{spec}}}

\examples{
data(tico)
spec <- meanspec(tico, f=22050, plot=FALSE)
specdB <- meanspec(tico, f=22050, dB="max0", plot=FALSE)
# all peaks
fpeaks(spec)
# 10 highest peaks
fpeaks(spec, nmax=10)
# highest peak (ie dominant frequency)
fpeaks(spec, nmax=1)
# peaks that are separated by more than 500 Hz
fpeaks(spec, freq=500)
# peaks with a left slope higher than 0.1
fpeaks(spec, amp=c(0.1,0))
# peaks with a right slope higher than 0.1
fpeaks(spec, amp=c(0,0.1))
# peaks with left and right slopes higher than 0.1
fpeaks(spec, amp=c(0.1,0.1))
# peaks above a 0.5 threshold
fpeaks(spec, threshold=0.5)
# peaks of a dB spectrum with peaks showing slopes higher than 3 dB
fpeaks(specdB, amp=c(3,3))
# comparing different parameter settings
meanspec(tico, f=22050)
col <- c("#ff000090","#0000ff75","#00ff00")
cex <- c(2,1.25,1.5)
pch <- c(19,17,4)
title(main="Peak detection \n (spectrum with values between 0 and 1)")
res1 <- fpeaks(spec, plot = FALSE)
res2 <- fpeaks(spec, amp=c(0.02,0.02), plot =FALSE)
res3 <- fpeaks(spec, amp=c(0.02,0.02), freq=200, plot = FALSE)
points(res1, pch=pch[1], col=col[1], cex=cex[1])
points(res2, pch=pch[2], col=col[2], cex=cex[2])
points(res3, pch=pch[3], col=col[3], cex=cex[3])
legend("topright", legend=c("all peaks","amp", "amp & freq"), pch=pch,
pt.cex=cex, col=col, bty="n")
# example with a cepstral spectrum
data(sheep)
res <- ceps(sheep,f=8000,at=0.4,wl=1024,plot=FALSE)
fpeaks(res, nmax=4, xlab="Quefrency (s)")
# melscale
require(tuneR)
mel <- melfcc(sheep, nbands = 256, dcttype = "t3", fbtype = "htkmel", spec_out=TRUE)
melspec.mean <- apply(mel$aspectrum, MARGIN=2, FUN=mean)
melspec.mean <- melspec.mean/max(melspec.mean) # [0,1] scaling 
fpeaks(melspec.mean, nmax=4, f=8000, mel=TRUE)
fpeaks(melspec.mean, freq=4, f=8000, mel=TRUE) # freq in Hz!
fpeaks(melspec.mean, threshold=0.3, f=8000, mel=TRUE)
fpeaks(melspec.mean, amp=c(0.1,0.1), f=8000, mel=TRUE)
}

\keyword{dplot}
\keyword{ts}
