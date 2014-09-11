\name{acoustat}

\alias{acoustat}

\title{Statistics on time and frequency STFT contours}

\description{This function returns statistics based on STFT time and
  frequency contours.}

\usage{
acoustat(wave, f, wl = 512, ovlp = 0, wn = "hanning",
tlim = NULL, flim = NULL,
aggregate = sum, fraction = 90,
plot = TRUE, type = "l", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{window length for the analysis
    (even number of points) (by default = 512).}
  \item{ovlp}{overlap between two successive windows (in \%).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default
    \code{"hanning"}).}
  \item{tlim}{modifications of the time limits of the analysis (in s).}
  \item{flim}{modifications of the frequency limits of the analysis (in kHz).} 
  \item{aggregate}{a character vector of length 1 specifying the function to be applied on the rows (time) and columns (frequency) of the STFT matrix. By default set to \code{sum}.}
  \item{fraction}{a numeric vector of length 1, specifying a particular
    fraction of the contours amplitude to be captured by the initial and
    terminal percentile values (in \%). See details.}
  \item{plot}{a logical, if \code{TRUE} a two-frame plot is returned
    with the time and frequency contours and percentiles displayed.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn. See \code{\link{plot}} for details (by default "l" for lines).}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{The principle of \code{acoustat} is as follows:
  \enumerate{
    \item Compute the short-term Fourier transform (STFT) with usual
    parameters (\code{wl} for window length, \code{ovlp} for overlap of
    successive windows, and \code{wn} for the name of window shape).
    \item This results in a time * frequency matrix.
    \item Compute an aggregation function (specified with the argument
    \code{aggregate} set by default to \code{sum}) accross rows and
    columns of time * frequency matrix.)
    \item This results in two components: (\emph{i}) the time contour,
    and (\emph{ii}) the frequency contour.
    \item Each contour is considered as a probability mass function
    (PMF) and transformed into a cumulated distribution function (CDF).
    \item Measures are extracted from each CDF: median (M), initial
    percentile (P1) value, terminal percentile (P2) value, interpercentile range (IPR). P1, P2 and IPR are defined using a fraction parameter
    (\code{fraction}) that sets the percent of the contour amplitude to be captured by the initial and terminal percentile values. A fraction of 50\% would result in the familiar quartiles and interquartile range. An energy fraction of 80\% would return the 10th and 90th percentile values, and the width of the range in between.
  }
}


\value{
    The function returns a list with 10 items:
    \item{time.contour}{the time contour as a two-column matrix, the
    first colum being time (s) and the second colum being the amplitude
    probability mass function (no scale).}
    \item{freq.contour}{the frequency contour as a two-column matrix,
    the first colum being frequency (kHz) and the second colum being the
    amplitude probability mass function (no scale).}
    \item{time.P1}{the time initial percentile}
    \item{time.M}{the time median}
    \item{time.P2}{the time terminal percentile}
    \item{time.IPR}{the time interpercentile range}
    \item{freq.P1}{the frequency initial percentile}
    \item{freq.M}{the frequency median}
    \item{freq.P2}{the frequency terminal percentile}
    \item{freq.IPR}{the frequency interpercentile range}
  }

  \references{Fristrup, K. M. and Watkins, W. A. 1992. Characterizing
    acoustic features of marine animal sounds.
    Woods Hole Oceanographic Institution Technical Report WHOI-92-04.\cr
    Details of the principles can be found at \url{http://www.birds.cornell.edu/brp/research/algorithm/automated-and-robust-measurement-of-signal-features}}

  \author{Jerome Sueur}

  \note{\code{acoustat} was originally developped in Matlab language by Kurt Fristrup and XXXX Watkins (1992) .\cr
    The \acronym{R} function was kindly checked by Kurt Fristrup.}

  \seealso{\code{\link{meanspec}}, \code{\link{specprop}}}

  \examples{
data(tico)
note <- cutw(tico, from=0.5, to=0.9, output="Wave")
## default setting
acoustat(note)
## change the percentile fraction
acoustat(note, fraction=50)
## change the STFT parameters
acoustat(note, wl=1024, ovlp=80)
## change the function to compute the aggregate contours
## standard deviation instead of sum   
acoustat(note, aggregate=sd)
## direct time and frequency selection     
acoustat(tico, tlim=c(0.5,0.9), flim=c(3,6))
## some useless graphical changes
acoustat(note, type="o", col="blue") 
}

\keyword{ts}
