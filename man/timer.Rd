\name{timer}

\alias{timer}

\title{Time measurements of a time wave}

\description{
  This function computes and shows the duration of signal periods, 
  pause periods and their ratio.
}

\usage{
timer(wave, f, threshold = 5, dmin = NULL, envt="abs",
power = 1, msmooth = NULL, ksmooth = NULL,
ssmooth = NULL, asmooth=NULL, tlim = NULL, plot = TRUE, plotthreshold = TRUE,
col = "black", colval = "red", 
xlab = "Time (s)", ylab = "Amplitude", ...)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{threshold}{amplitude threshold for signal detection (in \%), or
  alternatively a function to be applied on the waveform scaled between
  0 and 1. See examples.}
  \item{dmin}{time threshold (minimum duration) for signal detection (in s).}  
  \item{envt}{the type of envelope to be used: either "abs" for absolute
    amplitude envelope or "hil" for Hilbert amplitude envelope. See \code{\link{env}}.}
  \item{power}{a power factor applied to the amplitude
    envelope. Increasing \code{power} will reduce low amplitude
    modulations and increase high amplide modulations. This can be used
    to reduce background noise (by default equals to 1, \emph{i.e.} no change.}
  \item{msmooth}{a vector of length 2 to smooth the amplitude envelope with a 
    mean sliding window. The first component is the window length
    (in number of points). The second component is the overlap between
    successive windows (in \%).  See \code{\link{env}}.}
  \item{ksmooth}{kernel smooth for the amplitude enveloppe via \code{\link{kernel}}. See \code{\link{env}}.}
  \item{ssmooth}{sum smooth for the amplitude enveloppe. See \code{\link{env}}.}
  \item{asmooth}{autocorrelation smooth for the amplitude enveloppe. See \code{\link{env}}.}
  \item{tlim}{modifications of the time X-axis limits.}
  \item{plot}{logical, if \code{TRUE} plots the envelope and the measurements 
    (by default \code{TRUE}).}
  \item{plotthreshold}{logical, if \code{TRUE} plots the threshold
    as an horizontal line on the graph (by default \code{TRUE}).}  
  \item{col}{colour of the envelope.}
  \item{colval}{colour of plotted measurements.}
  \item{xlab}{title of the x-axis.}
  \item{ylab}{title of the y-axis.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\value{
  A list containing seven items:
  \item{s}{duration of signal period(s) in seconds}
  \item{p}{duration of pause period(s) in seconds}
  \item{r}{ratio between the signal and silence periods(s)}
  \item{positions}{a list containing four elements:}
  \item{s.start}{start position(s) of signal period(s)}
  \item{s.end}{end position(s) of signal period(s)}
  \item{first}{whether the first event detected is a pause or a signal}
  }

\section{Warning}{Setting to high values to \code{msmooth} or
  \code{ssmooth} might return inaccurate results. Double check your
  results if so.}

\author{Jerome Sueur}

\seealso{\code{\link{env}}, \code{\link{cutw}}, \code{\link{pastew}}.}
