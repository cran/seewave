\name{zapsilw}

\alias{zapsilw}

\title{Zap silence periods of a time wave}

\description{
This function simply delete the silence periods of a time wave. 
}

\usage{zapsilw(wave, f, threshold = 5, plot = TRUE, Sample = FALSE, ...)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{threshold}{amplitude threshold (in \%) between silence and signal.}
  \item{plot}{logical, if \code{TRUE} plots the orginal and the new oscillograms
      (by default \code{TRUE}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.      
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{Use the argument \code{threshold} to set the level of silence. See
the \code{examples}.} 

\seealso{\code{\link{afilter}}, \code{\link{oscillo}}}

\examples{
data(orni)
zapsilw(orni,f=22050,colwave="red")
# setting the threshold value
zapsilw(orni,f=22050,threshold=1)
}

\keyword{ts}
