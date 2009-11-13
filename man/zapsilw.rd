\name{zapsilw}

\alias{zapsilw}

\title{Zap silence periods of a time wave}

\description{
  This function simply deletes the silence periods of a time wave. 
}

\usage{zapsilw(wave, f, threshold = 5, plot = TRUE, output = "matrix", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{threshold}{amplitude threshold (in \%) between silence and signal.}
  \item{plot}{logical, if \code{TRUE} plots the orginal and the new oscillograms
    (by default \code{TRUE}).}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

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
