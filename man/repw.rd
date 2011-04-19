\name{repw}

\alias{repw}

\title{Repeat a time wave}

\description{This function repeats a time wave}

\usage{repw(wave, f, times = 2, plot = FALSE, output= "matrix", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{times}{a numeric of length 1 describing the number
    of times the wave has to be repeated.}
  \item{plot}{logical, if \code{TRUE} plots the repeated time wave.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{addsilw}}, \code{\link{cutw}},
  \code{\link{deletew}}, \code{\link{fadew}}, \code{\link{mutew}},
  \code{\link{pastew}}, \code{\link{revw}}, \code{\link{zapsilw}}}

\examples{
data(tico)
repw(tico,f=22050,plot=TRUE)
}

\keyword{dplot}

\keyword{ts}
