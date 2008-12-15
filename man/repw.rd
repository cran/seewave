\name{repw}

\alias{repw}

\title{Repeat a time wave}

\description{This function repeats a time wave}

\usage{repw(wave, f, times = 2, plot = FALSE, Sample = FALSE, ...)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}  \item{times}{a numeric of length 1 describing the number
  of times the wave has to be repeated.}
  \item{plot}{logical, if \code{TRUE} plots the repeated time wave}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

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
