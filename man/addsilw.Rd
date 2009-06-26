\name{addsilw}

\alias{addsilw}

\title{Add or insert a silence section}

\description{Add or insert a silence section to a time wave.}

\usage{addsilw(wave, f, at = "end", choose = FALSE, d = NULL,
plot = FALSE, Sample = FALSE,...)}

\arguments{
	\item{wave}{a \code{vector}, a \code{matrix} (first column),
	an object of class \code{ts}, \code{\link[sound]{Sample}} (left channel),
	or \code{\link[tuneR]{Wave}} (left channel).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is an object of class \code{ts},
	\code{\link[sound]{Sample}}, or \code{\link[tuneR]{Wave}}.}
  \item{at}{position where to add the silence section (in s).
  Can be also specified as \code{"start"}, \code{"middle"} or \code{"end"}.}
  \item{choose}{logical, if \code{TRUE} the point where silence will
	be added into \code{wave2} (=\code{at}) can be graphically chosen with a cursor.}
  \item{d}{duration of the silence section to add (in s).}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot
  of \code{wave} with the new silence section (by default \code{TRUE}).}
  \item{Sample}{if \code{TRUE} and \code{plot} is \code{FALSE}
  returns an object of class \code{\link[sound]{Sample}}}.
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned as a one-column matrix
or as a \code{\link[sound]{Sample}} object if \code{Sample} is \code{TRUE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{cutw}},\code{\link{deletew}},
\code{\link{fadew}},\code{\link{pastew}}, \code{\link{mutew}},\code{\link{revw}},
\code{\link{zapsilw}}
}

\examples{
data(tico)
addsilw(tico,f=22050,d=0.2)
addsilw(tico,f=22050,at="end",d=0.5)
addsilw(tico,f=22050,at=0.33,d=0.46)
}

\keyword{dplot}
\keyword{ts}

