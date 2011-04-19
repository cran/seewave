\name{addsilw}

\alias{addsilw}

\title{Add or insert a silence section}

\description{Add or insert a silence section to a time wave.}

\usage{addsilw(wave, f, at = "end", choose = FALSE, d = NULL,
plot = FALSE, output = "matrix", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{at}{position where to add the silence section (in s).
  Can be also specified as \code{"start"}, \code{"middle"} or \code{"end"}.}
  \item{choose}{logical, if \code{TRUE} the point where silence will
	be added into \code{wave2} (=\code{at}) can be graphically chosen with a cursor.}
  \item{d}{duration of the silence section to add (in s).}
  \item{plot}{logical, if \code{TRUE} returns an oscillographic plot
  of \code{wave} with the new silence section (by default \code{TRUE}).}
  \item{output}{character string, the class of the object to return, either
  \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
of the returned object is set with the argument \code{output}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{oscillo}}, \code{\link{cutw}},\code{\link{deletew}},
\code{\link{fadew}},\code{\link{pastew}}, \code{\link{mutew}},\code{\link{revw}},
\code{\link{zapsilw}}
}

