\name{read.audacity}

\alias{read.audacity}

\title{Audacity audio markers import}

\description{Read audio markers as exported by Audacity.}

\usage{read.audacity(file, format)}

\arguments{
  \item{file}{A .txt file produced by Audacity when exporting time or
    time x frequency markers.}
  \item{format}{The format of the file name that will appear in the
  value, that is in the first column of the data frame returned. if
  \code{"dir"} then the full path to the file is returned, if
  \code{"base"} only the base name of the file is returned.}
}

\details{Audacity opens the possibility to annotate sound files with a
  marker channel. These markers can be exported as .txt files. The
  function \code{read.audacity} import such .txt files whether they
  contain time markers or time x frequency markers. 
}

\value{A \code{data.frame}. The size of the \code{data.frame} differs
  whether the .txt file contains time markers or time x frequency
  markers.\cr
  For time markers, the \code{data.frame} contains 4 columns:
  \enumerate{
    \item{\code{file} returning the name of the input file either with
    the full path or with the base name only (see argument \code{format}),}
    \item{\code{label} the text label,}
    \item{\code{t1} the start time in seconds,}
    \item{\code{t2}the end time in seconds.}
  } 
  For time x frequency markers, the \code{data.frame} contains 6
  columns:
  \enumerate{
    \item{\code{file} returning the name of the input file either with
    the full path or with the base name only (see argument \code{format}),}
    \item{\code{label} the text label,}
    \item{\code{t1} the start time in seconds,}
    \item{\code{t2} the end time in seconds,}
    \item{\code{f1} the lower frequency in Hz,}
    \item{\code{f2} the upper frequency in Hz.}
  }
}

\references{Audacity is a free software
  distributed under the terms of the GNU General Public License.\cr
  Web site: \url{https://www.audacityteam.org/}
}

\author{Jerome Sueur}

\seealso{\code{\link{write.audacity}}}

\examples{
\dontrun{
## If 'markers.txt' is an export of Audacity markers
x <- read.audacity("markers.txt")
}
}

\keyword{input}
\keyword{data}
