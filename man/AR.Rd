\name{AR}

\alias{AR}

\title{Acoustic Richness index}

\description{This function computes the Acoustic Richness index based on
M and Ht indices}

\usage{AR(..., datatype = "objects", envt = "hil",
msmooth = NULL, ksmooth = NULL, ssmooth = NULL,
pattern = "[wav]$|[mp3]$")}

\arguments{
\item{\dots}{\code{Wave}, \code{WaveMC}, \code{audioSample} objects if \code{datatype="objects"}, or a path as a character string to a directory including \code{.wav} and/or \code{.mp3} files if \code{datatype="files"}.}
\item{datatype}{A character string to specify if inputs are either \code{R} objects (\code{datatype="objects"}, default) or files (\code{datatype="files"}).}
\item{envt}{the type of envelope to be returned: either \code{"abs"} for absolute amplitude envelope or \code{"hil"} for Hilbert (default) amplitude envelope. See \code{\link{env}}.}
\item{msmooth}{mean smooth. See \code{\link{env}}.}
\item{ksmooth}{kernel smooth via kernel. See \code{\link{env}}.}
\item{ssmooth}{sum smooth. See \code{\link{env}}.}
\item{pattern}{an optional regular expression. Only file names which match the regular expression will be returned when \code{datatype="files"}. By default \code{.wav} or \code{.mp3} files. See \code{\link{dir}}.}
}

\details{AR is ranked index based on the rank of the M and Ht indices obtained with the functions \code{\link{M}} and \code{\link{th}} respectively following:
  \deqn{AR = \frac{rank(M) \times rank(H_{t})}{n^2}}{%
      AR = rank(M)*rank(Ht) / n^2}
    with
    \deqn{0 \leq AR \leq 1}{%
      0 <= AR <= 1}
}

\value{A \code{data.frame} with three columns (M, Ht, AR) and n columns, with n the number of objects (respectively files) used as input.}

\references{Depraetere M, Pavoine S, Jiguet F, Gasc A, Duvail S, Sueur J (2012) Monitoring animal diversity using acoustic indices: implementation in a temperate woodland. \emph{Ecological Indicators}, \bold{13}, 46-54.}

\author{Jerome Sueur and Marion Depraetere}

\note{As a ranked index, the results returned by AR strongly depends with the set of objects (respectively files) used as   input. Comparaison between different data sets may be spurious. Computing AR on a set of a single object does not make any sense but is allowed.}

\seealso{\code{\link{M}}, \code{\link{th}}, \code{\link{env}}}

\examples{
## input as R objects
data(orni)
data(tico)
AR(orni, tico)
## give names to objects if you wish to have them as row names of the returned data.frame
AR(orni=orni, tico=tico)
## input as files, here stored in the working directory
## file names will be used as row names of the returned data.frame
\dontrun{
require(tuneR)
AR(getwd(), datatype="files")
}
}

\keyword{ts}
\keyword{index}