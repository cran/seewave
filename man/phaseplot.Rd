\name{phaseplot}

\alias{phaseplot}

\title{Phase-phase 2D or 3D plot of a time wave}

\description{This function returns a 2D or 3D representation of a time
  wave according to its first, second and possibly third derivatives.}

\usage{phaseplot(wave, f, dim = 3, plot = TRUE, type = "l",
xlab = "1st derivative",
ylab = "2nd derivative",
zlab = "3rd derivative", ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{dim}{a vector of lenght 1, the number of dimensions of the
  plot. Can be either 2 or 3.}
  \item{plot}{logical, if \code{TRUE} plots phase-phase plot (by default \code{TRUE}).}
  \item{type}{type of plot that should be drawn. See \code{\link{plot}} for details (by default "l" for lines).}
  \item{xlab}{title of the x axis.}
  \item{ylab}{title of the y axis.}
  \item{zlab}{title of the z axis.}
  \item{\dots}{other \code{\link{plot}} or \code{plot3d} graphical
  parameters of the package \pkg{rgl}.}
}

\value{If \code{plot} is \code{FALSE} then a 2 or 3 column matrix is
  returned. The position of the column is related to the order of the
  derivative (\emph{i}. \emph{e}. first colum = first derivative.}

\references{For use of such plots see:
  Rice AN, Land BR, Bass AH (2011) - Nonlinear acoustic complexity in a
  fish 'two-voice' system. \emph{Proceedings of the Royal Society B}, in press.
}

\note{Phase-phase plot can be used to test non-linearity.}

\author{Jerome Sueur}

\examples{
\dontrun{
require(rgl)
data(tico)
phaseplot(tico)
}}

\keyword{dplot}
\keyword{ts}
