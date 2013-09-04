\name{covspectro}

\alias{covspectro}

\title{Covariance between two spectrograms}

\description{
  This function tests the similarity between two spectrograms by returning
  their maximal covariance and the time shift related to it. 
}

\usage{
covspectro(wave1, wave2, f, wl = 512, wn = "hanning", n,
plot = TRUE, plotval = TRUE,
method = "spearman", col = "black", colval = "red", cexval = 1,
fontval = 1, xlab = "Time (s)",
ylab = "Normalised covariance (cov)", type = "l", pb = FALSE, ...)
}

\arguments{
  \item{wave1}{a first R object.}     
  \item{wave2}{a second R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{wl}{length of the window for the analysis
    (even number of points, by default = 512).}
  \item{wn}{window name, see \code{\link{ftwindow}} (by default \code{"hanning"}).}
  \item{n}{number of covariances computed between \code{wave1} and \code{wave2}
    when sliding \code{wave2} along \code{wave1}.}
  \item{plot}{logical, if \code{TRUE} plots r values against frequency shift
    (by default \code{TRUE}).}
  \item{plotval}{logical, if \code{TRUE} adds to the plot maximum R value
    and frequency offset (by default \code{TRUE}).}  
  \item{method}{a character string indicating which correlation coefficient is
    to be computed ("pearson", "spearman", or "kendall")
    (see \code{\link{cor}}).}
  \item{col}{colour of r values.}
  \item{colval}{colour of r max and frequency offset values.}
  \item{cexval}{character size of r max and frequency offset values.}
  \item{fontval}{font of r max and frequency offset values.}
  \item{xlab}{title of the frequency axis.}
  \item{ylab}{title of the r axis.}
  \item{type}{if \code{plot} is \code{TRUE}, type of plot that should be drawn.
    See \code{\link{plot}} for details (by default "l" for lines).}
  \item{pb}{if \code{TRUE} returns a text progress bar in the console.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
  Successive covariances between the spectrogram of \code{wave1} and
  the spectrogram of \code{wave2} are computed when regularly sliding
  forward and backward \code{wave2} along \code{wave1}.\cr
  The maximal covariance is obtained at a particular shift (time offset). 
  This shift may be positive or negative.\cr 
  \code{n} sets in how many steps \code{wave2} will be slided along \code{wave1}.
  Time process can be then decreased by setting low \code{n} value.\cr
  Inverting \code{wave1} and \code{wave2} may give slight different results.
}

\value{
  If \code{plot} is \code{FALSE}, \code{covspectro} returns a list containing
  three components:
  \item{cov}{the successive covariance values between \code{wave1} and \code{wave2}.}
  \item{covmax}{the maximum covariance between \code{wave1} and \code{wave2}.}
  \item{t}{the time offset corresponding to \code{cov}.}
}

\references{
  Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
    communication}. Springer, Berlin, Heidelberg.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\seealso{\code{\link{corspec}}, \code{\link{corenv}}, \code{\link{spectro}},
  \code{\link{cor}},}

\examples{
# covariance between two notes of a birdsong
\dontrun{
data(tico)
note1<-cutw(tico, f=22050, from=0.5, to=0.9)
note2<-cutw(tico, f=22050, from=0.9, to=1.3)
covspectro(note1,note2,f=22050,n=37)
}
}

\keyword{dplot}
\keyword{ts}
