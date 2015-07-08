\name{corenv}

\alias{corenv}

\title{Cross-correlation between two time wave envelopes}

\description{
  This function tests the similarity between two time wave envelopes by returning
  their maximal correlation and the time shift related to it. 
}
\usage{
corenv(wave1, wave2, f, envt="hil", msmooth = NULL, ksmooth = NULL,
plot = TRUE, plotval = TRUE,
method = "spearman", col = "black", colval = "red",
cexval = 1, fontval = 1, xlab = "Time (s)",
ylab = "Coefficient of correlation (r)", type = "l", pb = FALSE, ...)
}

\arguments{
  \item{wave1}{a first R object.}     
  \item{wave2}{a second R object.}
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{envt}{the type of envelope to be used: either "abs" for absolute
    amplitude envelope or "hil" for Hilbert amplitude envelope. See \code{\link{env}}.}
  \item{msmooth}{a vector of length 2 to smooth the amplitude envelope with a 
    mean sliding window. The first component is the window length
    (in number of points). The second component is the overlap between
    successive windows (in \%). See \code{\link{env}}.}
  \item{ksmooth}{kernel smooth via \code{\link{kernel}}. See \code{\link{env}}.}
  \item{plot}{logical, if \code{TRUE} plots r values against frequency shift
    (by default \code{TRUE}).}
  \item{plotval}{logical, if \code{TRUE} adds to the plot maximum r value
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
  Successive correlations between the envelopes of \code{wave1} and \code{wave2}
  are computed when regularly sliding forward and backward \code{wave2} along
  \code{wave1}.\cr
  The maximal correlation is obtained at a particular shift (time offset). 
  This shift may be positive or negative.\cr
  The higher \code{smooth} is set up,
  the faster will be the computation but less precise the results will be.\cr 
  The corresponding p value, obtained with \code{\link{cor.test}}, is  plotted.
  Inverting \code{wave1} and \code{wave2} may give slight different results.
}

\value{
  If \code{plot} is \code{FALSE}, \code{corenv} returns a list containing four
  components:
  \item{r}{a two-column matrix, the first colum corresponding to the time
    shift (frequency x-axis) and the second column corresponding to the successive
    r correlation values between \code{env1} and \code{env2} (correlation y-axis).}
  \item{rmax}{the maximum correlation value between \code{x} and \code{y}.}
  \item{p}{the p value corresponding to \code{rmax}.}
  \item{t}{the time offset corresponding to \code{rmax}.}
}

\author{Jerome Sueur}

\seealso{\code{\link{env}},\code{\link{corspec}},\code{\link{covspectro}},
  \code{\link{cor}},\code{\link{cor.test}}.}

\examples{
\dontrun{
data(orni)
# cross-correlation between two echemes of a cicada song
wave1<-cutw(orni,f=22050,from=0.3,to=0.4,plot=FALSE)
wave2<-cutw(orni,f=22050,from=0.58,to=0.68,plot=FALSE)
corenv(wave1,wave2,f=22050)
}
}

\keyword{dplot}
\keyword{ts}
