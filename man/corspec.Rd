\name{corspec}

\alias{corspec}

\title{Cross-correlation between two frequency spectra}

\description{
This function tests the similarity between two frequency spectra by returning
their maximal correlation and the frequency shift related to it. 
}

\usage{
corspec(x, y, range, plot = TRUE, plotval = TRUE,
method = "spearman", col = "black", colval = "red",
cexval = 1, fontval = 1, xlab = "Frequency (kHz)",
ylab = "Coefficient of correlation (r)", ...)
}

\arguments{
  \item{x}{a first data set resulting of a spectral analysis obtained
      with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}
  \item{y}{a second data set resulting of a spectral analysis obtained
      with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}
  \item{range}{range of \code{x} and \code{y} (in kHz).}
  \item{plot}{logical, if \code{TRUE} plots r values against frequency shift (by default \code{TRUE}).}
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
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{
It is important not to have data in dB.\cr
Successive correlations between \code{x} and \code{y} are computed when regularly
shifting \code{y} towards lower or higher frequencies.\cr
The maximal correlation is obtained at a particular shift (frequency offset). 
This shift may be positive or negative.\cr
The corresponding p value, obtained with \code{\link{cor.test}}, is plotted.\cr
Inverting \code{x} and \code{y} may give slight different results, see examples.
}

\value{
If \code{plot} is \code{FALSE}, \code{corspec} returns a list containing four
     components:
  \item{r}{the successive correlation values between \code{x} and \code{y}.}
  \item{rmax}{the maximum correlation value between \code{x} and \code{y}.}
  \item{p}{the p value corresponding to \code{rmax}.}
  \item{f}{the frequency offset corresponding to \code{rmax}.}
}

\references{
Hopp, S. L., Owren, M. J. and Evans, C. S. (Eds) 1998. \emph{Animal acoustic 
communication}. Springer, Berlin, Heidelberg.}

\author{Jérôme Sueur \email{jerome.sueur@ibaic.u-psud.fr}}

\seealso{\code{\link{spec}}, \code{\link{meanspec}}, \code{\link{corspec}},
        \code{\link{covspectro}}, \code{\link{cor}}, \code{\link{cor.test}}.}

\examples{
data(tico)
# compare the two first notes spectra
a<-spec(tico,f=22050,wl=512,at=0.2,plot=FALSE)
c<-spec(tico,f=22050,wl=512,at=1.1,plot=FALSE)
op<-par(mfrow=c(2,1), mar=c(4.5,4,3,1))
spec(tico,f=22050,wl=512,at=0.2,col="blue",type="l")
par(new=TRUE)
spec(tico,f=22050,wl=512,at=1.1,col="green",type="l")
legend(x=8,y=0.5,c("Note A", "Note C"),lty=1,col=c("blue","green"),bty="o")
par(mar=c(5,4,2,1))
corspec(a,c,range=c(0,11.025),type="l",
  ylim=c(-0.25,0.8),xaxs="i",yaxs="i",las=1)
par(op)
# different correlation methods give different results...
op<-par(mfrow=c(3,1))
corspec(a,c,range=c(0,11.025),
  type="l",xaxs="i",las=1, ylim=c(-0.25,0.8))
title("spearmann correlation (by default)")
corspec(a,c,range=c(0,11.025),
  type="l",xaxs="i",las=1,ylim=c(0,1),method="pearson")
title("pearson correlation")
corspec(a,c,range=c(0,11.025),
  type="l",xaxs="i",las=1,ylim=c(-0.23,0.5),method="kendall")
title("kendall correlation")
par(op)
# inverting x and y does not give exactly similar results
op<-par(mfrow=c(2,1),mar=c(2,4,3,1))
corspec(a,c,range=c(0,11.025),type="l")
corspec(c,a,range=c(0,11.025),type="l")
par(op)
}

\keyword{dplot}
\keyword{ts}