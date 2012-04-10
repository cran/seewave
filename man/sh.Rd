\name{sh}

\alias{sh}

\title{Shannon and Renyi spectral entropy}

\description{This function computes the Shannon or Renyi entropy of a frequency
  spectrum}

\usage{sh(spec, alpha = NULL)}

\arguments{
  \item{spec}{a data set resulting of a spectral analysis obtained
    with \code{\link{spec}} or \code{\link{meanspec}} (not in dB).}
  \item{alpha}{by default \code{NULL} to compute Shannon entropy,
  otherwise superior to 0 but different to 1 to comptute Renyi entropy.}
}

\details{
  Shannon spectral entropy is calculated according to:\cr
  \deqn{S = -\frac{\sum_{i = 1}^{N} y_i log_2(y_i)}{log_2(N)}}{% 
    S = -sum(ylogy)/log(N)}
  with:\cr
  \emph{y} = relative amplitude of the \emph{i} frequency,\cr
  and
  \deqn{\sum_{i = 1}^{N} y_i = 1}{% 
    sum(y) = 1}\cr 
  and \emph{N} = number of frequencies.\cr
  Renyi spectral entropy of order alpha is calucalted according to:\cr
  \deqn{R = \frac{1}{1-\alpha} \times log_2(\sum_{i = 1}^{N} y_{i}^{\alpha})}{% 
    R = (1/(1-alpha)) * log2(sum(y^alpha))}
  with \deqn{\alpha \geq 0}{% 
    alpha >= 0}
  and \deqn{\alpha \neq 1}{% 
   alpha !=1}
 }


\value{A numeric vector of length 1 is returned.
  The value has no unit.}

\references{
  Han, NC, Muniandy SV, Dayou J (2011) Acoustic classification of
  Australian anurans based on hybrid spectral-entropy approach.  \emph{Applied
  Acoustics}. \cr

  Nunes RR, Almeida de MP, Sleigh JW (2004) -
  Spectral entropy: a new method for anesthetic adequacy.
  \emph{Revista Brasileira de Anestesiologia}, \bold{54}, 413-422.\cr

  Renyi A (1961) - On measures of information and entropy. Proceedings of the 4th Berkeley Symposium on Mathematics, Statistics and Probability 1960. pp. 547-561.
}

\author{Jerome Sueur}

\note{
  The Shannon spectral entropy of a noisy signal will tend towards 1 whereas
  the Shannon spectral entropy of a pure tone signal will tend towards
  0. See Han \emph{et al}. for details regarding the Renyi entropy.}

\seealso{\code{\link{csh}},\code{\link{th}}, \code{\link{H}}, \code{\link{sfm}}}

\examples{
a<-synth(f=8000,d=1,cf=2000,plot=FALSE)
speca<-spec(a,f=8000,at=0.5,plot=FALSE)
##########################
# Shannon spectral entropy
##########################
sh(speca)
# [1] 0.2336412
b<-noisew(d=1,f=8000)
specb<-spec(b,f=8000,at=0.5,plot=FALSE)
sh(specb)
# close to 1
##########################
# Renyi spectral entropy
##########################
sh(speca, alpha=2)
sh(speca, alpha=3)
}

\keyword{ts}

