\name{env}

\alias{env}

\title{Amplitude envelope of a time wave}

\description{This function returns the absolute or Hilbert
amplitude envelope of a time wave.}

\usage{env(wave, f, envt = "hil", msmooth = NULL,
ksmooth = NULL, plot = TRUE, k = 1, j = 1, ...)}

\arguments{
  \item{wave}{data describing a time wave
  or a \code{\link[sound]{Sample}} object generated loading a wav file
  with \code{\link[sound]{loadSample}} (package \pkg{sound}).}
  \item{f}{sampling frequency of \code{wave} (in Hz).
  Does not need to be specified if \code{wave} is a \code{\link[sound]{Sample}} object.}
  \item{envt}{the type of envelope to be returned: either "abs" for absolute
    amplitude envelope or "hil" for Hilbert amplitude envelope. See Details section.}
  \item{msmooth}{a vector of length 2 to smooth the amplitude envelope with a 
    mean sliding window. The first component is the window length
    (in number of points). The second component is the overlap between
    successive windows (in \%). See examples.}
  \item{ksmooth}{kernel smooth via \code{\link{kernel}}. See examples.}
  \item{plot}{logical, if \code{TRUE} returns a plot
      of \code{wave} envelope (by default \code{TRUE}).}
  \item{k}{number of horizontal sections when \code{plot} is \code{TRUE}
    (by default =1).}
  \item{j}{number of vertical sections when \code{plot} is \code{TRUE}
    (by default =1).}
  \item{\dots}{other \code{\link{oscillo}} graphical parameters.}
}

\details{
When \code{envt} is set as "abs", the amplitude envelope returned
is the absolute value of \code{wave}.\cr
When \code{envt} is set as "hil", the amplitude envelope returned is the modulus
(\code{\link{Mod}}) of the analytical signal of \code{wave}
obtained through the Hilbert transform (\code{\link{hilbert}}).
}

\value{Data are returned as one-column matrix when \code{plot} is \code{FALSE}.}

\author{Jerome Sueur \email{sueur@mnhn.fr}}

\note{Be aware that smoothing with either \code{msmooth} or \code{ksmooth}
changes the original number of points describing \code{wave}.}

\seealso{\code{\link{oscillo}},\code{\link{hilbert}}}

\examples{
data(tico)
# Hilbert amplitude envelope
env(tico,f=22050)
# absolute amplitude envelope
env(tico,f=22050,envt="abs")
# smoothing with a 10 points and 50% overlaping mean sliding window
env(tico,f=22050,msmooth=c(10,50))
# smoothing kernel
env(tico,f=22050,ksmooth=kernel("daniell",10))
# overplot of oscillographic and envelope representations
oscillo(tico,f=22050)
par(new=TRUE)
env(tico,f=22050,colwave=2)
}

\keyword{dplot}
\keyword{ts}