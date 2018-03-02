\name{env}

\alias{env}

\title{Amplitude envelope of a time wave}

\description{This function returns the absolute or Hilbert
amplitude envelope of a time wave.}

\usage{env(wave, f, envt = "hil", 
msmooth = NULL, ksmooth = NULL, ssmooth = NULL,
fftw = FALSE, norm = FALSE,
plot = TRUE, k = 1, j = 1, ...)}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to be specified if embedded in \code{wave}.}
  \item{envt}{the type of envelope to be returned: either "abs" for absolute
    amplitude envelope or "hil" for Hilbert amplitude envelope. See Details section.}
  \item{msmooth}{a vector of length 2 to smooth the amplitude envelope with a 
    mean sliding window. The first component is the window length
    (in number of points). The second component is the overlap between
    successive windows (in \%). See examples.}
  \item{ksmooth}{kernel smooth via \code{\link{kernel}}. See examples.}
  \item{ssmooth}{sum smooth.}
  \item{fftw}{if \code{TRUE} calls the function \code{FFT} of the
    library \code{fftw} for faster computation for the Hilbert amplitude
  envelope (\code{envt="hil"}) and/or for kernell smoothing
  (\code{ksmooth}).}
  \item{norm}{a logical, if \code{TRUE} the amplitude of the envelope
		is normalised between 0 and 1.}
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

\author{Jerome Sueur. Implementation of 'fftw' argument by Jean Marchal and Francois Fabianek.}

\note{Be aware that smoothing with either \code{msmooth} or \code{ksmooth}
changes the original number of points describing \code{wave}.}

\seealso{\code{\link{oscillo}},\code{\link{hilbert}}}

\examples{
data(tico)
# Hilbert amplitude envelope
env(tico)
# absolute amplitude envelope
env(tico, envt="abs")
# smoothing with a 10 points and 50% overlaping mean sliding window
env(tico, msmooth=c(10,50))
# smoothing kernel
env(tico, ksmooth=kernel("daniell",10))
# sum smooth
env(tico, ssmooth=50)
# overplot of oscillographic and envelope representations
oscillo(tico)
par(new=TRUE)
env(tico, colwave=2)
}

\keyword{dplot}
\keyword{ts}
