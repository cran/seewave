\name{akamatsu}

\alias{akamatsu}

\title{Water tank minimum resonant and cutoff frequencies}

\description{This function computes the resonant and cutoff frequencies when
  recording in a given aquarium according to the criteria explained in Akamatsu et al. (2002)}

\usage{
akamatsu(Lx, Ly, Lz, mode = c(1,1,1),
         c = 148000,  plot = FALSE, xlab = "Frequency (kHz)",
         ylab = "Attenuation distance (cm)", ...)
}

\arguments{
  \item{Lx}{watertank length (in cm).}
  \item{Ly}{watertank width (in cm).}
  \item{Lz}{watertank height (in cm).}
  \item{mode}{mode, see details.}
  \item{c}{sound velocity in cm/s (by default 148000 cm/s in water).}
  \item{plot}{logical, if \code{TRUE} plots the attenuation distance in function of frequency.}
  \item{xlab}{title of the x axis if \code{plot} is \code{TRUE}.}
  \item{ylab}{title of the y axis if \code{plot} is \code{TRUE}.}
  \item{\dots}{other \code{\link{plot}} graphical parameters.}
}

\details{From Akamatsu et al. (2002): \cr \cr
  \emph{1. Resonant frequency} \cr \cr
  The calculated resonant frequencies of a rectangular glass tank with the dimension of Lx , Ly , and Lz (in centimeters) can
  be described by the following equation: \cr

  \deqn{f^{rectangular}_{lmn} = \frac{c}{2} \times \sqrt{\left(\frac{l}{L_{x}}\right)^2 +
      \left(\frac{m}{L_{y}}\right)^2  + \left(\frac{n}{L_{z}}\right)^2}}{% 
    f = 0.5*c*sqrt((l/Lx)^2 + (m/Ly)^2 + (n/Lz)^2)}

  where \emph{c} is the sound velocity (cm/s) and each \emph{l}, \emph{m}, \emph{n} reprents an
  integer, and the combination of these paramameters designates the
  'mode number'. The mode (1, 1, 1) represents the resonance wave of minimum
  frequency. The mode (2, 1, 1) represents one of the higher order of
  resonant component and has additional node of the soundpressure level
  at the middle of the X axis, \emph{i.e.}, Lx/2. \cr
  \cr
  \cr
  \emph{2. Cutoff frequency} \cr \cr
  The cutoff frequency can be calculated as follows: \cr
  \deqn{f^{rectangular}_{cutoff} = \frac{c}{2} \times \sqrt{
      \left(\frac{1}{L_{y}}\right)^2  + \left(\frac{1}{L_{z}}\right)^2}}{% 
    f = 0.5*c*sqrt((1/Ly)^2 + (1/Lz)^2)}
  \cr
  \cr
  \emph{3. Attenuation distance} \cr \cr
The theoretical attenuation distance \emph{D} can be expressed in function of the
cutoff frequency and the projected frequency following: \cr
  \deqn{D^{rectangular}(f) = 2 \times log_{10} \times \frac{c}{4 \pi
  f^{rectangular}_{cutoff}} \times \frac{1}{\sqrt{1-\left(\frac{f}{f^{rectangular}_{cutoff}}\right)^2}}}{% 
    f = (2*log(10)*c)/(4*pi*Fcut*sqrt(1-(f/Fcut)^2))}
}

\value{A list of two items:
  \item{res}{Resonant frequency (in Hz). See \code{Details}}
  \item{cut}{Cut frequency (in Hz). See \code{Details}}
}

\references{Akamatsu T, Okumura T, Novarini N, Yan HY (2002) Emprical refinements applicable to the recording of fish sounds in small tanks. \emph{Journal of the Acoustical Society of America}, 112, 3073-3082.}

\author{Camille Desjonqueres}

\examples{
akamatsu(60, 30, 40)
}

\keyword{ts}
