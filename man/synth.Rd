\name{synth}

\alias{synth}

\title{Synthesis of time wave (additive model)}

\description{
  This functions synthesizes pure or harmonic tone sound with amplitude modulation (am)
  and/or frequency modulation (fm).
}

\usage{
synth(f, d, cf, a = 1, signal = "sine", shape = NULL, p = 0,
am = c(0, 0, 0), fm = c(0, 0, 0, 0), harmonics = 1, 
plot = FALSE, listen = FALSE, output = "matrix",...)
}

\arguments{
  \item{f}{sampling frequency (in Hz).}
  \item{d}{duration (in s).}
  \item{cf}{carrier frequency (in Hz).}
  \item{a}{amplitude (linear scale, relative when adding different waves).}
  \item{signal}{a character vector specifying the shape of the signal,
    see \code{details}.}
  \item{shape}{modification of the whole amplitude shape of the wave,
    see \code{details}.}
  \item{p}{phase (in radians).}
  \item{am}{a numeric vector of length 3 describing amplitude modulation parameters,
    see \code{details}.}
  \item{fm}{a numeric vector of length 4 describing frequency modulation parameters,
    see \code{details}.}
  \item{harmonics}{a numeric specifying the number and the relative
    amplitude of harmonics, see \code{details}.}
  \item{plot}{if \code{TRUE} returns the spectrogram
    of the synthezised sound (by default \code{FALSE}).}
  \item{listen}{if \code{TRUE} the new sound is played back.}
  \item{output}{character string, the class of the object to return, either
    \code{"matrix"}, \code{"Wave"}, \code{"Sample"}, \code{"audioSample"} or \code{"ts"}.}
  \item{\dots}{other \code{\link{spectro}} graphical parameters.}
}

\details{

  \itemize{

    \item{\code{signal} is a character vector of length 1 that specifies
      the function used to synthesize the
      signal. There are three options:
      \enumerate{
	\item{"sine": for a sinus function}
	\item{"tria": for a triangle function}
	\item{"square": for a square function}
	\item{"saw": for a square function}
      }
    }
    
    \item{\code{shape} is a character vector of length 1 that allows
      to modify the whole amplitude shape of the wave. There
      are four options:
      \enumerate{
	\item{"incr": linear increase}
	\item{"decr": linear decrease} 
	\item{"sine": sinusoid-like shape}
	\item{"tria": triangular shape}
      }
    }

    \item{\code{am} is a numeric vector of length 3 including:
      \enumerate{
	\item{the amplitude modulation depth (in \%)}
	\item{the frequency of the amplitude modulation (in Hz),}
	\item{the phase of the amplitude modulation (in radian).}
      }
    }

    \item{\code{fm} is a numeric vector of length 3 including:
      \enumerate{
	\item{the maximum excursion of a sinusoidal frequency modulation (in
	  Hz),}
	\item{the frequency of a sinusoidal frequency modulation (in Hz),}
	\item{the maximum excursion of a linear frequency modulation (in
Hz).}
        \item{the phase of the frequency modulation (in radian).}
      }
    }

    \item{\code{harmonics} is a numeric vector that controls the number and the
      relative amplitude of harmonics synthesized. \cr
      By default \code{harmonics = 1} meaning that a pure tone made of a single
      harmonic (fundamental) will be produced.\cr
      To produce harmonics, the length of \code{harmonics} has to be
      greater than 1. The length of \code{harmonics} will set the number of
      harmonics, including the first one (fundamental). The value of
      each element of \code{harmonics} specify the relative ampltiude of
      each harmonic. The first value must equal to \code{1}. \cr
      Here are some examples:
      \itemize{
	\item{\code{harmonics = c(1, 0.5, 0.25)} will produce a sound
	  with three harmonics (fundamental + 2 harmonics), the second
	  harmonic having an amplitude half the fundamental ampltiude and the second
	  harmonic an amplitude a quarter of the fundamental amplitude.}
        \item{\code{harmonics = c(1, 0, 0.25)} will produce a sound with
	  two harmonics (fundamental + 1 harmonic) the second harmonic
	  having a null relative amplitude.}
	\item{\code{harmonics = rep(1,4)} will produce a sound with four
	  harmonics (fundamental + 3 harmonics) of equal amplitude.}
    }
  }
 }
}

\value{If \code{plot} is \code{FALSE}, a new wave is returned. The class
  of the returned object is set with the argument \code{output}.}

\references{Hartmann, W. M. 1998 \emph{Signals, sound and sensation}. New York: Springer.}

\author{Jerome Sueur and Laurent Lellouch.}

\seealso{\code{\link{synth2}}, \code{\link{noisew}}, \code{\link{pulse}}, \code{\link{echo}}}

\examples{
## You can use plot=TRUE and spectro() options
## to directly 'see' the new-built sounds
f <- 8000 # sampling frequency
d <- 1    # duration (1 s)
cf <- 440 # carrier frequecy (440 Hz, i.e. flat A tone)
# pure sinusoidal tone
s <- synth(f=f,d=d,cf=cf)
# pure triangular tone
s <- synth(f=f,d=d,cf=cf, signal="tria")
# pure tone with triangle overall shape
s <- synth(f=f,d=d,cf=cf,shape="tria")
# pure tones with am
s <- synth(f=f,d=d,cf=cf,am=c(50,10))
# pure tones with am
# and phase shift of pi radian (180 degrees)
s <- synth(f=f,d=d,cf=cf,am=c(50,10,pi))
# pure tone with +1000 Hz linear fm 
s <- synth(f=f,d=d,cf=cf,fm=c(0,0,1000))
# pure tone with sinusoidal fm
# (maximum excursion of 250 Hz, frequency of 10 Hz)
s <- synth(f=f,d=d,cf=cf,fm=c(250,10,0))
# pure tone with sinusoidal fm
# (maximum excursion of 250 Hz, frequency of 10 Hz,
# phase shift of pi radian (180 degrees))
s <- synth(f=f,d=d,cf=cf,fm=c(250,10,0, pi))
# pure tone with sinusoidal am
# (maximum excursion of 250 Hz, frequency of 10 Hz)
# and linear fm (maximum excursion of 500 Hz)
s <- synth(f=f,d=d,cf=cf,fm=c(250,10,500))
# the same with am
s <- synth(f=f,d=d,cf=cf,am=c(50,10), fm=c(250,10,250))
# the same with am and a triangular overall shape 
s <- synth(f=f,d=d,cf=cf,shape="tria",am=c(50,10), fm=c(250,10,250))
# an harmonic sound
s <- synth(f=f,d=d,cf=cf, harmonics=c(1, 0.5, 0.25))
# a clarinet-like sound
clarinet <- c(1, 0, 0.5, 0, 0.14, 0, 0.5, 0, 0.12, 0, 0.17)
s <- synth(f=f, d=d, cf = 235.5, harmonics=clarinet)
# inharmonic FM sound built 'manually'
fm <- c(250,5,0)
F1<-synth(f=f,d=d,cf=cf,fm=fm)
F2<-synth(f=f,d=d,a=0.8,cf=cf*2,fm=fm)
F3<-synth(f=f,d=d,a=0.6,cf=cf*3.5,fm=fm)
F4<-synth(f=f,d=d,a=0.4,cf=cf*6,fm=fm)
final1<-F1+F2+F3+F4
spectro(final1,f=f,wl=512,ovlp=75,scale=FALSE)
}

\keyword{datagen}
\keyword{ts}
