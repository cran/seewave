\name{wav2leq}

\alias{wav2leq}

\title{Convert wave amplitude to equivalent continuous sound pressure level (Leq)}

\description{This function computes the Leq value  of a time wave based
  on the properties of the recording chain, if known.}

\usage{
wav2leq(wave, f, channel = 1, gain, dt = 1, sensitivity = -35, Vadc = 2, pref = 2 * 10^-5)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{channel}{channel of the R object, by default left channel (1).}
  \item{gain}{total gain applied to the sound (preamplifer +
  amplifier), in dB.}
  \item{dt}{integration time step, in s (by default 1 s).}
  \item{sensitivity}{sensitivity of the microphone, in dB/V (by default - 35 dB/V).}
  \item{Vadc}{maximal voltage (peak to peak) converted by the analog to
  digital convertor ADC, in V (by default 2 V).}
  \item{pref}{sound pressure reference in the medium, in Pa (by default =
  2*10^-5 Pa in air).}
}

\note{This function can be used to turn a recorder, as a Songmeter or an
  Audiomoth, in a kind of soundlevel meter.}

\value{A numeric vector of length 1 returning the Leq value, in dB.}

\author{Sylvain Haupert}

\seealso{\code{\link{wav2dBSPL}}}

\examples{
data(tico)
## with a 1 s time of integration (that is for the first time window)
wav2leq(tico, gain=18)
## with a 0.5 s time of integration (that is for 3 successive time windows)
wav2leq(tico, dt=0.5, gain=18)
## for the complete sound
wav2leq(tico, dt=duration(tico), gain=18)
## which is equivalent to
meandB(wav2dBSPL(tico, gain=18))
}

\keyword{input}
