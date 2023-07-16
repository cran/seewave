\name{wav2dBSPL}

\alias{wav2dBSPL}

\title{Convert wave amplitude to instantaneous sound pressure level (dB SPL)}

\description{This function applies an absolute dB SPL scale to the amplitude
of a time wave based on the properties of the recording chain, if
known.}

\usage{
wav2dBSPL(wave, f, channel = 1, gain, sensitivity = -35, Vadc = 2, pref = 2 * 10^-5)
}

\arguments{
  \item{wave}{an R object.}     
  \item{f}{sampling frequency of \code{wave} (in Hz). Does not need to
    be specified if embedded in \code{wave}.}
  \item{channel}{channel of the R object, by default left channel (1).}
  \item{gain}{total gain applied to the sound (preamplifer +
  amplifier), in dB.}
  \item{sensitivity}{sensitivity of the microphone, in dB/V (by default - 35 dB/V).}
  \item{Vadc}{maximal voltage (peak to peak) converted by the analog to
  digital convertor ADC, in V (by default 2 V).}
  \item{pref}{sound pressure reference in the medium, in Pa (by default =
  2*10^-5 Pa in air).}
}

\value{A numeric vector of the same legnth as \code{wave} containing the sound waveform in dB SPL (Sound Pressure level in dB)}

\note{This function can be used to turn a recorder, as a Songmeter or an
  Audiomoth, in a kind of soundlevel meter.}

\author{Sylvain Haupert}

\seealso{\code{\link{wav2leq}}}

\examples{
## note that the original gain, sensitivity and Vadc of tico
## are unknown, this is therefore a kind of 'unreal' example
data(tico)
res <- wav2dBSPL(tico, gain=18)
plot(x=seq(0,duration(tico),length=length(tico@left)), y=res, type="l", xlab="Time (s)",
ylab="dB SPL", las=1)
## dB SPL approximation for the complete sound
meandB(res)
}

\keyword{input}