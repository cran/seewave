pkgname <- "seewave"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('seewave')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("H")
### * H

flush(stderr()); flush(stdout())

### Name: H
### Title: Total entropy
### Aliases: H
### Keywords: ts

### ** Examples

data(orni)
H(orni,f=22050)
# changing the spectral parameter (wl)
H(orni,f=22050,wl=1024)
# changing the temporal parameter (msmooth)
H(orni,f=22050,msmooth=c(20,0))



cleanEx()
nameEx("Q")
### * Q

flush(stderr()); flush(stdout())

### Name: Q
### Title: Resonance quality factor of a frequency spectrum
### Aliases: Q
### Keywords: dplot ts

### ** Examples

# bird song
data(tico)
t<-spec(tico,f=22050,at=1.1,plot=FALSE,dB="max0")
op<-par(mfrow=c(2,1),las=1)
Q(t,type="l")
Q(t,type="l",xlim=c(3.8,4.2),ylim=c(-60,0))
title("zoom in")
par(op)
# cricket, changing the dB level
data(pellucens)
p<-spec(pellucens,f=11025,at=0.5,plot=FALSE,dB="max0")
op<-par(mfrow=c(3,1))
Q(p,type="l",xlim=c(1.8,2.6),ylim=c(-70,0))
title("level = - 3 (default value)",col.main="red")
Q(p,type="l",level=-6,
    xlim=c(1.8,2.6),ylim=c(-70,0),colval="blue")
title("level = - 6",col.main="blue")
Q(p,type="l",level=-9,
    xlim=c(1.8,2.6),ylim=c(-70,0),colval="green")
title("level = - 9",col.main="green")
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("afilter")
### * afilter

flush(stderr()); flush(stdout())

### Name: afilter
### Title: Amplitude filter
### Aliases: afilter
### Keywords: ts

### ** Examples

data(orni)
op<-par(mfrow=c(2,1))
afilter(orni,f=22050)
title(main = "threshold level = 5")
afilter(orni,f=22050,threshold=0.5,colwave="blue")
title(main = "threshold level = 0.5")
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ama")
### * ama

flush(stderr()); flush(stdout())

### Name: ama
### Title: Amplitude modulation analysis of a time wave
### Aliases: ama
### Keywords: dplot ts

### ** Examples

data(orni)
# detection of 2 main amplitude modulations in a cicada song:
# one with a 0.020 kHz frequency (due to signal/silence periodicity)
# one with a 0.258 kHz frequency (due to pulses in the echemes)
# one with a 2.369 kHz frequency (fundamental frequency)
ama(orni,f=22050,wl=1024)
# these amplitude modulations can be identify with a cursor:
ama(orni,f=22050,wl=1024,identify=TRUE)



cleanEx()
nameEx("attenuation")
### * attenuation

flush(stderr()); flush(stdout())

### Name: attenuation
### Title: Generate sound intensity attenuation data
### Aliases: attenuation
### Keywords: data

### ** Examples

# theoretical attenuation up to 150 m of a 100 dB/1m sound source
attenuation(lref=100,dstop=150,n=200)



cleanEx()
nameEx("autoc")
### * autoc

flush(stderr()); flush(stdout())

### Name: autoc
### Title: Short-term autocorrelation of a time wave
### Aliases: autoc
### Keywords: dplot ts

### ** Examples

data(sheep)
# fundamental frequency of a sheep
res<-autoc(sheep, f=8000, threshold=5, fmin=100, fmax=700)
spectro(sheep, f=8000, ovlp=75, scale=FALSE)
points(res, pch=20)
legend(0.5, 3.6, "Fundamental frequency", pch=20, bty=0, cex=0.7)



cleanEx()
nameEx("ccoh")
### * ccoh

flush(stderr()); flush(stdout())

### Name: ccoh
### Title: Continuous coherence function between two time waves
### Aliases: ccoh
### Keywords: dplot ts

### ** Examples

wave1<-synth(d=1,f=4000,cf=500)
wave2<-synth(d=1,f=4000,cf=800)
ccoh(wave1,wave2,f=4000)



cleanEx()
nameEx("ceps")
### * ceps

flush(stderr()); flush(stdout())

### Name: ceps
### Title: Cepstrum or real cepstrum
### Aliases: ceps
### Keywords: dplot ts

### ** Examples

data(sheep)
ceps(sheep,f=8000,at=0.4,wl=1024)



cleanEx()
nameEx("cepstro")
### * cepstro

flush(stderr()); flush(stdout())

### Name: cepstro
### Title: 2D-cepstrogram of a time wave
### Aliases: cepstro
### Keywords: dplot ts

### ** Examples

data(sheep)
sheepc <- cutw(sheep, f=8000, from = 0.19, to = 2.3)
cepstro(sheepc,f=8000)



cleanEx()
nameEx("coh")
### * coh

flush(stderr()); flush(stdout())

### Name: coh
### Title: Coherence between two time waves
### Aliases: coh
### Keywords: dplot ts

### ** Examples

wave1<-synth(d=1,f=4000,cf=500)
wave2<-synth(d=1,f=4000,cf=800)
coh(wave1,wave2,f=4000)



cleanEx()
nameEx("convSPL")
### * convSPL

flush(stderr()); flush(stdout())

### Name: convSPL
### Title: Convert sound pressure level in other units
### Aliases: convSPL
### Keywords: math

### ** Examples

# conversion of two SPL measurements taken at 0.5 m from the source
convSPL(c(80,85),d=0.5) 



cleanEx()
nameEx("corenv")
### * corenv

flush(stderr()); flush(stdout())

### Name: corenv
### Title: Cross-correlation between two time wave envelopes
### Aliases: corenv
### Keywords: dplot ts

### ** Examples

data(orni)
# cross-correlation between two echemes of a cicada song
wave1<-cutw(orni,f=22050,from=0.3,to=0.4,plot=FALSE)
wave2<-cutw(orni,f=22050,from=0.58,to=0.68,plot=FALSE)
corenv(wave1,wave2,f=22050)



cleanEx()
nameEx("corspec")
### * corspec

flush(stderr()); flush(stdout())

### Name: corspec
### Title: Cross-correlation between two frequency spectra
### Aliases: corspec
### Keywords: dplot ts

### ** Examples

data(tico)
# compare the two first notes spectra
a<-spec(tico,f=22050,wl=512,at=0.2,plot=FALSE)
c<-spec(tico,f=22050,wl=512,at=1.1,plot=FALSE)
op<-par(mfrow=c(2,1), mar=c(4.5,4,3,1))
spec(tico,f=22050,at=0.2,col="blue")
par(new=TRUE)
spec(tico,f=22050,at=1.1,col="green")
legend(x=8,y=0.5,c("Note A", "Note C"),lty=1,col=c("blue","green"),bty="o")
par(mar=c(5,4,2,1))
corspec(a,c, ylim=c(-0.25,0.8),xaxs="i",yaxs="i",las=1)
par(op)
# different correlation methods give different results...
op<-par(mfrow=c(3,1))
corspec(a,c,xaxs="i",las=1, ylim=c(-0.25,0.8))
title("spearmann correlation (by default)")
corspec(a,c,xaxs="i",las=1,ylim=c(0,1),method="pearson")
title("pearson correlation")
corspec(a,c,xaxs="i",las=1,ylim=c(-0.23,0.5),method="kendall")
title("kendall correlation")
par(op)
# inverting x and y does not give exactly similar results
op<-par(mfrow=c(2,1),mar=c(2,4,3,1))
corspec(a,c)
corspec(c,a)
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("covspectro")
### * covspectro

flush(stderr()); flush(stdout())

### Name: covspectro
### Title: Covariance between two spectrograms
### Aliases: covspectro
### Keywords: dplot ts

### ** Examples

# covariance between two notes of a birdsong
data(tico)
note1<-cutw(tico, f=22050, from=0.5, to=0.9)
note2<-cutw(tico, f=22050, from=0.9, to=1.3)
covspectro(note1,note2,f=22050,n=37)



cleanEx()
nameEx("crest")
### * crest

flush(stderr()); flush(stdout())

### Name: crest
### Title: Crest factor and visualization
### Aliases: crest
### Keywords: ts

### ** Examples

data(tico)
crest(tico, f=22050)
# see the crest location and change the default graphical parameters
crest(tico, f=22050, plot=TRUE, sym="-")



cleanEx()
nameEx("csh")
### * csh

flush(stderr()); flush(stdout())

### Name: csh
### Title: Continuous spectral entropy
### Aliases: csh
### Keywords: ts

### ** Examples

data(orni)
csh(orni,f=22050,wl=512,ovlp=50)
# using the threshold argument can lead to some edge effets
# here sh=1 at the end of echemes
csh(orni,f=22050,wl=512,ovlp=50,threshold=5)



cleanEx()
nameEx("cutspec")
### * cutspec

flush(stderr()); flush(stdout())

### Name: cutspec
### Title: Cut a frequency spectrum
### Aliases: cutspec
### Keywords: ts

### ** Examples

data(orni)
a<-meanspec(orni,f=22050,plot=FALSE)
b<-cutspec(a,flim=c(4,8))
# quick check with a plot
plot(b,type="l")
# effects on spectral entropy
sfm(a)
sfm(b)



cleanEx()
nameEx("cutw")
### * cutw

flush(stderr()); flush(stdout())

### Name: cutw
### Title: Cut a section of a time wave
### Aliases: cutw
### Keywords: dplot ts

### ** Examples

# a 0.4 s section in a bird song
data(tico)
a<-cutw(tico,f=22050,from=0.5,to=0.9)
oscillo(a,22050)
# a direct way to see what has been cut
cutw(tico,f=22050,from=0.5,to=0.9,plot=TRUE)



cleanEx()
nameEx("dBscale")
### * dBscale

flush(stderr()); flush(stdout())

### Name: dBscale
### Title: dB colour scale for a spectrogram display
### Aliases: dBscale
### Keywords: dplot ts

### ** Examples

data(pellucens)
# place the scale on the left and not on the right as spectro() does
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1, 2), nc = 2), widths = c(1, 5))
par(mar=c(5,3,4,2))
dBscale(collevels=seq(-30,0,1),side=2)
par(mar=c(5,4,4,2))
spectro(pellucens, f=22050,wl=512,scale=FALSE)
par(def.par)
# place the scale on the top and not on the right as spectro() does
def.par <- par(no.readonly = TRUE)
layout(matrix(c(0,1,2,2), nc = 2, byrow=TRUE),widths=c(1,2),heights=(c(1,5.5)))
par(mar=c(0.5,3,4,2))
dBscale(collevels=seq(-30,0,1), textlab = "",side=3)
mtext("Amplitude (dB)",side=2,line = 1,at=0.6,cex=0.75)
par(mar=c(5,4,0.5,2))
spectro(pellucens, f=22050,wl=512,scale=FALSE)
par(def.par)   



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("dBweight")
### * dBweight

flush(stderr()); flush(stdout())

### Name: dBweight
### Title: dB weightings
### Aliases: dBweight
### Keywords: ts math

### ** Examples

# weight for a 50 Hz frequency
dBweight(f=50)
# A weight for the 1/3 Octave centre frequencies.
dBweight(f=c(20,25,31.5,40,50,63,80,100,125,160,200,250,
315,400,500,630,800,1000,1500,
1600,2000,2500,3150,4000,5000,
6300,8000,10000,12500,16000,20000))$A
# correction for a 50 Hz sound emitted at 100 dB
dBweight(f=50, dB=100)
# weighting curves plot
f <- seq(10,20000,by=10)
par(las=1)
plot(f, dBweight(f)$A, type="n", log="x",
xlim=c(10,10^5),ylim=c(-80,20),xlab="",ylab="",xaxt="n",yaxt="n")
abline(v=c(seq(10,100,by=10),seq(100,1000,by=100),
seq(1000,10000,by=1000),seq(10000,100000,by=10000),
c(100,1000,10000,100000)),col="lightgrey",lty=2)
abline(v=c(100,1000,10000,100000),col="grey")
abline(h=seq(-80, 20, 20),col="grey")
par(new=TRUE)
plot(f, dBweight(f)$A, type="l", log="x",
xlab="Frequency (Hz)", ylab="dB",lwd=2, col="blue", xlim=c(10,10^5),ylim=c(-80,20))
title(main="Acoustic weighting curves (10 Hz -20 kHz)")
lines(x=f, y=dBweight(f)$B, col="green",lwd=2)
lines(x=f, y=dBweight(f)$C, col="red",lwd=2)
lines(x=f, y=dBweight(f)$D, col="black",lwd=2)
legend("bottomright",legend=c("dB(A)","dB(B)","dB(C)","dB(D)"),
lwd=2,col=c("blue","green","red","black"),bty="o",bg="white")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("deletew")
### * deletew

flush(stderr()); flush(stdout())

### Name: deletew
### Title: Delete a section of a time wave
### Aliases: deletew
### Keywords: dplot ts

### ** Examples

# deletion a 0.4 s section in a bird song
data(tico)
a<-deletew(tico,f=22050,from=0.5,to=0.9)
oscillo(a,22050)
# a direct way to see what has been cut
deletew(tico,f=22050,from=0.5,to=0.9,plot=TRUE)



cleanEx()
nameEx("dfreq")
### * dfreq

flush(stderr()); flush(stdout())

### Name: dfreq
### Title: Dominant frequency of a time wave
### Aliases: dfreq
### Keywords: dplot ts

### ** Examples

data(tico)
f <- 22050
# default
dfreq(tico,f)
# using the amplitude threshold and changing the graphical output
dfreq(tico, f, ovlp=50,threshold=5, type="l", col=2)
# using 'at' argument for specific positions along the time axis
dfreq(tico, f, at=c(0.25, 0.75, 1.2, 1.6))
dfreq(tico, f, at=seq(0.5, 1.4, by=0.005), threshold=5)
# a specific number of measures on a single note
dfreq(tico, f, at=seq(0.5, 0.9, len=100), threshold=5, xlim=c(0.5,0.9))
# overlap on spectrogram
# and use of 'clip' argument to better track the dominant frequency
# in noisy conditions
op <- par()
ticon <- tico@left/max(tico@left) + noisew(d=length(tico@left)/f, f)
spectro(ticon, f)
res <- dfreq(tico, f, clip=0.2, plot=FALSE)
points(res, col=2, pch =13)
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("diffenv")
### * diffenv

flush(stderr()); flush(stdout())

### Name: diffenv
### Title: Difference between two amplitude envelopes
### Aliases: diffenv
### Keywords: dplot ts

### ** Examples

data(tico) ; tico <- tico@left
data(orni) ; orni <- orni@left
# selection in tico of two waves with similar duration
tico2<-tico[1:length(orni)]
diffenv(tico2,orni,f=22050,plot=TRUE)
# smoothing the envelope gives a better graph but slightly changes the result
diffenv(tico2,orni,f=22050,msmooth=c(20,0),plot=TRUE)



cleanEx()
nameEx("diffspec")
### * diffspec

flush(stderr()); flush(stdout())

### Name: diffspec
### Title: Difference between two frequency spectra
### Aliases: diffspec
### Keywords: dplot ts

### ** Examples

a<-noisew(f=8000,d=1)
b<-synth(f=8000,d=1,cf=2000)
c<-synth(f=8000,d=1,cf=1000)
d<-noisew(f=8000,d=1)
speca<-spec(a,f=8000,wl=512,at=0.5,plot=FALSE)
specb<-spec(b,f=8000,wl=512,at=0.5,plot=FALSE)
specc<-spec(c,f=8000,wl=512,at=0.5,plot=FALSE)
specd<-spec(d,f=8000,wl=512,at=0.5,plot=FALSE)
diffspec(speca,speca,f=8000)
#[1] 0 => similar spectra of course !
diffspec(speca,specb)
diffspec(speca,specc,plot=TRUE)
diffspec(specb,specc,plot=TRUE)
diffspec(speca,specd,plot=TRUE)



cleanEx()
nameEx("diffwave")
### * diffwave

flush(stderr()); flush(stdout())

### Name: diffwave
### Title: Difference between two time waves
### Aliases: diffwave
### Keywords: ts

### ** Examples

data(tico) ; tico <- tico@left
data(orni) ; orni <- orni@left
# selection in tico to have two waves of similar duration (length)
tico <- tico[1:length(orni)]
diffwave(tico,orni,f=22050)
# changing the frequency parameter (wl)
diffwave(tico,orni,f=22050,wl=1024)
# changing the temporal parameter (msmooth)
diffwave(tico,orni,f=22050,msmooth=c(20,0))



cleanEx()
nameEx("discrets")
### * discrets

flush(stderr()); flush(stdout())

### Name: discrets
### Title: Time series discretisation
### Aliases: discrets
### Keywords: ts

### ** Examples

# a random variable
discrets(rnorm(30))
discrets(rnorm(30),symb=3)
# a frequency spectrum
data(tico)
spec1<-spec(tico,f=22050,at=0.2,plot=FALSE)
discrets(spec1[,2])



cleanEx()
nameEx("drawenv")
### * drawenv

flush(stderr()); flush(stdout())

### Name: drawenv
### Title: Draw the amplitude envelope of a time wave
### Aliases: drawenv
### Keywords: datagen ts

### ** Examples

a<-synth(d=1,f=22050,cf=1000)
# drawenv(a,f=22050,plot=TRUE)
# choose points on the oscillogram view to draw a new enveloppe
# stop (ESC on Windows; right mouse button on Linux)
# check the result on the second graphics device opened thanks to plot=TRUE



cleanEx()
nameEx("dynoscillo")
### * dynoscillo

flush(stderr()); flush(stdout())

### Name: dynoscillo
### Title: Dynamic oscillogram
### Aliases: dynoscillo
### Keywords: dplot ts

### ** Examples

require(rpanel)
data(tico)
dynoscillo(tico, wn=4)



cleanEx()
nameEx("dynspec")
### * dynspec

flush(stderr()); flush(stdout())

### Name: dynspec
### Title: Dynamic sliding spectrum
### Aliases: dynspec
### Keywords: dplot ts

### ** Examples

data(sheep)
require(rpanel)
dynspec(sheep,f=8000,wl=1024,ovlp=50,osc=TRUE)
dev.off()



cleanEx()
nameEx("echo")
### * echo

flush(stderr()); flush(stdout())

### Name: echo
### Title: Echo generator
### Aliases: echo
### Keywords: datagen ts

### ** Examples

# generation of the input wave
a<-synth(f=11025,d=1,cf=2000,shape="tria",am=c(50,10),fm=c(1000,10,1000))
# generation of three echoes
# with respectively a relative amplitude of 0.8, 0.4, and 0.2
# and with a delay of 1s, 2s, and 3s  from the beginning of the input wave
aecho<-echo(a,f=11025,amp=c(0.8,0.4,0.2),delay=c(1,2,3))
# oscillographic output to see what we have generated
op<-par(mfrow=c(2,1))
oscillo(a,f=11025,title="Input signal")
oscillo(aecho,f=11025,colwave="blue",title="Signal with echoes",coltitle="blue")
par(op)
# another echo with time delays overlapping with the input wave
echo(a,f=11025,amp=c(0.4,0.2,0.4),delay=c(0.6,0.8,1.5),plot=TRUE,listen=TRUE)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("env")
### * env

flush(stderr()); flush(stdout())

### Name: env
### Title: Amplitude envelope of a time wave
### Aliases: env
### Keywords: dplot ts

### ** Examples

data(tico)
# Hilbert amplitude envelope
env(tico,f=22050)
# absolute amplitude envelope
env(tico,f=22050,envt="abs")
# smoothing with a 10 points and 50% overlaping mean sliding window
env(tico,f=22050,msmooth=c(10,50))
# smoothing kernel
env(tico,f=22050,ksmooth=kernel("daniell",10))
# sum smooth
env(tico,f=22050,ssmooth=50)
# overplot of oscillographic and envelope representations
oscillo(tico,f=22050)
par(new=TRUE)
env(tico,f=22050,colwave=2)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("export")
### * export

flush(stderr()); flush(stdout())

### Name: export
### Title: Export sound data
### Aliases: export
### Keywords: IO

### ** Examples

a<-synth(f=8000,d=2,cf=2000,plot=FALSE)
export(a,f=8000)
unlink("a.txt")



cleanEx()
nameEx("fadew")
### * fadew

flush(stderr()); flush(stdout())

### Name: fadew
### Title: Fade in and fade out of a time wave
### Aliases: fadew
### Keywords: dplot ts

### ** Examples

a<-noisew(d=5,f=4000)
op<-par(mfrow=c(3,1))
fadew(a,f=4000,din=1,dout=2,plot=TRUE,title="Linear",cexlab=0.8)
fadew(a,f=4000,din=1,dout=2,shape="exp",plot=TRUE,title="Exponential shape",
    colwave="blue",coltitle="blue",cexlab=0.8)
fadew(a,f=4000,din=1,dout=2,shape="cos",plot=TRUE,title="Cosinus-like shape",
    colwave="red",coltitle="red",cexlab=0.8)
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("fbands")
### * fbands

flush(stderr()); flush(stdout())

### Name: fbands
### Title: Frequency bands plot
### Aliases: fbands
### Keywords: dplot ts

### ** Examples

data(sheep)
spec <- meanspec(sheep, f=8000, plot=FALSE)
# default plot
fbands(spec)
# setting a specific number of bands
fbands(spec, bands=6)
#setting specific regular bands limits
fbands(spec, bands=seq(0,4,by=0.25))
# some plot tuning
op <- par(las=1)
fbands(spec, bands=seq(0,4,by=0.1),
       horiz=TRUE, col=heat.colors(41),
       xlab="", ylab="",
       cex.axis=0.75, cex.names = 0.75,
       axes=FALSE)
par(op)
# showing or not the width of the bands
oct <- octaves(440,3)/1000
op <- par(mfrow=c(2,1))
fbands(spec, bands=oct, col="blue")
fbands(spec, bands=oct, width = TRUE, col="red")
par(op)
# kind of horizontal zoom
op <- par(mfrow=c(2,1))
fbands(spec, bands=seq(0,4,by=0.2), col=c(rep(1,10),
   rep("orange",5),rep(1,5)), main="all frequency range")
fbands(spec, bands=seq(2,3,by=0.2),
   col="orange", main="a subset or zoom in")
par(op)
# kind of dynamic frequency bands
specs <- dynspec(sheep, f=8000, plot= FALSE)$amp
out <- apply(specs, f=8000, MARGIN=2, FUN = fbands, bands = seq(0,4,by=0.2), col = 1, ylim=c(0,max(specs))) 



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("fdoppler")
### * fdoppler

flush(stderr()); flush(stdout())

### Name: fdoppler
### Title: Doppler effect
### Aliases: fdoppler
### Keywords: math

### ** Examples

# a 400 Hz source moving toward or away from the observer at 85 m/s
fdoppler(f=400,vs=85)
# [1] 533.3333
fdoppler(f=400,vs=85,movs="away")
# [1] 320
# use wasp() if you wish to have exact sound speed at a specific temperature
fdoppler(f=wasp(f=400,t=25)$c, vs=85)
# [1] 461.8667
# Doppler effect at different source speeds
f<-seq(1,10,by=1); lf<-length(f)
v<-seq(10,300,by=20); lv<-length(v)
res<-matrix(numeric(lf*lv),ncol=lv)
for(i in 1:lv) res[,i]<-fdoppler(f=f,vs=v[i])
op<-par(bg="lightgrey")
matplot(x=f,y=res,type="l",lty=1,las=1,col= spectro.colors(lv),
xlab="Source frequency (kHz)", ylab="Altered frequency (kHz)")
legend("topleft",legend=paste(as.character(v),"m/s"),
lty=1,col= spectro.colors(lv))
title(main="Doppler effect at different source speeds")
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ffilter")
### * ffilter

flush(stderr()); flush(stdout())

### Name: ffilter
### Title: Frequency filter
### Aliases: ffilter
### Keywords: ts

### ** Examples

a<-noisew(f=8000,d=1)
# low-pass
b<-ffilter(a,f=8000,to=1500)
spectro(b,f=8000,wl=512)
# high-pass
c<-ffilter(a,f=8000,from=2500)
spectro(c,f=8000,wl=512)
# band-pass
d<-ffilter(a,f=8000,from=1000,to=2000)
spectro(d,f=8000,wl=512)
# band-stop
e<-ffilter(a,f=8000,from=1500,to=2500,bandpass=FALSE)
spectro(e,f=8000,wl=512)
# custom
myfilter1<-rep(c(rep(0,32),rep(1,32)),4)
g<-fir(a,f=8000,custom=myfilter1)
spectro(g,f=8000)



cleanEx()
nameEx("field")
### * field

flush(stderr()); flush(stdout())

### Name: field
### Title: Near field and far field limits
### Aliases: field
### Keywords: ts

### ** Examples

# 1 kHz near field at 1 cm from the source
field(f=1000,d=0.01)
# playing with distance from source and sound frequency
op<-par(bg="lightgrey")
D<-seq(0.01,0.5,by=0.01); nD<-length(D)
F<-seq(100,1000,by=25); nF<-length(F)
a<-matrix(numeric(nD*nF),nrow=nD)
for(i in 1:nF) a[,i]<-field(f=F[i],d=D)$kd
matplot(x=D,y=a,type="l",lty=1,col= spectro.colors(nF),
  xlab="Distance from the source (m)", ylab="k*d")
title("Variation of the product k*d with distance and frequency")
text(x=c(0.4,0.15),y=c(0.02,1), c("Near Field","Far Field"),font=2)
legend(x=0.05,y=1.4,c("100 Hz","1000 Hz"),lty=1,
  col=c(spectro.colors(nF)[1],spectro.colors(nF)[nF]),bg="grey")
abline(h=0.1)
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("fir")
### * fir

flush(stderr()); flush(stdout())

### Name: fir
### Title: Finite Impulse Response filter
### Aliases: fir
### Keywords: ts

### ** Examples

a<-noisew(f=8000,d=1)
# low-pass
b<-fir(a,f=8000,to=1500)
spectro(b,f=8000)
# high-pass
c<-fir(a,f=8000,from=2500)
spectro(c,f=8000)
# band-pass
d<-fir(a,f=8000,from=1000,to=2000)
spectro(d,f=8000)
# band-stop
e<-fir(a,f=8000,from=1500,to=2500,bandpass=FALSE)
spectro(e,f=8000)
# custom filter manually generated
myfilter1<-rep(c(rep(0,32),rep(1,32)),4)
g<-fir(a,f=8000,custom=myfilter1)
spectro(g,f=8000)
# custom filter generated using spec()
data(tico)
myfilter2<-spec(tico,f=22050,at=0.7,wl=512,plot=FALSE)
b<-noisew(d=1,f=22050)
h<-fir(b,f=22050,custom=myfilter2)
spectro(h,f=22050)



cleanEx()
nameEx("fma")
### * fma

flush(stderr()); flush(stdout())

### Name: fma
### Title: Frequency modulation analysis
### Aliases: fma
### Keywords: dplot ts

### ** Examples

# a sound with a 1 Khz sinusoid FM
a<-synth(d=1, f=8000, cf=1500, fm=c(1000,1000,0))
fma(a,f=8000)



cleanEx()
nameEx("fpeaks")
### * fpeaks

flush(stderr()); flush(stdout())

### Name: fpeaks
### Title: Frequency peak detection
### Aliases: fpeaks
### Keywords: dplot ts

### ** Examples

data(tico)
spec <- meanspec(tico, f=22050, plot=FALSE)
specdB <- meanspec(tico, f=22050, dB="max0", plot=FALSE)
# all peaks
fpeaks(spec)
# 10 highest peaks
fpeaks(spec, nmax=10)
# highest peak (ie dominant frequency)
fpeaks(spec, nmax=1)
# peaks that are separated by more than 500 Hz
fpeaks(spec, freq=500)
# peaks with a left slope higher than 0.1
fpeaks(spec, amp=c(0.1,0))
# peaks with a right slope higher than 0.1
fpeaks(spec, amp=c(0,0.1))
# peaks with left and right slopes higher than 0.1
fpeaks(spec, amp=c(0.1,0.1))
# peaks above a 0.5 threshold
fpeaks(spec, threshold=0.5)
# peaks of a dB spectrum with peaks showing slopes higher than 3 dB
fpeaks(specdB, amp=c(3,3))
# comparing different parameter settings
meanspec(tico, f=22050)
col <- c("#ff000090","#0000ff75","#00ff00")
cex <- c(2,1.25,1.5)
pch <- c(19,17,4)
title(main="Peak detection \n (spectrum with values between 0 and 1)")
res1 <- fpeaks(spec, plot = FALSE)
res2 <- fpeaks(spec, amp=c(0.02,0.02), plot =FALSE)
res3 <- fpeaks(spec, amp=c(0.02,0.02), freq=200, plot = FALSE)
points(res1, pch=pch[1], col=col[1], cex=cex[1])
points(res2, pch=pch[2], col=col[2], cex=cex[2])
points(res3, pch=pch[3], col=col[3], cex=cex[3])
legend("topright", legend=c("all peaks","amp", "amp & freq"), pch=pch,
pt.cex=cex, col=col, bty="n")
# example with a cepstral spectrum
data(sheep)
res <- ceps(sheep,f=8000,at=0.4,wl=1024,plot=FALSE)
fpeaks(res, nmax=4, xlab="Quefrency (s)")



cleanEx()
nameEx("ftwindow")
### * ftwindow

flush(stderr()); flush(stdout())

### Name: ftwindow
### Title: Fourier transform windows
### Aliases: ftwindow
### Keywords: ts IO

### ** Examples

a<-ftwindow(512)
b<-ftwindow(512,wn="bartlett")
c<-ftwindow(512,wn="blackman")
d<-ftwindow(512,wn="flattop")
e<-ftwindow(512,wn="hanning")
f<-ftwindow(512,wn="rectangle")
all<-cbind(a,b,c,d,e,f)
matplot(all,type="l",col=1:6,lty=1:6)
legend(legend=c("hamming","bartlett","blackman","flattop","hanning","rectangle"),
x=380,y=0.95,col=1:6,lty=1:6,cex=0.75)



cleanEx()
nameEx("fund")
### * fund

flush(stderr()); flush(stdout())

### Name: fund
### Title: Fundamental frequency track
### Aliases: fund
### Keywords: dplot ts

### ** Examples

data(sheep)
fund(sheep,f=8000,fmax=300,type="l")
# with 50% overlap between successive sliding windows, time zoom and 
# amplitude filter (threshold)
fund(sheep,f=8000,fmax=300,type="b",ovlp=50,threshold=5,ylim=c(0,1),cex=0.5)
# overlaid on a spectrogram
spectro(sheep,f=8000,ovlp=75,zp=16,scale=FALSE,palette=rev.gray.colors.2)
par(new=TRUE)
fund(sheep,f=8000,fmax=300,type="p",pch=24,ann=FALSE,
  xaxs="i",yaxs="i",col="black",bg="red",threshold=6)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("hilbert")
### * hilbert

flush(stderr()); flush(stdout())

### Name: hilbert
### Title: Hilbert transform and analytic signal
### Aliases: hilbert
### Keywords: ts

### ** Examples

a<-synth(f=8000, d=1, cf=1000)
aa<-hilbert(a, f=8000)



cleanEx()
nameEx("ifreq")
### * ifreq

flush(stderr()); flush(stdout())

### Name: ifreq
### Title: Instantaneous frequency
### Aliases: ifreq
### Keywords: ts dplot

### ** Examples

# generate a sound with sine and linear frequency modulations
a<-synth(d=1, f=8000, cf=1500, fm=c(200,10,1000))
# plot on a single graphical device the instantaneous frequency and phase
op<-par(mfrow=c(2,1))
ifreq(a,f=8000,main="Instantaneous frequency")
ifreq(a,f=8000,phase=TRUE,main="Instantaneous phase")
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("itakura.dist")
### * itakura.dist

flush(stderr()); flush(stdout())

### Name: itakura.dist
### Title: Itakuro-Saito distance
### Aliases: itakura.dist
### Keywords: distribution ts

### ** Examples

# Comparison of two spectra
data(tico)
tico1 <- spec(tico, at=0.65, plot=FALSE)
tico2 <- spec(tico, at=1.1, plot=FALSE)
itakura.dist(tico1, tico2) 



cleanEx()
nameEx("kl.dist")
### * kl.dist

flush(stderr()); flush(stdout())

### Name: kl.dist
### Title: Kullback-Leibler distance
### Aliases: kl.dist
### Keywords: distribution ts

### ** Examples

# Comparison of two spectra
data(tico)
tico1 <- spec(tico, at=0.65, plot=FALSE)
tico2 <- spec(tico, at=1.1, plot=FALSE)
kl.dist(tico1, tico2)    # log2 (binary logarithm)
kl.dist(tico1, tico2, base=exp(1))  # ln (natural logarithm)



cleanEx()
nameEx("ks.dist")
### * ks.dist

flush(stderr()); flush(stdout())

### Name: ks.dist
### Title: Kolmogorov-Smirnov distance
### Aliases: ks.dist
### Keywords: distrubtion ts

### ** Examples

# Comparison of two spectra and plot of the cumulated spectra with the K-S distance
data(tico)
tico1 <- spec(tico, at=0.65, plot=FALSE)
tico2 <- spec(tico, at=1.1, plot=FALSE)
ks.dist(tico1, tico2, plot=TRUE)



cleanEx()
nameEx("lfs")
### * lfs

flush(stderr()); flush(stdout())

### Name: lfs
### Title: Linear Frequency Shift
### Aliases: lfs
### Keywords: ts

### ** Examples

data(orni)
a<-lfs(orni,f=22050,shift=1000)
spectro(a,f=22050)
# to be compared with the original signal
spectro(orni,f=22050)



cleanEx()
nameEx("listen")
### * listen

flush(stderr()); flush(stdout())

### Name: listen
### Title: Play a sound wave
### Aliases: listen
### Keywords: ts

### ** Examples

## NOT RUN
# data(tico)
# listen(tico,f=22050)
# listen(tico,f=22050,from=0.5,to=1.5)
# listen(noise(d=1,f=8000,Wave=TRUE))
## change f to play the sound a different speed
# data(sheep)
## normal
# listen(sheep,f=8000)
## two times faster
# listen(sheep,f=8000*2)
## two times slower
# listen(sheep,f=8000/2)



cleanEx()
nameEx("localpeaks")
### * localpeaks

flush(stderr()); flush(stdout())

### Name: localpeaks
### Title: Local maximum frequency peak detection
### Aliases: localpeaks
### Keywords: dplot ts

### ** Examples

data(sheep)
spec <- meanspec(sheep, f=8000)
# a specific number of bands with all the same size
localpeaks(spec, bands=5)
# bands directly specified  with a regular sequence
localpeaks(spec, bands=seq(0,8/2,by=0.5))
# bands directly specified  with an irregular sequence
localpeaks(spec, bands=c(0,0.5,1,1.5,3,4))
# Amaj octave bands, note that there is no peak detection
# in the higher part of the spectrum as sequence stops at 3520 Hz
localpeaks(spec, bands=octaves(440, below=3, above=3)/1000) 



cleanEx()
nameEx("logspec.dist")
### * logspec.dist

flush(stderr()); flush(stdout())

### Name: logspec.dist
### Title: Log-spectral distance
### Aliases: logspec.dist
### Keywords: distribution ts

### ** Examples

# Comparison of two spectra
data(tico)
tico1 <- spec(tico, at=0.65, plot=FALSE)
tico2 <- spec(tico, at=1.1, plot=FALSE)
logspec.dist(tico1, tico2)



cleanEx()
nameEx("meandB")
### * meandB

flush(stderr()); flush(stdout())

### Name: meandB
### Title: mean of dB values
### Aliases: meandB
### Keywords: math

### ** Examples

meandB(c(89,90,95))



cleanEx()
nameEx("meanspec")
### * meanspec

flush(stderr()); flush(stdout())

### Name: meanspec
### Title: Mean frequency spectrum of a time wave
### Aliases: meanspec
### Keywords: dplot ts

### ** Examples

data(orni)
# compute the mean spectrum of the whole time wave
meanspec(orni,f=22050)
# compute the mean spectrum of a time wave section (from 0.32 s to 0.39 s)
meanspec(orni,f=22050,from=0.32,to=0.39)
# different window lengths
op<-par(mfrow=c(3,1))
meanspec(orni,f=22050,wl=256)
title("wl=256")
meanspec(orni,f=22050,wl=1024)
title("wl=1024")
meanspec(orni,f=22050,wl=4096)
title("wl=4096")
par(op)
# different overlap values (almost no effects here...)
op<-par(mfrow=c(3,1))
meanspec(orni,f=22050)
title("ovlp=0")
meanspec(orni,f=22050,ovlp=50)
title("ovlp=50")
meanspec(orni,f=22050,ovlp=95)
title("ovlp=95")
par(op)
# use of flim to zoom in
op<-par(mfrow=c(2,1))
meanspec(orni,f=22050)
title("zoom in")
meanspec(orni,f=22050,wl=512,flim=c(4,6))
par(op)
# comparaison of spectrum and mean spectrum
op<-par(mfrow=c(2,1))
spec(orni,f=22050)
title("spec()")
meanspec(orni,f=22050)
title("meanspec()")
par(op)
# log scale on frequency axis
meanspec(orni, f=22050, log="x")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("mel")
### * mel

flush(stderr()); flush(stdout())

### Name: mel
### Title: Hertz / Mel conversion
### Aliases: mel
### Keywords: math

### ** Examples

x<-seq(0,10000,by=50)
y<-mel(x)
plot(x,y,type="l",xlab = "f (hertz)", ylab = "f (mel)",
  main = "Mel scale", col="red")



cleanEx()
nameEx("micsens")
### * micsens

flush(stderr()); flush(stdout())

### Name: micsens
### Title: Microphone sensitivity and conversion
### Aliases: micsens
### Keywords: math

### ** Examples

# conversion of a sensitivity of 2 mV/Pa
micsens(2)
# conversion of a sensitivity of -54 dB re 1V/Pa
micsens(-54,inverse=TRUE)



cleanEx()
nameEx("moredB")
### * moredB

flush(stderr()); flush(stdout())

### Name: moredB
### Title: Addition of dB values
### Aliases: moredB
### Keywords: math

### ** Examples

# two sources of 60 dB give an intensity or pressure level of 63 dB
moredB(c(60,60))
# addition of three sources
moredB(c(89,90,95))



cleanEx()
nameEx("mutew")
### * mutew

flush(stderr()); flush(stdout())

### Name: mutew
### Title: Replace time wave data by 0 values
### Aliases: mutew
### Keywords: dplot ts

### ** Examples

data(tico)
mutew(tico,f=22050,from=0.5,to=0.9)



cleanEx()
nameEx("noisew")
### * noisew

flush(stderr()); flush(stdout())

### Name: noisew
### Title: Generate noise
### Aliases: noisew
### Keywords: datagen ts

### ** Examples

# add noise to a synthetic signal
a<-noisew(d=1,f=8000)
b<-synth(f=8000,d=1,cf=2000,plot=FALSE)
c<-a+b
spectro(c,f=8000)



cleanEx()
nameEx("notefreq")
### * notefreq

flush(stderr()); flush(stdout())

### Name: notefreq
### Title: Frequency of a muscical note
### Aliases: notefreq
### Keywords: maths

### ** Examples

# Some notes frequency (use apply-like functions when dealing with character strings)
sapply(c("C", "A", "Gb"), notefreq)

# C major scale plot
n <- 1:12
freq <- notefreq(n)
names <- c("C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B") 
plot(n, freq, pch=19, cex=1.5,
     xlab = "Note name",
     ylab = "Frequency (Hz)",
     xaxt="n", las=1, main="Third octave")
axis(side=1, at=n, labels=names)
abline(h=freq, col="lightgrey")

# C major scale sound
f <- 2000 # sampling rate
s <- NULL
for (i in 1:length(freq))
  {
    tmp <- synth(d=0.5, f=f, cf=freq[i])
    s <- pastew(s, tmp, at="start", f)
  }
spectro(s, f, ovlp=75)



cleanEx()
nameEx("octaves")
### * octaves

flush(stderr()); flush(stdout())

### Name: octaves
### Title: Octave values
### Aliases: octaves
### Keywords: maths

### ** Examples

names <- c("C","D","E","F","G","A","B")
values <- c(261.63, 293.66, 329.64, 349.23, 392, 440, 493.88)
res <- sapply(values, FUN=octaves)/1000
op <- par(las=1,mfrow=c(2,1))
par(mar=c(0,4,1,1))
matplot(x=1:7, y=res, t="o", pch=names, xlab="",
    ylab="Frequency (kHz) [linear scale]", col=rainbow(7), xaxt="n")
par(mar=c(4.5,4,0,1))
matplot(x=1:7, y=res, t="o", pch=names, xlab="Octave",
    ylab="Frequency (kHz) [log scale]", col=rainbow(7), ylog=TRUE, log="y")
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("orni")
### * orni

flush(stderr()); flush(stdout())

### Name: orni
### Title: Song of the cicada Cicada orni
### Aliases: orni
### Keywords: datasets

### ** Examples

data(orni)
oscillo(orni,f=22050)



cleanEx()
nameEx("oscillo")
### * oscillo

flush(stderr()); flush(stdout())

### Name: oscillo
### Title: Show a time wave as an oscillogram
### Aliases: oscillo
### Keywords: dplot ts

### ** Examples

data(tico)
# a simple oscillogram of a bird song
oscillo(tico,f=22050)
# zoom in
op<-par(mfrow=c(4,1),mar=c(4.5,4,2,2))
oscillo(tico,22050,cexlab=0.75)
oscillo(tico,22050,from=0.5,to=0.9,cexlab=0.75)
oscillo(tico,22050,from=0.65,to=0.75,cexlab=0.75)
oscillo(tico,22050,from=0.68,to=0.70,cexlab=0.75)
par(op)
# the same divided in four lines
oscillo(tico,f=22050,k=4,j=1)
# the same divided in different numbers of lines and columns
oscillo(tico,f=22050,k=4,j=4)
oscillo(tico,f=22050,k=2,j=2,byrow=TRUE)
oscillo(tico,f=22050,k=2,j=2,byrow=FALSE)
# overplot of oscillographic and envelope representations
oscillo(tico,f=22050)
par(new=TRUE)
env(tico,f=22050,colwave=2)
# full colour modifications in a two-frame oscillogram
op<-par(bg="grey")
oscillo(tico,f=22050,k=4,j=1,title=TRUE,colwave="black",
    coltitle="yellow",collab="red",colline="white",
    colaxis="blue",coly0="grey50")
par(op)
# change the title
data(orni)
oscillo(orni,f=22050,title="The song of a famous cicada")
# move along the signal using scroll
require(rpanel)
oscillo(tico,f=22050,scroll=8)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("oscilloST")
### * oscilloST

flush(stderr()); flush(stdout())

### Name: oscilloST
### Title: Show a stereo time wave as oscillograms
### Aliases: oscilloST
### Keywords: dplot ts

### ** Examples

a<-synth(f=8000,d=1,cf=2000,am=c(50,10),plot=FALSE)
b<-synth(f=8000,d=1,cf=1000,fm=c(0,0,2000),plot=FALSE)
oscilloST(a,b,f=8000)



cleanEx()
nameEx("pastew")
### * pastew

flush(stderr()); flush(stdout())

### Name: pastew
### Title: Paste a time wave to another one
### Aliases: pastew
### Keywords: dplot ts

### ** Examples

data(tico)
# double a data set describing a bird song
a<-pastew(tico,tico,f=22050)
oscillo(a,f=22050)
# a direct way to see what has been pasted
pastew(tico,tico,f=22050,plot=TRUE)
# cut a section and then paste it at the beginning
a<-cutw(tico, f=22050, from=0.5, to=0.9)
pastew(a,tico,f=22050,at="start",plot=TRUE)
# or paste it at a specific location
pastew(a,tico,f=22050,at=1.4,plot=TRUE)
# setting the argument 'join' to TRUE might be useful
# to smooth pasting when some phase problem occur
# generate two sine waves
a <- synth(cf=50, f=400, d=0.1)
b <- synth(cf=100, f=400, d=0.1)
# paste it with 'join' turned to FALSE
# there is a click at the junction between the two waves
pastew(a, b, f=400, plot=TRUE)
# that can be removed by setting 'join' to TRUE
pastew(a, b, f=400, join=TRUE, plot=TRUE)



cleanEx()
nameEx("peewit")
### * peewit

flush(stderr()); flush(stdout())

### Name: peewit
### Title: Song of the bird Vanellus vanellus
### Aliases: peewit
### Keywords: datasets

### ** Examples

data(peewit)
oscillo(peewit,f=22050)



cleanEx()
nameEx("pellucens")
### * pellucens

flush(stderr()); flush(stdout())

### Name: pellucens
### Title: Calling song of the tree cricket Oecanthus pellucens
### Aliases: pellucens
### Keywords: datasets

### ** Examples

data(pellucens)
oscillo(pellucens,f=11025)



cleanEx()
nameEx("phaseplot")
### * phaseplot

flush(stderr()); flush(stdout())

### Name: phaseplot
### Title: Phase-phase 2D or 3D plot of a time wave
### Aliases: phaseplot
### Keywords: dplot ts

### ** Examples

require(rgl)
data(tico)
phaseplot(tico)



cleanEx()
nameEx("pulse")
### * pulse

flush(stderr()); flush(stdout())

### Name: pulse
### Title: Generate rectangle pulse
### Aliases: pulse
### Keywords: datagen ts

### ** Examples

pulse(dbefore=0.5,dpulse=0.1,dafter=0.3,f=8000,plot=TRUE)



cleanEx()
nameEx("repw")
### * repw

flush(stderr()); flush(stdout())

### Name: repw
### Title: Repeat a time wave
### Aliases: repw
### Keywords: dplot ts

### ** Examples

data(tico)
repw(tico,f=22050,plot=TRUE)
# use 'join' for smooth pasting
par(mfrow=c(2,1))
a <- synth(cf=50, f=400, d=0.1)
repw(a, f=400, plot=TRUE)
title(main="join is FALSE")
points(x=0.1, y=0, cex=2, col=2)
repw(a, f=400, join=TRUE, plot=TRUE)
title(main="join is TRUE")
points(x=0.1, y=0, cex=2, col=2)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("resamp")
### * resamp

flush(stderr()); flush(stdout())

### Name: resamp
### Title: Resample a time wave
### Aliases: resamp
### Keywords: dplot ts

### ** Examples

data(peewit)
# downsampling
a<-resamp(peewit,f=22050,g=11025)
# oversampling
b<-resamp(peewit,f=22050,g=44100)



cleanEx()
nameEx("revw")
### * revw

flush(stderr()); flush(stdout())

### Name: revw
### Title: Time reverse of a time wave
### Aliases: revw
### Keywords: dplot ts

### ** Examples

data(tico)
# simple reverse
revw(tico,f=22050,plot=TRUE)
# envelope reverse only
revw(tico,f=22050,ifreq=FALSE, plot=TRUE)
# instantaneous frequency reverse only
revw(tico,f=22050,env=FALSE, plot=TRUE)



cleanEx()
nameEx("rmam")
### * rmam

flush(stderr()); flush(stdout())

### Name: rmam
### Title: Remove the amplitude modulations of a time wave
### Aliases: rmam
### Keywords: ts

### ** Examples

# generate a new sound with amplitude modulation
a<-synth(f=8000, d=1, cf=1500, am=c(50,10))
# remove the amplitude modulation and plot the result
rmam(a,f=8000,plot=TRUE)



cleanEx()
nameEx("rmnoise")
### * rmnoise

flush(stderr()); flush(stdout())

### Name: rmnoise
### Title: Remove noise
### Aliases: rmnoise
### Keywords: ts

### ** Examples

# synthesis of a 440 Hz sound with background noise
n <- noisew(d=1,f=8000)
s <- synth(d=1,f=8000,cf=440)
ns <- n+s
# remove noise (but low frequency content still there)
a <- rmnoise(ns,f=8000)



cleanEx()
nameEx("rmoffset")
### * rmoffset

flush(stderr()); flush(stdout())

### Name: rmoffset
### Title: Remove the offset of a time wave
### Aliases: rmoffset
### Keywords: dplot ts

### ** Examples

data(tico)
# artifically generates an offset
tico2<-tico+0.1
# see the wave with an offset
oscillo(tico2,f=22050)
# remove the offset
rmoffset(tico2,f=22050,plot=TRUE)



cleanEx()
nameEx("rms")
### * rms

flush(stderr()); flush(stdout())

### Name: rms
### Title: Root Mean Square
### Aliases: rms
### Keywords: ts

### ** Examples

# simple rms
rms(1:10)
# rms of a normalized envelope
data(sheep)
env <- env(sheep, f=8000)
rms(env)



cleanEx()
nameEx("roughness")
### * roughness

flush(stderr()); flush(stdout())

### Name: roughness
### Title: Roughness or total curvature
### Aliases: roughness
### Keywords: ts

### ** Examples

data(tico)
spec <- meanspec(tico, plot=FALSE)[,2]
roughness(spec) 



cleanEx()
nameEx("rugo")
### * rugo

flush(stderr()); flush(stdout())

### Name: rugo
### Title: Rugosity of a time wave
### Aliases: rugo
### Keywords: ts

### ** Examples

data(tico) ; tico <-tico@left
# rugosity of the original recording
rugo(tico/max(tico))
# add artificially some background noise
noise <- noisew(d=length(tico)/22050, f=22050)
ticon1 <- tico/max(tico) + 0.5*noise
# new rugosity (higher)
rugo(ticon1/max(ticon1))



cleanEx()
nameEx("savewav")
### * savewav

flush(stderr()); flush(stdout())

### Name: savewav
### Title: Save a .wav file
### Aliases: savewav
### Keywords: IO

### ** Examples

require(tuneR)
a<-synth(f=8000,d=2,cf=2000,plot=FALSE)
# the name of the file is automatically the name of the object
# here: "a.wav"
savewav(a,f=22050)
unlink("a.wav")
# if you wish to change the name, use 'file' argument
savewav(a,f=22050,file="b.wav")
unlink("b.wav")



cleanEx()
nameEx("seedata")
### * seedata

flush(stderr()); flush(stdout())

### Name: seedata
### Title: A quick look at quantitative data
### Aliases: seedata
### Keywords: dplot

### ** Examples
seedata(rnorm(1000))


cleanEx()
nameEx("setenv")
### * setenv

flush(stderr()); flush(stdout())

### Name: setenv
### Title: Set the amplitude envelope of a time wave to another one
### Aliases: setenv
### Keywords: datagen ts

### ** Examples

data(tico)
a<-synth(d=1,f=22050,cf=1000)
# apply 'tico' ammplitude envelope to 'a' that has a square amplitude envelope
setenv(a,tico,f=22050,plot=TRUE)
# the same but with smoothing the envelope
setenv(a,tico,f=22050,ksmooth=kernel("daniell",50),plot=TRUE)



cleanEx()
nameEx("sfm")
### * sfm

flush(stderr()); flush(stdout())

### Name: sfm
### Title: Spectral Flatness Measure
### Aliases: sfm
### Keywords: ts

### ** Examples

a<-synth(f=8000,d=1,cf=2000,plot=FALSE)
speca<-spec(a,f=8000,at=0.5,plot=FALSE)
sfm(speca)
# [1] 0
b<-noisew(d=1,f=8000)
specb<-spec(b,f=8000,at=0.5,plot=FALSE)
sfm(specb)
# [1] 0.8233202



cleanEx()
nameEx("sh")
### * sh

flush(stderr()); flush(stdout())

### Name: sh
### Title: Shannon and Renyi spectral entropy
### Aliases: sh
### Keywords: ts

### ** Examples

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



cleanEx()
nameEx("sheep")
### * sheep

flush(stderr()); flush(stdout())

### Name: sheep
### Title: Sheep bleat
### Aliases: sheep
### Keywords: datasets

### ** Examples

data(sheep)
oscillo(sheep,f=8000)



cleanEx()
nameEx("simspec")
### * simspec

flush(stderr()); flush(stdout())

### Name: simspec
### Title: Similarity between two frequency spectra
### Aliases: simspec
### Keywords: dplot ts

### ** Examples

a<-noisew(f=8000,d=1)
b<-synth(f=8000,d=1,cf=2000)
c<-synth(f=8000,d=1,cf=1000)
d<-noisew(f=8000,d=1)
speca<-spec(a,f=8000,at=0.5,plot=FALSE)
specb<-spec(b,f=8000,at=0.5,plot=FALSE)
specc<-spec(c,f=8000,at=0.5,plot=FALSE)
specd<-spec(d,f=8000,at=0.5,plot=FALSE)
simspec(speca,speca)
simspec(speca,specb)
simspec(speca,specc,plot=TRUE)
simspec(specb,specc,plot=TRUE)
#[1] 12.05652
simspec(speca,specd,plot=TRUE)



cleanEx()
nameEx("smoothw")
### * smoothw

flush(stderr()); flush(stdout())

### Name: smoothw
### Title: A function to tentavily smooth a time wave
### Aliases: smoothw
### Keywords: ts

### ** Examples

# An example to show that smoothw() may change
# the frequency content of your sound
data(orni)
orni2 <- smoothw(orni, wl=2, out="Wave")
orni10 <- smoothw(orni, wl=10, out="Wave")
orni50 <- smoothw(orni, wl=50, out="Wave")
orni100 <- smoothw(orni, wl=100, out="Wave")
meanspec(orni)
lines(meanspec(orni2, plot=FALSE), col=2)
lines(meanspec(orni10, plot=FALSE), col=3)
lines(meanspec(orni50, plot=FALSE), col=4)
lines(meanspec(orni100, plot=FALSE), col=5)
legend("topright", col=1:5, lty=1, legend=c("original","wl=2","wl=10","wl=50","wl=100"))



cleanEx()
nameEx("spec")
### * spec

flush(stderr()); flush(stdout())

### Name: spec
### Title: Frequency spectrum of a time wave
### Aliases: spec
### Keywords: dplot ts

### ** Examples

data(tico)
# spectrum of the whole signal, in absolute or dB amplitude,
# horizontaly or vertically
op<-par(mfrow=c(2,2))
spec(tico,f=22050)
spec(tico,f=22050,col="red",plot=2)
spec(tico,f=22050,dB="max0",col="blue")
spec(tico,f=22050,dB="max0",col="green",plot=2)
par(op)
# an indirect way to compare spectra 
a<-spec(tico,f=22050,wl=512,at=0.2,plot=FALSE)
b<-spec(tico,f=22050,wl=512,at=0.7,plot=FALSE)
c<-spec(tico,f=22050,wl=512,at=1.1,plot=FALSE)
d<-spec(tico,f=22050,wl=512,at=1.6,plot=FALSE)
all<-cbind(a[,2],b[,2],c[,2],d[,2])
matplot(x=a[,1],y=all,yaxt="n",
    xlab="Frequency (kHz)",ylab="Amplitude",xaxs="i",type="l")
legend(8,0.8,c("Note A","Note B", "Note C", "Note D"),bty="o",
    lty=c(1:4),col=c(1:4))
# spectrum from a particular position to another one
op<-par(mfrow=c(2,1))
oscillo(tico,f=22050)
abline(v=c(0.5,0.9),col="red",lty=2)
spec(tico,f=22050,wl=512,from=0.5,to=0.9,col="red")
title("Spectrum of the note B")
par(op)
# spectrum and spectrogram
data(orni)
orni1<-cutw(orni,f=22050,from=0.32,to=0.39)
layout(matrix(c(1,2),nc=2),widths=c(3,1))
par(mar=c(5,4,3,0.5))
spectro(orni1,f=22050,wl=128,zp=8,ovlp=85,scale=FALSE)
par(mar=c(5,1,3,0.5))
spec(orni1,f=22050,col="red",plot=2,flab="",yaxt="n")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("specprop")
### * specprop

flush(stderr()); flush(stdout())

### Name: specprop
### Title: Spectral properties
### Aliases: specprop
### Keywords: dplot ts

### ** Examples

data(orni)
a<-meanspec(orni,f=22050,plot=FALSE)
specprop(a,f=22050)
# to get a single measure of the list
specprop(a,f=22050)$mode
# to get the results structured
specprop(a,f=22050,str=TRUE)
# to limit the analysis between 4 and 6 kHz
specprop(a,f=22050,flim=c(4,6),str=TRUE)
# plots
specprop(a,f=22050,plot=1)
specprop(a,f=22050,plot=2)



cleanEx()
nameEx("spectro")
### * spectro

flush(stderr()); flush(stdout())

### Name: spectro
### Title: 2D-spectrogram of a time wave
### Aliases: spectro
### Keywords: dplot ts

### ** Examples

data(tico)
data(pellucens)
# simple plots
spectro(tico,f=22050)
spectro(tico,f=22050,osc=TRUE)
spectro(tico,f=22050,scale=FALSE)
spectro(tico,f=22050,osc=TRUE,scale=FALSE)
# change the dB scale by setting a different dB reference value (20 microPa)
spectro(tico,f=22050, dBref=2*10e-5)
# manipulating wl
op<-par(mfrow=c(2,2))
spectro(tico,f=22050,wl=256,scale=FALSE)
title("wl = 256")
spectro(tico,f=22050,wl=512,scale=FALSE)
title("wl = 512")
spectro(tico,f=22050,wl=1024,scale=FALSE)
title("wl = 1024")
spectro(tico,f=22050,wl=4096,scale=FALSE)
title("wl = 4096")
par(op)
# vertical zoom using flim
spectro(tico,f=22050, ylim=c(2,6))
spectro(tico,f=22050, ylimd=c(2,6))
# a full plot
pellu2<-cutw(pellucens,f=22050,from=1,plot=FALSE)
spectro(pellu2,f=22050,ovlp=85,zp=16,osc=TRUE,
    cont=TRUE,contlevels=seq(-30,0,20),colcont="red",
    lwd=1.5,lty=2,palette=rev.terrain.colors)
# black and white spectrogram 
spectro(pellu2,f=22050,ovlp=85,zp=16,
    palette=rev.gray.colors.1)
# colour modifications
data(sheep)
spectro(sheep,f=8000,palette=temp.colors,collevels=seq(-115,0,1))
spectro(pellu2,f=22050,ovlp=85,zp=16,
palette=rev.cm.colors,osc=TRUE,colwave="orchid1") 
spectro(pellu2,f=22050,ovlp=85,zp=16,osc=TRUE,palette=rev.heat.colors,
colbg="black",colgrid="white", colwave="white",colaxis="white",collab="white")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("spectro3D")
### * spectro3D

flush(stderr()); flush(stdout())

### Name: spectro3D
### Title: 3D-spectrogram of a time wave
### Aliases: spectro3D
### Keywords: dplot ts

### ** Examples

require(rgl)
data(tico)
spectro3D(tico,f=22050,wl=512,ovlp=75,zp=16,maga=4,palette=rev.terrain.colors)



cleanEx()
nameEx("symba")
### * symba

flush(stderr()); flush(stdout())

### Name: symba
### Title: Symbol analysis of a numeric (time) series
### Aliases: symba
### Keywords: ts

### ** Examples

# analysis of a frequency spectrum
data(tico)
spec1<-spec(tico,f=22050,at=0.2,plot=FALSE)
symba(spec1[,2],plot=TRUE)
# it might be better to round the values
symba(round(spec1[,2],2),plot=TRUE)
# in that case the symbol entropy is almost similar to the spectral entropy
symba(round(spec1[,2],2),entrop="rel")$h1
sh(spec1)
# to compare two frequency spectra
spec2<-spec(tico,f=22050,wl=512,at=1.1,plot=FALSE)
symba(round(spec1[,2],2),round(spec2[,2],2),plot=TRUE)



cleanEx()
nameEx("synth")
### * synth

flush(stderr()); flush(stdout())

### Name: synth
### Title: Synthesis of time wave
### Aliases: synth
### Keywords: datagen ts

### ** Examples

# pure tone
synth(f=22050,d=1,cf=4000,plot=TRUE)
# pure tone with sinusoid-like overall shape
synth(f=22050,d=1,cf=4000,shape="sine",plot=TRUE,osc=TRUE)
# pure tones with am
synth(f=22050,d=1,cf=4000,am=c(50,10),plot=TRUE,osc=TRUE)
# pure tone with +2000 Hz linear fm 
synth(f=22050,d=1,cf=4000,fm=c(0,0,2000),plot=TRUE)
# pure tone with sinusoidal fm
# (maximum excursion of 1000 Hz, frequency of 10 Hz)
synth(f=22050,d=1,cf=4000,fm=c(1000,10,0),plot=TRUE,wl=256,ovlp=75)
# pure tone with sinusoidal am
# (maximum excursion of 1000 Hz, frequency of 10 Hz)
# and linear fm (maximum excursion of 1000 Hz)
synth(f=22050,d=1,cf=4000,fm=c(1000,10,1000),plot=TRUE,wl=256,ovlp=75)
# the same with am
synth(f=22050,d=1,cf=4000,am=c(50,10),
    fm=c(1000,10,1000),plot=TRUE,wl=256,ovlp=75,osc=TRUE)
# the same with am and a triangular overall shape 
synth(f=22050,d=1,cf=4000,shape="tria",am=c(50,10),
    fm=c(1000,10,1000),plot=TRUE,wl=256,ovlp=75,osc=TRUE)   
# more complex sound
F1<-synth(f=22050,cf=2000,d=1,fm=c(500,5,0))
F2<-synth(f=22050,a=0.8,cf=4000,d=1,fm=c(500,5,0))
F3<-synth(f=22050,a=0.6,cf=6000,d=1,fm=c(500,5,0))
F4<-synth(f=22050,a=0.4,cf=8000,d=1,fm=c(500,5,0))
final1<-F1+F2+F3+F4
spectro(final1,f=22050,wl=512,ovlp=75,scale=FALSE)



cleanEx()
nameEx("th")
### * th

flush(stderr()); flush(stdout())

### Name: th
### Title: Temporal entropy
### Aliases: th
### Keywords: ts

### ** Examples

# Temporal entropy of a cicada song
data(orni)
envorni<-env(orni,f=22050,plot=FALSE)
th(envorni)
# Smoothing the envelope might slightly change the result.
envorniS<-env(orni,f=22050,smooth=c(50,0),plot=FALSE)
th(envorniS)
# If we mute a part of the cicada song, the temporal entropy decreases
orni2<-mutew(orni,f=22050,from=0.3,to=0.55,plot=FALSE)
envorni2<-env(orni2,f=22050,plot=FALSE)
th(envorni2)
# The temporal entropy of noise tends towards 1
a<-noisew(d=1,f=8000)
enva<-env(a,f=8000,plot=FALSE)
th(enva)
# But be aware that the temporal entropy
# of a sustained sound also tends towards 1
b<-synth(f=8000,d=1,cf=2000,plot=FALSE)
envb<-env(b,f=8000,plot=FALSE)
th(envb)



cleanEx()
nameEx("tico")
### * tico

flush(stderr()); flush(stdout())

### Name: tico
### Title: Song of the bird Zonotrichia capensis
### Aliases: tico
### Keywords: datasets

### ** Examples

data(tico)
oscillo(tico,f=22050)



cleanEx()
nameEx("timer")
### * timer

flush(stderr()); flush(stdout())

### Name: timer
### Title: Time measurements of a time wave
### Aliases: timer
### Keywords: dplot ts

### ** Examples

data(tico)
timer(tico,f=22050,threshold=5,msmooth=c(50,0))
# to compare with an oscillographic representation
data(orni)
op<-par(mfrow=c(2,1))
timer(orni,f=22050,threshold=5,msmooth=c(40,0),tck=0.05,
        bty="l",colval="blue")
title(main="A cicada song made of five echemes",col="blue")
oscillo(orni,f=22050,k=1,j=1)
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("wasp")
### * wasp

flush(stderr()); flush(stdout())

### Name: wasp
### Title: WAve length and SPeed of sound
### Aliases: wasp
### Keywords: math

### ** Examples

# wavelength (m) of a 2000 Hz air-borne sound at 20 degrees Celsius
wasp(f=2000)$l
# [1] 0.1717

# sound speed in sea at 0 and -500 m for a respective temperature of 22 degrees Celcius and 11 degrees Celsius
wasp(f=1000,s=30,d=c(0,500),t=c(22,11),medium="sea")$c
# [1] 1521.246 1495.414

# wavelength (m) of a 1000 Hz sound in a medium unspecified where c = 1497 m/s
wasp(f=1000,c=1497)$l
# [1] 1.497

# variation of wavelength according to frequency and air temperature
op<-par(bg="lightgrey")
a<-seq(1000,20000,by=100) ; na<-length(a)
b<-seq(-20,40,by=10) ; nb<-length(b)
res<-matrix(numeric(na*nb),nrow=na)
for(i in 1:nb) res[,i]<-wasp(a,t=b[i])$l
matplot(x=a,y=res,type="l",lty=1,col= spectro.colors(nb),
  xlab="Frequency (Hz)",ylab="Wavelength (m)")
title("Wavelength of air-borne sound at different temperatures (deg. C)")
legend(x=15000,y=0.3,c("-20","-10","0","10","20","30","40"),
  lty=1,col= spectro.colors(nb),bg="grey")
par(op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("wav2flac")
### * wav2flac

flush(stderr()); flush(stdout())

### Name: wav2flac
### Title: wav-flac file conversion
### Aliases: wav2flac
### Keywords: IO

### ** Examples

if(nzchar(Sys.which("flac"))) # check that FLAC is installed on your system
{
# synthesis of a 1kHz sound
a<-synth(d=10,f=8000,cf=1000)
# save it as a .wav file in the default working directory
savewav(a,f=8000)
# compress it to FLAC format and overwrite on the file a.wav
wav2flac("a.wav", overwrite=TRUE)
# back to .wav format
wav2flac("a.flac", reverse=TRUE)
# remove the files
unlink(c("a.wav","a.flac"))
}



cleanEx()
nameEx("wf")
### * wf

flush(stderr()); flush(stdout())

### Name: wf
### Title: Waterfall display
### Aliases: wf
### Keywords: dplot ts

### ** Examples

data(tico)
wf(tico,f=22050)
# changing the display parameters
jet.colors <- colorRampPalette(c("blue", "green"))
wf(tico,f=22050, hoff=0, voff=2, col=jet.colors, border = NA)
# matrix input instead of a time wave and transparent lines display
m <- numeric()
for(i in seq(-pi,pi,len=40)) {m <- cbind(m,10*(sin(seq(0,2*pi,len=100)+i)))}
wf(x=m, lines=TRUE, col="#0000FF50",xlab="Time", ylab="Amplitude",
main="waterfall display")



cleanEx()
nameEx("zapsilw")
### * zapsilw

flush(stderr()); flush(stdout())

### Name: zapsilw
### Title: Zap silence periods of a time wave
### Aliases: zapsilw
### Keywords: ts

### ** Examples

data(orni)
zapsilw(orni,f=22050,colwave="red")
# setting the threshold value
zapsilw(orni,f=22050,threshold=1)



cleanEx()
nameEx("zc")
### * zc

flush(stderr()); flush(stdout())

### Name: zc
### Title: Instantaneous frequency of a time wave by zero-crossing
### Aliases: zc
### Keywords: dplot ts

### ** Examples

data(pellucens)
pellu1<-cutw(pellucens,f=22050,from=0,to=1,plot=FALSE)
# without interpolation
zc(pellu1,f=22050,threshold=5,pch=20)
# with interpolation
zc(pellu1,f=22050,threshold=5,interpol=20,pch=20)
# a way to plot with a line and to filter low frequencies
pellu2<-zc(pellu1,f=22050,threshold=5,interpol=20,plot=FALSE)
pellu3<-na.omit(pellu2[,2])
pellu4<-pellu3[pellu3>3]
plot(x=seq(0,nrow(pellu1)/22050,length.out=length(pellu4)),
    y=pellu4,type="l",xlab="Time(s)",ylab="Frequency(kHz)")



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
