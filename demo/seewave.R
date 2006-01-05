require(seewave); data(tico); data(orni); data(pellucens); data(alauda)

op1<-par(ask=TRUE)

# different oscillograms of a tropical sparrow song
oscillo(tico,f=22050)
oscillo(tico,f=22050,k=2,j=2,byrow=TRUE)
oscillo(tico,f=22050,k=4,j=1,title=TRUE,colwave="green",colbg="grey",
    coltitle="yellow",collab="red",colline="white",
    colaxis="blue",coly0="grey50")

# overplot of oscillographic and envelope representation
env<-oscillo(tico,f=22050,env=TRUE,plot=FALSE,smooth=40)
ticonorm<-tico/max(tico)
envnorm<-env/max(env)
oscillo(ticonorm,f=22050)
par(new=TRUE)
plot(envnorm,type="l",col="red",xaxs="i",yaxs="i",ann=FALSE,xaxt="n",yaxt="n",
    ylim=range(ticonorm),bty="l",lwd=2)
legend(x=4, y=1,"smoothed envelope", col="red",lty=1,lwd=2,bty="n",cex=0.75)

# temporal automatic measurements
timer(orni,f=22050,threshold=5,smooth=40,
        bty="l",xaxs="i",colval="blue")
title(main="Timer() for automatic time measurements",col="blue")

# comparaison of a full spectrum and a mean spectrum of a cicada song
op<-par(mfrow=c(2,1))
spec(orni,f=22050,type="l")
title("spec()")
meanspec(orni,f=22050,wl=512,type="l")
title("meanspec()")
par(op)

# basic 2D spectrogram of a bird song
spectro(tico,f=22050,wl=512,ovlp=50,zp=16,collevels=seq(-40,0,0.5))

# spectrogram and dominant frequency overlaid of a bird song
spectro(alauda,f=22050,wl=512,ovlp=75,zp=8,palette=rev.gray.colors,scale=FALSE)
par(new=TRUE)
dfreq(alauda,f=22050,wl=512,ovlp=50,threshold=6,type="l",col="red",lwd=2,
    ann=FALSE,xaxs="i",yaxs="i")

# 2D spectrogram of a cricket song with colour modifications
op<-par(bg="black",col="white")
pellu2<-cutw(pellucens,f=22050,from=1,to=nrow(pellucens)/22050,plot=FALSE)
spectro(pellu2,f=22050,wl=512,ovlp=85,collevels=seq(-25,0,1),osc=TRUE,palette=rev.heat.colors,
colgrid="white", colwave="white",colaxis="white",collab="white",
colline="white")
par(op)

# sound synthesis
F1<-synth(f=22050,am=c(50,10),cf=2000,d=1,fm=c(500,5,0),plot=FALSE)
F2<-synth(f=22050,a=0.8,cf=4000,am=c(50,10),d=1,fm=c(500,5,0),plot=FALSE)
F3<-synth(f=22050,a=0.6,cf=6000,am=c(50,10),d=1,fm=c(500,5,2000),plot=FALSE)
F4<-synth(f=22050,a=0.4,cf=8000,am=c(50,10),d=1,fm=c(500,5,2000),plot=FALSE)
final1<-F1+F2+F3+F4
spectro(final1,f=22050,wl=512,ovlp=75,osc=TRUE)
title(main="synthesis of a AM/FM sound")

# 3D spectrogram of a tropical sparrow  song
spectro3D(tico,f=22050,wl=512,ovlp=75,zp=16,maga=2)

par(op1)