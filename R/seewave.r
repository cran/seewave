################################################################################
# Seewave v. 1.0. by Jérôme Sueur, Caroline Simonis-Sueur & Thierry Aubin
# Acknowledgements: Michel Baylac, Emmanuel Paradis, Martin Maechler 
# Requires R>=2.2.0, rgl
################################################################################

################################################################################
#                                AMA                                         
################################################################################

ama<-function(
wave,
f,
wl = wl,
plot = TRUE,
...)

{
env<-oscillo(wave, f=f, env = TRUE, plot = FALSE)
meanspec(env, f=f, wl=wl, plot=plot,...)
}


################################################################################
#                                ATTENUATION                                         
################################################################################

attenuation<-function
(
lref,
dref=1,
dstop,
n,
plot = TRUE,
xlab = "Distance (m)",
ylab = "dB",
...
)

{
data<-numeric(n)
step<-seq(dref,dstop,length.out=n)
{for(i in step){data[which(step==i)]<-lref-(20*log10(i/dref))}}

plot(x=step,y=data,xlab=xlab,ylab=ylab,...)
if (plot == FALSE) return(data)
}


################################################################################
#                                AUTOC                                         
################################################################################

autoc<-function(
wave,
f,
wl,
fmin,
threshold = FALSE,
plot = TRUE,
xlab = "Time (s)",
ylab = "Frequency (kHz)",
ylim = c(0,f/2000),
...)

{
cat("please wait...")
if (.Platform$OS.type == "windows") flush.console()
  
n<-nrow(wave)
step<-seq(1,n-wl,wl-(wl/100))
fmini<-round(wl*(fmin/(f/2)))

if (threshold)
  {
  thres<-max(abs(wave))*(threshold/100)
  wave1<-ifelse(abs(wave)<=thres,yes=0,no=1)
  wave2<-as.matrix(wave[,1]*wave1[,1])
  wave<-as.matrix(wave2)
  }

# discards the two last windows because of the lag
N<-length(step)-1  
R<-matrix(data=numeric(wl*N),wl,N)

for (j in 1:N)
  {
  for (i in 1:wl)
    {  
    R[i,j]<-(1/(2*wl)+1)*sum(wave[step[j]:(step[j]+wl-1),1]
                              *wave[(step[j]+i):(step[j]+wl+i-1),1])
    }
  }

tfond<-numeric(N)
for (k in 1:N) {tfond[k]<-which.max(R[-c(fmini:wl),k])}
y0<-f/tfond/1000
y<-ifelse(y0==f/1000,yes=NA,no=y0)
  
if (plot == TRUE)
  {
  x<-seq(0,n/f,length.out=N)
  plot(x=x, y=y,
  xlab = xlab,
  ylab = ylab, ylim = ylim,
  las =1,
  ...)
  }

else return(y)

}



################################################################################
#                                CONVSPL                                         
###############################################################################

convSPL<-function
(
x,
d = 1,
Iref = 10^-12,
pref = 2*10^-5 
)

{
P<-4*pi*(d^2)*Iref*(10^(x/10))
I<-Iref*(10^(x/10))
p<-pref*(10^(x/20))
conv<-list(P = P, I = I, p = p)
return(conv)
}


################################################################################
#                                CORENV                                         
################################################################################

corenv<-function(
wave1,
wave2,
f,
smooth =20,
plot = TRUE,
plotval = TRUE,
method = "spearman",
col = "black",
colval = "red",
cexval = 1,
fontval = 1,
xlab = "Time (s)",
ylab = "Coefficient of correlation (r)",
...)

{
n<-nrow(wave1)

if (smooth == 0 | smooth == FALSE) smooth<-1

if (smooth<10)
  {
  cat("please wait...")
  if (.Platform$OS.type == "windows") flush.console()
  }
  
x<-oscillo(wave=wave1,f=f,env=TRUE,smooth=smooth,plot=FALSE)
y<-oscillo(wave=wave2,f=f,env=TRUE,smooth=smooth,plot=FALSE)

nx<-length(x)
ny<-length(y)

if (nx != ny) stop("'wave 1' and 'wave 2' must have the same length")

meanx<-mean(x)
meany<-mean(y)
diffx<-x-meanx
r1<-numeric(nx)
r2<-numeric(nx)

for (i in 0:(nx-1))
{
r1[i+1]<-cor(x=x,y=c(y[(i+1):ny],rep(0,i)),method = method)
}

for (i in 0:(nx-1))
{
r2[i+1]<-cor(x=x,y=c(rep(0,i),y[1:(ny-i)]),method = method)
}

r2<-r2[-1]
r<-c(rev(r1),r2)
rmax<-max(r)
offset<-which.max(r)
if (offset<=(length(r)/2)) {offsetp<-which.max(r)} else {offsetp<-which.max(r)-1}
offsett<-((offsetp*smooth)-n)/f

if (offsetp < nx)
      {
      p<-cor.test(x=x, y=c(y[(nx-offsetp+1):ny],rep(0,nx-offsetp)),
      method = method)
      }
else
      {
      p<-cor.test(x=x, y=c(rep(0,offsetp-nx),y[1:(ny-(offsetp-nx))]),
      method = method)
      }
p<-p$p.value


if (plot == TRUE)
  {
  X<-seq(-n/f,n/f,length.out=2*nx-1)
  plot(x = X, y = r, xlab = xlab, ylab = ylab, col = col,...)
  if (plotval==TRUE)
    {
    mtext(paste(
      "rmax = ", as.character(round(rmax,2)),
      ", offset = ", as.character(round(offsett,3)), "s", sep=" "),
      side=3, line=-2, col=colval, cex=cexval, font=fontval)
    segments(x0=round(offsett,3),y0=min(r), x1=round(offsett,3), y1=rmax, col=colval, lty=2)
    segments(x0=-n/f,y0=rmax, x1=round(offsett,3), y1=rmax, col=colval, lty=2)
    points(x=round(offsett,3),y=rmax,pch=19,cex=1, col=colval)
  }
  }

else
  {
  corr<-list(r = r, r = rmax, p = p, t = offsett)
  return(corr)
  }

}



################################################################################
#                                CORSPEC                                         
################################################################################

corspec<-function(
x,
y,
range,
plot = TRUE,
plotval = TRUE,
method = "spearman",
col = "black",
colval = "red",
cexval = 1,
fontval = 1,
xlab = "Frequency (kHz)",
ylab = "Coefficient of correlation (r)",
...)

{
nx<-length(x)
ny<-length(y)

if (nx != ny) stop("x and y must have the same length")
if(any(x < 0) | any(x < 0)) stop("data does not have to be in dB")

meanx<-mean(x)
meany<-mean(y)
diffx<-x-meanx
r1<-numeric(nx)
r2<-numeric(nx)

for (i in 0:(nx-1))
{
r1[i]<-cor(x=x,y=c(y[(i+1):ny],rep(0,i)),method = method)
}

for (i in 0:(nx-1))
{
r2[i+1]<-cor(x=x,y=c(rep(0,i),y[1:(ny-i)]),method = method)
}

r2<-r2[-1]
r<-c(rev(r1),r2)
rmax<-max(r)
offset<-which.max(r)
if (offset<=(length(r)/2)) {offsetp<-which.max(r)} else {offsetp<-which.max(r)-1}
offsetf<-((range[2]-range[1])*(offsetp-nx))/nx

if (offsetp < nx)
      {
      p<-cor.test(x=x, y=c(y[(nx-offsetp+1):ny],rep(0,nx-offsetp)),
      method = method)
      }
else
      {
      p<-cor.test(x=x, y=c(rep(0,offsetp-nx),y[1:(ny-(offsetp-nx))]),
      method = method)
      }
p<-p$p.value    


if (plot == TRUE)
  {
  X<-seq(-range[2],range[2],length.out=2*nx-1)
  plot(x = X, y = r, xlab = xlab, ylab = ylab, col = col,...)
  if (plotval==TRUE)
    {
    mtext(paste(
      "rmax = ", as.character(round(rmax,2)),
      ", offset = ", as.character(round(offsetf,2)), "kHz", sep=" "),
      side=3, line=-2, col=colval, cex=cexval, font=fontval)
    segments(x0=offsetf,y0=min(r), x1=offsetf, y1=rmax, col=colval, lty=2)
    segments(x0=-range[2],y0=rmax, x1=offsetf, y1=rmax, col=colval, lty=2)
    points(x=offsetf,y=rmax,pch=19,cex=1, col=colval)
  }
  }
  
else  
  {
  corr<-list(r = r, rmax = rmax, p = p, f = offsetf)
  return(corr)
  }

}



################################################################################
#                                COVSPECTRO                                         
################################################################################


covspectro<-function
(
wave1,
wave2,
f,
wl,
n,
plot = TRUE,
plotval = TRUE,
method = "spearman",
col = "black",
colval = "red",
cexval = 1,
fontval = 1,
xlab = "Time (s)",
ylab = "Normalised covariance (cov)",
...
)

{
if (n>21)
  {
  cat("please wait...")
  if (.Platform$OS.type == "windows") flush.console()
  }
  
if((n <- as.integer(n)) %% 2 != 1) stop("'n' must be odd")

n<-(n%/%2)+1

n1<-nrow(wave1)
n2<-nrow(wave2)

if (n1 != n2) stop("'wave 1' and 'wave 2' must have the same length")

step1<-seq(1,n1-wl,wl); lstep1<-length(step1)
step2<-round(seq(1,n2,length.out=n))

# wave not time shifted
spectro1<-sspectro(wave=wave1,f=f,wl=wl)

spectro2a<-array(numeric((wl/2)*lstep1*n),dim=c((wl/2),lstep1,n))
spectro2b<-array(numeric((wl/2)*lstep1*n),dim=c((wl/2),lstep1,n))
cov1<-numeric(n)
cov2<-numeric(n)

# wave time shifted
# successive spectrograms
# covariance of spectrogram1/spectra1 with spectrogram2/spectra1,
# spectrogram1/spectra2 with spectrogram2/spectra2 and so on
# diagonal of the cov matrix and mean of this diagonal
# one mean cov for comparaison between 2 spectrogramsfor (i in step2)
for (i in step2)
{
spectro2a[,,which(step2==i)]<-sspectro(wave=as.matrix(c(wave2[i:n2],rep(0,i-1))),f=f,wl=wl)
spectro2a<-ifelse(spectro2a=="NaN",yes=0,no=spectro2a)
cov1[which(step2==i)]<-mean(diag(cov(x=spectro1,y=spectro2a[,,which(step2==i)],method = method)))
}

for (i in step2)
{
spectro2b[,,which(step2==i)]<-sspectro(wave=as.matrix(c(rep(0,i),wave2[1:(n2-i)])),f=f,wl=wl)
spectro2b<-ifelse(spectro2b=="NaN",yes=0,no=spectro2b)
cov2[which(step2==i)]<-mean(diag(cov(x=spectro1,y=spectro2b[,,which(step2==i)],method = method)))
}

# to normalise the covariance we need covmax that is the autocovariance of spectro1
covmax<-mean(diag(cov(x=spectro1,y=spectro1,method = method)))

# discard the first value of cov2 that is already computed in cov1
cov2<-cov2[-1]
cov3<-c(rev(cov1),cov2)
cov4<-cov3/covmax
cov4max<-max(cov4)
offset<-which.max(cov4)
offsetp<-which.max(cov4)
offsett<-(((offsetp*n1)/n)-n1)/f

if (plot == TRUE)
  {
  x<-seq(-n1/f,n1/f,length.out=(2*n)-1)
  plot(x = x, y = cov4, xlab = xlab, ylab = ylab, col = col,...)
  if (plotval == TRUE)
    {
    mtext(paste(
      "covmax = ", as.character(round(cov4max,2)),
      ", offset = ", as.character(round(offsett,3)), "s", sep=" "),
      side=3, line=-2, col=colval, cex=cexval, font=fontval)
    segments(x0=round(offsett,3),y0=min(cov4), x1=round(offsett,3), y1=cov4max, col=colval, lty=2)
    segments(x0=-n1/f,y0=cov4max, x1=round(offsett,3), y1=cov4max, col=colval, lty=2)
    points(x=round(offsett,3),y=cov4max,pch=19,cex=1, col=colval)
    }
  }

else
  {
  covar<-list(cov = cov4, covmax = cov4max, t = offsett)
  return(covar)
  }

}



################################################################################
#                                CUTW                                         
################################################################################

cutw<-function(
wave,
f,
from=FALSE,
to=FALSE,
plot = FALSE,
marks = TRUE,
...)

{
if (from == 0) {a<-1; b<-round(to*f)}
if (from == FALSE) {a<-1; b<-round(to*f);from<-0}
if (to == FALSE) {a<-round(from*f); b<-nrow(wave);to<-nrow(wave)/f}
else {a<-round(from*f); b<-round(to*f)}
wavecut<-as.matrix(wave[a:b,])
  
if (plot == TRUE)
  {
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))
  par(mfrow=c(2,1))
  oscillo(wave,f=f,k=1,j=1,...)
  title(main="original")
    if (marks == TRUE)
      {
      abline(v=from, col="red", lty=2)
      abline(v=to, col="red", lty=2)
      }
  oscillo(wavecut,f=f,k=1,j=1,...)
  title(main="cut")
  }
else return(wavecut)
}

 
################################################################################
#                                DBSCALE                                        
################################################################################

dBscale<-function
(
collevels,
palette = spectro.colors,
side = 4,
textlab = "Amplitude\n(dB)",
cexlab = 0.75,
fontlab = 1,
collab = "black",
...
)

{
plot.new()
levels<-collevels
col <- palette(length(collevels) - 1)
par(las=1,...)
    
if (side == 2 | side == 4)
    {    
    plot.window(xlim = c(0, 1), ylim = range(collevels), xaxs = "i",
        yaxs = "i")
    mtext(textlab, side=3, outer=FALSE, line=1.5, adj=0, font=fontlab, cex=cexlab, col=collab)
    rect(xleft=0, ybottom=levels[-length(levels)], xright=0.95, ytop=levels[-1], col = col, lty=0)
    segments(x0=0,y0=0,x1=0.95,y1=0)
    segments(x0=0,y0=min(collevels),x1=0.95,y1=min(collevels))
          
    if (side == 2)
      {
      axis(2,...)
      abline(v=c(0,0.95))
      }
      
    if (side == 4)
      {
      axis(4,pos=0.95,...)
      abline(v=c(0,0.95))
      }
    }

if (side == 1  | side == 3)
    {    
    plot.window(xlim = range(collevels), ylim = c(0, 1), xaxs = "i",
        yaxs = "i")
    mtext(textlab, side=3, outer=FALSE, line=1.5, adj=0, font=fontlab, cex=cexlab, col=collab)
    rect(xleft=levels[-length(levels)], ybottom=0, xright=levels[-1], ytop=0.95, col = col, lty=0)
    segments(x0=0,y0=0,x1=0,y1=0.95)
    segments(x0=min(collevels),y0=0,x1=min(collevels),y1=0.95)
          
    if (side == 1)
      {
      axis(1,...)
      abline(h=0.95)
      }
    if (side == 3)
      {
      axis(3,pos=0.95,...)
      abline(h=0)
      }
    }    
}    


################################################################################
#                                DELETEW
################################################################################

deletew<-function(
wave,
f,
from=FALSE,
to=FALSE,
plot = FALSE,
marks = TRUE,
...)

{
if (from == 0) {a<-1; b<-round(to*f)}
if (from == FALSE) {a<-1; b<-round(to*f);from<-0}
if (to == FALSE) {a<-round(from*f); b<-nrow(wave);to<-nrow(wave)/f}
else {a<-round(from*f); b<-round(to*f)}
wavecut<-as.matrix(wave[-(a:b),])

if (plot == TRUE)
  {
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))
  par(mfrow=c(2,1))
  oscillo(wave,f=f,k=1,j=1,...)
  title(main="original")
    if (marks == TRUE)
      {
      abline(v=from, col="red", lty=2)
      abline(v=to, col="red", lty=2)
      }
  oscillo(wavecut,f=f,k=1,j=1,...)
  title(main="after deletion")
  }
else return(wavecut)
}


################################################################################
#                                DFREQ                                         
################################################################################

dfreq<-function(
wave,
f,
wl,
ovlp = 0,
threshold = FALSE,
plot = TRUE,
xlab = "Times (s)",
ylab = "Frequency (kHz)",
ylim = c(0,f/2000),
...)

{
n<-nrow(wave)
step<-seq(1,n-wl,wl-(ovlp*wl/100))

if (threshold)
  {
  thres<-max(abs(wave))*(threshold/100)
  wave1<-ifelse(abs(wave)<=thres,yes=0,no=1)
  wave2<-as.matrix(wave[,1]*wave1[,1])
  wave3<-ifelse(wave2==0,yes=1e-6,no=wave2)
  wave4<-as.matrix(wave3)
  }
  
else wave4<-wave

y1<-matrix(data=numeric(wl*length(step)),wl,length(step))
for(i in step)
{y1[,which(step==i)]<-Mod(fft(wave4[i:(wl+i-1),]*hanning.w(wl)))}
		
y2<-y1[1:(wl/2),]				

y3<-matrix(data=numeric(length(step)*2),length(step),2)
for (i in 1:length(step))
# [1,1] is to keep only the first line and firt column of the results (=c(1,1))
# when there is no max value, this happens when the fft is totatally flat
# (e. g. signal =0)
{y3[i,]<-as.numeric(which(y2==max(y2[,i]),arr.ind=TRUE)[1,1])}
y3<-(f*y3[,1])/(1000*wl)
# discards max results when signal = 0, i. e. when which.max = c(1,1)
y<-ifelse(y3==f/(wl*1000), yes=NA, no=y3)

if (plot==TRUE)
  {
	x<-seq(0,n/f,length.out=length(step))
	plot(x, y,
	xaxs = "i", xlab = xlab,
  yaxs = "i", ylab = ylab, ylim = ylim,
  ...)
	}
else
return(y)
}



################################################################################
#                               EXPORT                                        
################################################################################

export<-function(
wave,
f,
file = "newwave.txt",
...)

{
n<-nrow(wave)
header<-paste("[ASCII ",f,"Hz, Channels: 1, Samples: ",n,", Flags: 0]", sep="")
write.table(x=wave, file=file, row.names=FALSE, col.names=header, ...)
}



################################################################################
#                               MEANSPEC                                        
################################################################################

meanspec<-function(
wave,
f,
wl,
ovlp = 0,
PSD =  FALSE,
dB = FALSE,
from = FALSE,
to = FALSE,
peaks = FALSE,
identify = FALSE,
col = "black",
cex = 1,
colpeaks = "red",
cexpeaks = 1,
fontpeaks = 1,
plot = 1,
flab = "Frequency (kHz)",
alab = "Amplitude",
flim = c(0,f/2000),
alim = NULL,
...)

{

  
if (from|to)
  {
  if (from == 0) {a<-1; b<-round(to*f)}
  if (from == FALSE) {a<-1; b<-round(to*f);from<-0}
  if (to == FALSE) {a<-round(from*f); b<-nrow(wave);to<-nrow(wave)/f}
  else {a<-round(from*f); b<-round(to*f)}
  wave<-as.matrix(wave[a:b,])
  }

n<-nrow(wave)
N<-wl
step<-seq(1,n-wl,wl-(ovlp*wl/100))		# FT windows
y1<-matrix(data=numeric((wl)*length(step)),wl,length(step))

for(i in step)
   {y1[,which(step==i)]<-
            Mod(fft(c(wave[i:(wl+i-1),])*hanning.w(wl)))}
  
y2<-y1[1:(wl/2),]	# to keep only the relevant frequencies (half of the FT)
y3<-y2/max(y2)					      # to get only values between 0 and 1
y4<-apply(y3,MARGIN=1,mean)   # mean computation (by rows)
y5<-y4/max(y4)
# replaces 0 values in spectra that can't be processed by the following log10())
y<-ifelse(y5==0,yes=1e-6,no=y5)

if (PSD == TRUE) y<-y^2

if (peaks)
  {
  check.pks(y)
  p<-peaks(y,peaks)
  respeaks<-seq(y)[p]*f/N/1000
  }

x<-seq((f/1000)/N,f/2000,length.out=N/2)  
if (is.null(alim)==TRUE)
    {if (dB==TRUE) alim<-c(min(20*log10(y)),20) else alim<-c(0,1.1)}

if(plot == 1) # plots x-y graph with Frequency as X-axis
	{
    if(dB == TRUE)
	  {
    y<-20*log10(y)	
	  plot(x,y,
		xaxs = "i", xlab = flab, xlim = flim,
		yaxs = "i", yaxt = "s", ylab = alab, ylim = alim,
		cex = cex, col = col,
    las=1,...)
    }
    else
    {
    plot(x,y,
		xaxs = "i", xlab = flab, xlim = flim,
		yaxs = "i", yaxt = "n", ylab = alab, ylim = alim,
		cex = cex, col = col,
    las = 1,...)
    }

 	  if(identify == TRUE)
    {
    cat("choose points on the spectrum")
    if (.Platform$OS.type == "windows") flush.console()
    id<-identify(x=x,y=y,labels=round(x,2),tolerance=0.15)
    return(((id*f)/N)/1000)
    }     	
    
    if(peaks)                              
    {
    if (dB == TRUE)
    text(seq(y)[p]*f/N/1000, y[p]+5,
            as.character(round(seq(y)[p]*f/N/1000,3)),
            col = colpeaks, cex = cexpeaks, font = fontpeaks)
    else  
    text(seq(y)[p]*f/N/1000, y[p]+0.05,
            as.character(round(seq(y)[p]*f/N/1000,3)),
            col = colpeaks, cex = cexpeaks, font = fontpeaks)
    }
  }

if(plot == 2) # plots x-y graph with Frequency as Y-axis
	{
    if(dB == TRUE)
	  {
    y<-20*log10(y)	
	  plot(y,x,
		xaxs = "i", xlab = alab, xlim = alim,
		yaxs = "i", yaxt = "s", ylab = flab, ylim = flim,
    cex = cex, col = col,
    las = 1,...)
    }
    else
    {
    plot(y,x,
		xaxs = "i", xaxt = "n", xlab = alab, xlim = c(0,1.2),
		yaxs = "i", ylab = flab, ylim = flim,
    cex = cex, col = col,
    las = 1,...)
    }
	  
    if(identify == TRUE)
    {
    cat("choose points on the spectrum")
    if (.Platform$OS.type == "windows") flush.console()
    id<-identify(x=y,y=x,labels=round(x,2),tolerance=0.15)
    return(((id*f)/N)/1000)
    }    
    		
    if(peaks)
    {
    if (dB == TRUE)
    text(y[p]+10, seq(y)[p]*f/N/1000,
          as.character(round(seq(y)[p]*f/N/1000,3)),
          col = colpeaks, cex = cexpeaks)
    else  
    text(y[p]+0.1, seq(y)[p]*f/N/1000,
          as.character(round(seq(y)[p]*f/N/1000,3)),
          col = colpeaks, cex = cexpeaks)
    }
  }

if(plot == FALSE) 
  {
  if(dB == TRUE) y<-20*log10(y)
  spec<-y[1:length(y)]	
  if (peaks)
      {
      results<-list(spec = spec ,peaks = respeaks)
      return(results)
      }
  else return(spec)
  }
}


################################################################################
#                                MOREDB                                        
################################################################################

moredB<-function
(
x
)

{
if (is.matrix(x)==TRUE) n<-nrow(x)
if (is.vector(x)==TRUE) n<-length(x)
if (is.numeric(x)==TRUE) n<-length(x)

data1<-as.numeric(x/10)
data2<-as.numeric(n)
for(i in seq(1,n)) {data2[i]<-10^data1[i]}
data3<-10*log10(sum(data2))
return(data3)
}


################################################################################
#                                MUTE                                        
################################################################################

mute<-function
(wave,
f,
from = FALSE,
to = FALSE,
plot = TRUE,
...
)

{

n<-nrow(wave)
wave.muted<-as.matrix(rep(0,n))

if (from | to)
  {
  if (from == FALSE)
      {
      b<-round(to*f)
      wave.muted<-as.matrix(c(rep(0,b),wave[(b+1):n,]))
      }

  if (to == FALSE)
      {
      a<-round(from*f)
      wave.muted<-as.matrix(c(wave[1:(a-1),],rep(0,length(a:n))))
      }

  if (from != FALSE & to!= FALSE)
      {
      if (from > to) stop("'from' must be inferior to 'to'")
      if (from == 0) {a<-1; b<-round(to*f)}
      else {
      a<-round(from*f)
      b<-round(to*f)}
      wave.muted<-as.matrix(c(wave[1:(a-1),],rep(0,length(a:b)),wave[(b+1):n,]))
    }
  }
if (plot == TRUE) {oscillo(wave.muted,f=f,...)} else return(wave.muted)

}


################################################################################
#                                NOISE                                        
################################################################################

noise<-function(
f,
d,
...
)

{
as.matrix(rnorm(d*f))
}


################################################################################
#                                OSCILLO                                        
################################################################################

oscillo<-function
(
wave,
f,
from = FALSE,
to = FALSE,
zoom = FALSE,
k=1,
j=1,
labels = TRUE,
byrow = TRUE, 
env = FALSE,
smooth= 0,
identify = FALSE,
plot = TRUE,
colwave = "black",
colbg = "white",
coltitle = "black",
collab = "black",
cexlab = 1,
fontlab = 1,
colline = "black",
colaxis = "black",
coly0 = "grey47",
title = FALSE,
bty = "l"
)

{
p<-k*j
  
if (from|to)
  {
  if (from == 0) {a<-1; b<-round(to*f)}
  if (from == FALSE) {a<-1; b<-round(to*f);from<-0}
  if (to == FALSE) {a<-round(from*f); b<-nrow(wave);to<-nrow(wave)/f}
  else {a<-round(from*f); b<-round(to*f)}
  wave<-as.matrix(wave[a:b,])
  n<-nrow(wave)
  }

else {n<-nrow(wave); a<-0; from<-0; to<-n/f; wave<-as.matrix(wave)}

if (env == TRUE)
  {
    if (smooth)
    {
    z0<-abs(wave[0:n,])
    step<-seq(1,n-smooth,smooth)
    z1<-numeric(length(step)) 
    for(i in step)
      {z1[which(step==i)]<-mean(z0[i:(i+smooth)])}
    wave<-as.matrix(z1)
    n<-nrow(wave)
    x<-n%/%p
    f<-f/smooth
    }
    else
    wave<-as.matrix(abs(wave[0:n,]))
  if (plot == FALSE) return (wave)
  }

# to get a single window view
if (plot == TRUE)
{
if (k == 1 & j == 1)
{
  if (zoom == TRUE)
    {	
    par(tcl=0.5, bg=colbg, col.axis=colaxis, col=colline,las=0)
    plot(x=seq(from,to,length.out=n), y=wave,
          col=colwave, type="l",
		      xaxs="i", yaxs="i",
		      xlab="", ylab="",
		      yaxt="n",
		      cex.lab=0.8, font.lab=2,
          bty=bty
		      )
    if (bty == "l" | bty == "o")
          {axis(side=1, col=colline,labels=FALSE)
          axis(side=2, at=max(wave,na.rm=TRUE), col=colline,labels=FALSE)}
	   mtext("Time (s)",col=collab, font=fontlab,side=1,line=3,cex=cexlab)
	   mtext("Amplitude",col=collab, font=fontlab, cex=cexlab,side=2,line=2.5)
	   abline(h=0,col=coly0,lty=2)
    
    cat("choose start and end positions on the wave")
    if (.Platform$OS.type == "windows") flush.console()
    coord<-locator(n=2)
    from<-coord$x[1]; c<-from*f-a
    to<-coord$x[2]; d<-to*f-a
    if (d<c) {c<-d; d<-c}
    wave<-as.matrix(wave[c:d,1])
    n<-nrow(wave)
    }
  
  par(tcl=0.5, bg=colbg, col.axis=colaxis, col=colline,las=0)
	plot(x=seq(from,to,length.out=n), y=wave,
		col=colwave, type="l",
		xaxs="i", yaxs="i",
		xlab="", ylab="",
		yaxt="n",
		cex.lab=0.8, font.lab=2,
    bty=bty
		)
	
  if (bty == "l" | bty == "o")
      {axis(side=1, col=colline,labels=FALSE)
      axis(side=2, at=max(wave,na.rm=TRUE), col=colline,labels=FALSE)}
	
	if (labels == TRUE)
      { 
      mtext("Time (s)",col=collab, font=fontlab,side=1,line=3,cex=cexlab)
      mtext("Amplitude",col=collab, font=fontlab, cex=cexlab,side=2,line=3)
	    }
	    
  abline(h=0,col=coly0,lty=2)

  if (identify == TRUE)
      {
      cat("choose points on the wave")
      if (.Platform$OS.type == "windows") flush.console()
      x<-seq(from=from,to=to,length.out=n)
      y<-wave
      id<-identify(x=x, y=y, labels=round(x,3), col = "red", plot= TRUE)
      abline(v=(id/f)+from,col="red")
      return (round((id/f)+from,3))
      }
}

# to get a multi-window view
else
{
  if (zoom == TRUE)
  stop ("'zoom' does work with a single-frame window only ('k'=1 and 'j'=1)")
  if (identify == TRUE)
  stop ("'identify' does work with a single-frame window only ('k'=1 and 'j'=1)")
  x<-n%/%p
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))
  m<-matrix(1:p,k,j,byrow=byrow)
	layout(m)
	par(tcl=0.5,oma=c(3,2,2,0.5),
      mar=rep(0,4)+0.8, mgp=c(0,0.15,0),
      bg=colbg, col.axis=colaxis, col=colline, las=0)

# plots the first window
	wave1<-as.matrix(wave[0:x,]); n1<-nrow(wave1)
	plot(x=seq(from,from+(x/f),length.out=n1), y=wave1,
		col=colwave, type="l",
		xaxs="i", yaxs="i",
		xlab="", ylab="", ylim=range(wave),
		yaxt="n",bty=bty)
	axis(side=1, col=colline,labels=FALSE)
	if (bty == "l" | bty == "o")
        {axis(side=2, at=max(wave), col=colline,labels=FALSE)
        axis(side=1, col=colline,labels=FALSE)}
  abline(h=0,col=coly0,lty=2)
  
# title
if (title == FALSE)
	{title <- paste("")}
else 	
	{
	title<-paste("Window time =",
                as.character(round(n/(p*f),3)),"s - Total time =",
                as.character(round(n/f,3)), "s - Fs =",
                as.character(f),"Hz")
	mtext(title, col=coltitle, side=3,line=0.4, cex=1, outer=TRUE)
	}

# X-Y labels
if (labels == TRUE)
  {
  mtext("Time (s)",col=collab, side=1,line=1.5, font=fontlab,cex=cexlab,outer=TRUE)
	mtext("Amplitude",col=collab, side=2, font=fontlab,cex=cexlab,
        line=0.4,outer=TRUE)
  }
  
# plots following windows
for(i in 1:(p-1))
  {
	xx<-((i*n)%/%p)+1
	yy<-((i+1)*n)%/%p
	wave2<-as.matrix(wave[xx:yy,]); n2<-nrow(wave2)
	plot(x=seq(from+(xx/f),from+(yy/f),length.out=n2), y=wave2,
		col=colwave, type="l",
		xaxs="i", yaxs="i",
		xlab="", ylab="", ylim = range(wave),
		yaxt="n", bty=bty)

	if (bty == "l" | bty == "o")
        {axis(side=2, at = max (wave), col=colline,labels=FALSE)
       	axis(side=1, col=colline,labels=FALSE)}
  abline(h=0,col=coly0,lty=2)	
  }
}
}
else return (wave)
}


################################################################################
#                                PASTEW                                         
################################################################################

pastew<-function(
wave1,
wave2,
f,
at = FALSE,
plot = FALSE,
marks = TRUE,
...)

{
wave1<-as.matrix(wave1)
wave2<-as.matrix(wave2)

if (at)
  {
  pos<-round(at*f)
  wave2a<-as.matrix(wave2[c(1:pos),1])
  wave2b<-as.matrix(wave2[c(pos:nrow(wave2)),1])
  wave3<-rbind(wave2a,wave1,wave2b)
  }
else
  {wave3<-rbind(wave1,wave2)}

if (plot == TRUE)
  {
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))
  par(mfrow=c(3,1))
  oscillo(wave1,f=f,k=1,j=1)
  title(main="signal to be pasted")
  oscillo(wave2,f=f,k=1,j=1)
  title(main="signal to be completed")  
  oscillo(wave3,f=f,k=1,j=1)
  title(main="resulting signal")
    if (marks == TRUE)
    {
    abline(v=at, col="red", lty=2)
    abline(v=at+(nrow(wave1))/f, col="red", lty=2)
    }
  }
else return(wave3)
}


################################################################################
#                                Q                                         
################################################################################

Q<-function(
x,
range,
level = -3,
plot = TRUE,
colval = "red",
cexval = 1,
fontval = 1,
flab = "Frequency (kHz)",
alab = "Relative amplitude (dB)",
...)

{
spectrum<-x
if (max(spectrum) == 1) stop ("data must be in dB")
if (which.max(spectrum) == 1) 
    stop ("maximal peak does not have to be the first value of spectrum") 

n0<-length(spectrum)

spec1<-approx(spectrum,n=102400)
spec1<-as.matrix(spec1$y)
n1<-nrow(spec1)
level2<-round(max(spec1[,1]),1)+level

f0<-which.max(spec1[,1])
specA<-as.matrix(spec1[1:f0,1])
nA<-nrow(specA)
specB<-as.matrix(spec1[f0:nrow(spec1),1])
f1<-which(round(specA,1) == level2)
f1khz<-((f1[length(f1)]/n1)*(range[2]-range[1]))+range[1]
f2<-which(round(specB,1) == level2)+(nA-1)
f2khz<-((f2[1]/n1)*(range[2]-range[1]))+range[1]

Q<-f0/(f2[1]-f1[length(f1)])

# plot based on orignal data (=> spectrum)
if (plot == TRUE)
	{
	x<-seq(range[1],range[2],length.out=n0)
	y<-spectrum
	plot(x=x,y=y,xlab=flab,ylab=alab,...)
	arrows(f1khz,level2,f2khz,level2,length=0.1,col=colval,code=3,angle=15)
	text(paste("Q =",as.character(round(Q,2))),x=f2khz,y=level2,pos=4,
      col=colval, cex=cexval, font=fontval)
	}

else return(Q)

}


################################################################################
#                                RESAMP                                         
################################################################################


resamp<-function(
wave,
f,
g
)

{
n<-nrow(wave)
if (g==f) stop ("'f' and 'g' must be different")
if (g<f)  {r<-f/g; wave1<-wave[seq(1,n,by=r),1]}
if (g>f)  {s<-(n*g)/f; wave1<-approx(wave,n=s)$y}
return(as.matrix(wave1))
}


################################################################################
#                                SPEC                                         
################################################################################

spec<-function(
wave,
f,
wl = FALSE,
PSD = FALSE,
dB = FALSE,
at = FALSE,
from = FALSE,
to = FALSE,
peaks = FALSE,
identify = FALSE,
col = "black",
cex = 1,
colpeaks = "red",
cexpeaks = 1,
fontpeaks = 1,
plot = 1,
flab = "Frequency (kHz)",
alab = "Amplitude",
flim = c(0,f/2000),
alim = NULL,
...)

{
if (from|to)
  {
  if (from == 0) {a<-1; b<-round(to*f)}
  if (from == FALSE) {a<-1; b<-round(to*f);from<-0}
  if (to == FALSE) {a<-round(from*f); b<-nrow(wave);to<-nrow(wave)/f}
  else {a<-round(from*f); b<-round(to*f)}
  wave<-as.matrix(wave[a:b,])
  }
  
if (at)
  {
  c<-round(at*f)
  wl2<-wl%/%2
  wave<-as.matrix(wave[(c-wl2):(c+wl2),])
  n<-nrow(wave)
  }

n<-nrow(wave)
x<-seq((f/1000)/n,f/1000,length.out=n)
wave<-wave*hanning.w(n)
y1<-Mod(fft(wave[,1]))
y2<-y1/max(y1)
# replaces 0 values in spectra that can't be processed by the following log10()
y<-ifelse(y2==0,yes=1e-6,no=y2)

if (PSD == TRUE) y<-y^2

if (peaks)
  {
  check.pks(y)
  p<-peaks(y,peaks)
  respeaks<-seq(y)[p]*f/n/1000
  respeaks<-respeaks[1:(length(respeaks)/2)]
  }

if (is.null(alim)==TRUE)
  {if (dB==TRUE) alim<-c(min(20*log10(y)),20) else alim<-c(0,1.1)}
    
if(plot == 1) # plots x-y graph with Frequency as X-axis
	{
    if(dB == TRUE)
	  {
    y<-20*log10(y)	
	  plot(x=x,y=y,
		xaxs = "i", xlab = flab, xlim = flim,
		yaxs = "i", yaxt = "s", ylab = alab, ylim = alim,
		col = col, cex = cex,
    las = 1,
    ...)
    }
    else
    {
    plot(x=x,y=y,
		xaxs="i", xlab=flab, xlim = flim,
		yaxs="i", yaxt="n", ylab = alab, ylim=alim,
    col = col, cex = cex,
    las = 1,
    ...)
    }
	   	   
    if(identify == TRUE)
    {
    cat("choose points on the spectrum")
    if (.Platform$OS.type == "windows") flush.console()
    id<-identify(x=x,y=y,labels=round(x,2),tolerance=0.15)
    return(((id*f)/n)/1000)
    } 
	
    if(peaks)
    {
    if (dB == TRUE)
    text(seq(y)[p]*f/n/1000, y[p]+5,
              as.character(round(seq(y)[p]*f/n/1000,3)),
              col = colpeaks, cex = cexpeaks, font = fontpeaks)
    else  
    text(seq(y)[p]*f/n/1000, y[p]+0.05,
              as.character(round(seq(y)[p]*f/n/1000,3)),
              col = colpeaks, cex = cexpeaks, font = fontpeaks)
    }
  }

if(plot == 2) # plots x-y graph with Frequency as Y-axis
	{
    if(dB == TRUE)
	  {
    y<-20*log10(y)	
	  plot(x=y,y=x,
		xaxs = "i", xlab = alab, xlim = alim,
		yaxs = "i", yaxt = "s", ylab = flab, ylim = flim,
    col = col, cex = cex,
    las = 1,
    ...)
    }
    else
    {
    plot(x=y,y=x,
		xaxs = "i", xaxt = "n", xlab = alab, xlim = alim,
		yaxs = "i", ylab = flab, ylim = flim,
    col = col, cex = cex,
    las = 1,
    ...)
    }
		

    if(identify == TRUE)
    {
    cat("choose points on the spectrum")
    if (.Platform$OS.type == "windows") flush.console()
    id<-identify(x=y,y=x,labels=round(x,2),tolerance=0.15)
    return(((id*f)/n)/1000)
    } 
		
    if(peaks)
    {
    if (dB == TRUE)
    text(y[p]+10, seq(y)[p]*f/n/1000,
              as.character(round(seq(y)[p]*f/n/1000,3)),
              col = colpeaks, cex = cexpeaks, font= fontpeaks)
    else  
    text(y[p]+0.1, seq(y)[p]*f/n/1000,
              as.character(round(seq(y)[p]*f/n/1000,3)),
              col = colpeaks, cex = cexpeaks, font= fontpeaks)
    }
  }

if(plot == FALSE) 
  {
  if(dB == TRUE) y<-20*log10(y)
  spec<-y[1:(n/2)]	
  if (peaks)
      {
      results<-list(spec = spec ,peaks = respeaks)
      return(results)
      }
  else return(spec)
  }
}



################################################################################
#                                SPECTRO                                        
################################################################################

spectro<-function(
wave,
f,
wl,
zp=0,
ovlp=0,
plot = TRUE,
grid = TRUE,
osc = FALSE,
scale = TRUE,
cont = FALSE,
collevels = seq(-30,0,1),
palette = spectro.colors,
contlevels = seq (-30,0,10),
colcont = "black",
colgrid = "black",
plot.title =
    title(main = "", xlab = "Time (s)",
    ylab = "Frequency (kHz)"),
scalelab = "Amplitude\n(dB)",
scalefontlab = 1,
axisX = TRUE,
axisY = TRUE,
...)

{
n<-nrow(wave)
step<-seq(1,n-wl,wl-(ovlp*wl/100))		# FT windows

z1<-matrix(data=numeric((wl+(zp))*length(step)),wl+zp,length(step))
zpl<-zp/2
if(zpl==0)
  {
  for(i in step)
  {z1[,which(step==i)]<-Mod(fft(wave[i:(wl+i-1),]*hanning.w(wl)))}
  }

else
  {
  for(i in step)
  {z1[,which(step==i)]<-
  Mod(fft(c(1:zpl,wave[i:(wl+i-1),],1:zpl)*hanning.w(wl+zp)))}
  }	

# to keep only the relevant frequencies (half of the FT)
z2<-z1[1:((wl+zp)/2),]	
# to get only values between 0 and 1
z3<-z2/max(z2)					
# replaces 0 values in spectra (that can't be processed by the following log10())
z4<-ifelse(z3==0,yes=1e-6,no=z3)
# to get dB values
z<-20*log10(z4)[-1,]					  

# settings of X, Y, Z plot axis
X<-seq(0,n/f,length.out=length(step))
Y<-seq(0,f/2000,length.out=((wl+zp)/2)-1)
Z<-t(z)
   
if (plot==TRUE)
 	{
  Zlim<-range(Z, finite = TRUE) 
  if (osc==TRUE & scale==TRUE)
    {
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    layout(matrix(c(3, 2 ,1, 0), nc = 2, byrow=TRUE),
              widths = c(6, 1), heights=c(3,1))
    par(mar=c(5,3.7,1,0), las=0)
    soscillo(wave=wave,f=f,...)
    par(mar=c(5,4,3,0), las=0)
    filled.contour.modif1(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
			plot.title=plot.title,color.palette=palette,scalelab=scalelab,
      scalefontlab=scalefontlab,axisX=axisX,axisY=axisY)		
  	if (grid == TRUE) grid(nx=NA, ny=NULL, col=colgrid)
    }
  
  if (osc==FALSE & scale==TRUE)
    {
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    layout(matrix(c(2, 1), nc = 2), widths = c(5, 1))
    par(las=1)
    filled.contour.modif3(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
			plot.title=plot.title,color.palette=palette,scalelab=scalelab,
      scalefontlab=scalefontlab,axisX=axisX, axisY=axisY)		
    if (grid == TRUE) grid(nx=NA, ny=NULL, col=colgrid)
    }

  if (osc==TRUE & scale==FALSE)
    {
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    layout(matrix(c(2,1), nr = 2, byrow=TRUE), heights=c(3,1))
    par(mar=c(5,4,1,2)+0.1, las=0)
    soscillo(wave=wave,f=f,...)
    par(mar=c(2,4,2,2)+0.1)
    filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
			plot.title=plot.title, color.palette=palette,axisX=axisX , axisY=axisY)		
    if (grid == TRUE) grid(nx=NA, ny=NULL, col=colgrid)
    }
  
  if (osc==FALSE & scale==FALSE)
   {
   par(las=1)
   filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
			plot.title=plot.title, color.palette=palette, axisX=axisX , axisY=axisY)		
   if (grid == TRUE) grid(nx=NA, ny=NULL, col=colgrid)
   }

  if (cont==TRUE) 
	 {
   contour(X,Y,Z,add=TRUE,
	 levels=contlevels,nlevels=5,col=colcont,...)
	 }  
  }
else return(z)
}


################################################################################
#                                SPECTRO3D                                        
################################################################################

spectro3D<-function(
wave,
f,
wl,
zp=0,
ovlp=0,
plot=TRUE,
magt=10,
magf=10,
maga=2,
palette=rev.terrain.colors
)

{
n<-nrow(wave)
step<-seq(1,n-wl,wl-(ovlp*wl/100))		# FT windows

z1<-matrix(data=numeric((wl+(zp))*length(step)),wl+zp,length(step))
zpl<-zp/2
if(zpl==0)
{
for(i in step)
{z1[,which(step==i)]<-Mod(fft(wave[i:(wl+i-1),]*hanning.w(wl)))}
}

else
{
for(i in step)
{z1[,which(step==i)]<-Mod(fft(c(1:zpl,wave[i:(wl+i-1),],1:zpl)*hanning.w(wl+zp)))}
}		
					
# to keep only the relevant frequencies (half of the FT)
z2<-z1[1:((wl+zp)/2),]	
# to get only values between 0 and 1
z3<-z2/max(z2)		      
# replaces 0 values in spectra (that cannot be processed by the following log10())
z4<-ifelse(z3==0,yes=1e-6,no=z3) 
# to get dB values
z<-20*log10(z4)			      
	
if (plot == FALSE)
	return(z)
else
  {
	X<-magt*(1:ncol(z))
	Y<-magf*(1:nrow(z))
	Z<-maga*z
  
  # x,y,z ticks positions (for 5 ticks)
  Xat<-seq(magt,magt*ncol(z),by=(magt*ncol(z))/4)  
  Yat<-seq(magf,magf*nrow(z),by=(magf*nrow(z))/4)     
  Zat<-seq(min(Z),maga,by=abs(min(Z))/4)               
  
  # x,y,z, labels for 5 ticks
  Xlab<-as.character(round(seq(0,n/f,by=n/(4*f)),1))    
  Ylab<-as.character(round(seq((f/1000)/(wl+zp),f/2000,by=f/(4*2000)),1))
  Zlab<-as.character(round(seq(min(Z)/maga,0,by=abs(min(Z))/(4*maga)),1))
  
  Zlim <- range(Z)
  Zlen <- Zlim[2]-Zlim[1]+1
  colorlut <- palette(Zlen)
  col <- colorlut[Z-Zlim[1]+1] 
  
  library(rgl)
	rgl.clear()
  	rgl.bbox(color="white", emission="gray8", specular="gray",
    shininess=50, alpha=0.8,
    xat = Yat, xlab = Ylab, xunit=0, 
    yat = Zat, ylab = Zlab, yunit=0, 
    zat = Xat, zlab = Xlab, zunit=0) 
  rgl.surface(Y,X,Z,color=col, back = "line")
  text3d(x=1,z=magt*ncol(Z)/2,y=min(Z),text="Time (s)",color="white")
  text3d(z=1,x=magf*nrow(Z)/2,y=min(Z),text="Frequency (kHz)",color="white")
  text3d(x=1,z=0,y=0,text="Amplitude (dB)",color="white")
  }
}



################################################################################
#                                SYNTH                                         
################################################################################

synth<-function(
f, 
d, 
cf,  
a=1,  
p = 0,  
am=c(0,0),
fm=c(0,0,0), 
plot=TRUE,
wl = 512,
ovlp = 50,
...
)

{
amp<-am[1]/100  # AM modulation percentage
amf<-am[2]  # AM modulation frequency
fme<-fm[1]  # FM sinusoidal excursion
fmf<-fm[2]  # FM sinusoidal frequency
fmE<-fm[3]  # FM linear excursion

t <- seq(0, d*2*pi, length.out = f*d)


if (fme==0 & fmf!=0)          stop("FM sinusoidal excursion has to be set")
if (fme!=0 & fmf==0 & fmE==0) stop("FM sinusoidal frequency or FM linear excursion has to be set")
if (fme!=0 & fmf==0 & fmE!=0) stop("FM sinusoidal frequency has to be set")

if (fmE>0) freq<-seq(0,fmE/2,length.out=f*d) else freq<-rev(seq(fmE/2,0,length.out=f*d))

if (fme==0 & fmf==0) {sound0<-a*(1+amp*cos(amf*t))*sin((cf*t)+(freq*t)+p)}

if (fme!=0 & fmf!=0)
  {
  if (fmE == 0)       sound0<-a*(1+amp*cos(amf*t))*sin(cf*t+(fme/fmf)*sin(fmf*t+p)+p)
  else                sound0<-a*(1+amp*cos(amf*t))*sin(cf*t+(fme/fmf)*sin(fmf*t+p)+(freq*t)+p)
  }

sound1<-as.matrix(sound0)

if (plot == TRUE) spectro(sound1, f=f, wl=wl, ovlp=ovlp,...) else return(sound1)
}


################################################################################
#                                TIMER                                        
################################################################################

timer<-function(
wave,
f,
threshold,
smooth=0,
plot = TRUE,
plotthreshold = TRUE,
col = "black",
colval = "red",
xlab = "Time (s)",
ylab = "Amplitude",
...)


{
n<-nrow(wave)
thres<-max(abs(wave))*(threshold/100)

if (smooth)
  {
  z0<-abs(wave[0:n,])
  step<-seq(1,n-smooth,smooth)
  z1<-numeric(length(step))
  for(i in step)
    {z1[which(step==i)]<-mean(z0[i:(i+smooth)])}
  data<-as.matrix(z1)
  n<-nrow(data)
  thres<-max(data)*(threshold/100)
  f<-f/smooth
  wave1<-ts(z1,start=0,end=n/f,freq=f)}
else
  wave1<-ts(abs(wave[0:n,]),start=0,end=n/f,freq=f)

wave2<-ifelse(wave1<=thres,1,2)
wave3<-ts(wave2,start=0,end=n/f,freq=f)

# add successive values in wave1,
# values of 3 corresponds to the end of a silence or signal period
wave4<-numeric(n)  
for (i in 1:(n-1))  {wave4[i]<-wave2[i]+wave2[i+1]}

# sets a value of 3 at the first and last point of the signal
wave4[1]<-3
wave4[n]<-3

# gives the indeces of the end of a pause or signal period
wave5<-which(wave4==3)

# calculates the interval index between two successive ZC
nn<-length(wave5)
wave6<-numeric(nn-1)
for (i in 1:(nn-1)) {wave6[i]<-wave5[i+1]-wave5[i]}
# calculates signal and pause durations
y<-wave6/f

if (wave2[1]==1)  # if the file starts with a pause
  {
  pause<-y[seq(1,nn-1,by=2)]
  signal<-y[seq(2,nn-1,by=2)]
  }
  
else             # if the file starts with a signal
  {
  signal<-y[seq(1,nn-1,by=2)]
  pause<-y[seq(2,nn-1,by=2)]
  }
  
# computes the signal/pause ratio = duty cycle
ratio<-sum(signal)/sum(pause)

if (plot == TRUE)
  {
  plot(wave1/max(wave1),xlab=xlab,ylab=ylab,yaxt="n",ylim=c(0,1+0.1),col=col,...)
  if (plotthreshold == TRUE)
      {
      abline(h=thres, col=colval,lty=2)
      mtext(paste(as.character(threshold),"%"),
          side=2,line=0.5,at=thres,las=1,col=colval,cex=0.8)
      }
  par(new=TRUE)
  plot(wave3,xlab="",ylab="",yaxt="n",type="l",col=colval,ylim=c(1,2+0.1),...)

  wave8<-numeric(nn-1)
  for (i in 2:nn) {wave8[i]<-((wave5[i]-wave5[i-1])/2)+wave5[i-1]}
      
  if (wave2[1]==1)  # if the file starts with a pause
  {
  wave8.1<-wave8[seq(2,nn,by=2)]/f
  wave8.2<-wave8[seq(3,nn,by=2)]/f
  }
  else              # if the file starts with a signal
  {
  wave8.2<-wave8[seq(2,nn,by=2)]/f
  wave8.1<-wave8[seq(3,nn,by=2)]/f
  }
      
  ypl<-as.character(round(pause,2))
  ysl<-as.character(round(signal,2))
  text(x=wave8.1,y=1.075,ypl,col=colval,cex=0.8)
  text(x=wave8.2,y=2.075,ysl,col=colval,cex=0.8)
  
  }
else
  {
  timer<-list(s = signal,p = pause, r = ratio)
  return(timer)
  }
  
}


################################################################################
#                                ZC                                        
################################################################################

zc<-function(
wave,
f,
plot = TRUE,
interpol = 1,
threshold = FALSE,
xlab = "Time (s)",
ylab = "Frequency (kHz)",
ylim = c(0,f/2000),
...)

{
if (interpol > 5)
  {
  cat("please wait...")
  if (.Platform$OS.type == "windows") flush.console()
  }

n<-nrow(wave)

if (threshold)
  {
  thres<-max(abs(wave))*(threshold/100)
  wave1<-ifelse(abs(wave)<=thres,yes=0,no=1)
  wave2<-as.matrix(wave[,1]*wave1[,1])
  wave<-as.matrix(wave2)
  }
  
if (interpol > 1)
  {
  waveinterpol<-approx(wave,n=n*interpol)
  wave<-as.matrix(waveinterpol$y)
  F<-f*interpol
  }
else F<-f

# replaces null or positive values by 1 and negative values by 2  
wave1<-ifelse(wave>=0,1,2) 

# adds successive values in wave1, values of 3 corresponds to ZC
wave2<-numeric(n*interpol)  
for (i in 1:((n*interpol)-1))  {wave2[i]<-wave1[i]+wave1[i+1]}
# replaces 2 by 0
wave3<-ifelse(wave2==2, yes=0, no=wave2)
# replaces 4 by 0
wave4<-ifelse(wave3==4, yes=0, no=wave3)
# replaces 3 by their index
wave5<-replace(wave4,which(wave4==3),which(wave4==3))

# computes the period T between two successive zc 
wave6<-which(wave2==3) 
nn<-length(wave6) 
wave7<-numeric(nn)
for (i in 2:(nn-1)) {wave7[i]<-wave6[i+1]-wave6[i-1]}

# replaces indeces by T
wave8<-replace(wave5,which(wave5!=0),wave7)

# calculates the frequency
wave9<-F/(wave8)/1000 
y<-replace(wave9,which(wave9==Inf),NA)

if (plot == TRUE)
{
#if (args(type) == "l") stop("type "l" is not allowed")
x<-seq(0,n/f,length.out=n*interpol)
plot(x = x, y = y, xlab=xlab, ylab=ylab, las=1, ylim = ylim,...)
}
else return(y)
}



  
                       ###########################
                       ###########################
                       ### ACCESSORY FUNCTIONS ###
                       ###########################
                       ###########################

 
################################################################################
#                                FILLED.CONTOUR.MODIF1                                        
################################################################################ 
# modification of filled.contour in graphics by Ross Ihaka

filled.contour.modif1<-
function (x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)),
    z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE),
    zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels),
    nlevels = 20, color.palette = rev.gray.colors,
    col = color.palette(length(levels) - 1),
    plot.title, plot.axes, key.title, key.axes, scalelab, scalefontlab,
    asp = NA, xaxs = "i", yaxs = "i", axisX = TRUE, axisY = TRUE, axes = TRUE,
    las=1,...)
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq(0, 1, len = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
        stop("increasing 'x' and 'y' values expected")
    plot.new()
    par (mar=c(2,0.5,7,4)+0.1,las=1,cex=0.75)
        plot.window(xlim = c(0, 0.2), ylim = range(levels), xaxs = "i",
        yaxs = "i")
    mtext(text=scalelab,side=3,outer=FALSE,line=1.3,adj=0,cex=0.75,font=scalefontlab)
    rect(0, levels[-length(levels)], 0.2, levels[-1], col = col, lty=0)

    if (missing(key.axes)) {
        if (axes)
            axis(4,tcl=-0.5)
    }
    else key.axes
    box()
    if (!missing(key.title))
        key.title
    plot.new()
    par(mar=c(2,4,4,0)+0.1,cex=0.75)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
        stop("no proper 'z' matrix specified")
    if (!is.double(z))
        storage.mode(z) <- "double"
    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels),
        col = col))
    if (missing(plot.axes))
      {
        if(axisX)
            {
            title(main="", xlab="",ylab="")
            axis(1)
            }
        if(axisY)
            {
            title(main="", xlab="",ylab="")
            axis(2)
            }
      }
    box()
    if (missing(plot.title))
        title(...)
    else plot.title
    invisible()
}


################################################################################
#                                FILLED.CONTOUR.MODIF2                                        
################################################################################ 
# modification of filled.contour in graphics by Ross Ihaka

filled.contour.modif2<-function (x = seq(0, 1, len = nrow(z)),
    y = seq(0, 1, len = ncol(z)), z, xlim = range(x, finite = TRUE),
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
    col = color.palette(length(levels) - 1), plot.title, plot.axes, key.title,
    asp = NA, xaxs = "i", yaxs = "i", las = 1, axisX = TRUE, axisY = TRUE, 
    ...) 
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq(0, 1, len = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
        stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
        storage.mode(z) <- "double"
    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
        col = col))
    if (missing(plot.axes))
      {
        if(axisX)
            {
            title(main="", xlab="",ylab="")
            axis(1)
            }
        if(axisY)
            {
            title(main="", xlab="",ylab="")
            axis(2)
            }
      }
    box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
}


################################################################################
#                                FILLED.CONTOUR.MODIF3                                        
################################################################################ 
# modification of filled.contour in graphics by Ross Ihaka

filled.contour.modif3<-
function (x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)),
    z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE),
    zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels),
    nlevels = 20, color.palette = rev.gray.colors,
    col = color.palette(length(levels) - 1), plot.title, plot.axes, key.title,
    key.axes, scalelab, scalefontlab, asp = NA,
    xaxs = "i", yaxs = "i", axisX = TRUE, axisY = TRUE, axes = TRUE, las=1,
    ...)
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq(0, 1, len = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
        stop("increasing 'x' and 'y' values expected")
    plot.new()
    par (mar=c(3.8,0.5,6,4)+0.1,las=1,cex=0.75)
        plot.window(xlim = c(0, 0.2), ylim = range(levels), xaxs = "i",
        yaxs = "i")
    mtext(text=scalelab,side=3,outer=FALSE,line=1.6,adj=0,cex=0.75,font=scalefontlab)
    rect(0, levels[-length(levels)], 0.2, levels[-1], col = col, lty=0)

    if (missing(key.axes)) {
        if (axes)
            axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title))
        key.title
    plot.new()
    par(mar=c(5,4,4,0)+0.1,cex=0.75)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
        stop("no proper 'z' matrix specified")
    if (!is.double(z))
        storage.mode(z) <- "double"
    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels),
        col = col))
    if (missing(plot.axes))
      {
        if(axisX)
            {
            title(main="", xlab="",ylab="")
            axis(1)
            }
        if(axisY)
            {
            title(main="", xlab="",ylab="")
            axis(2)
            }
      }
    box()
    if (missing(plot.title))
        title(...)
    else plot.title
    invisible()
}


################################################################################
#                                HAMMING.W                                        
################################################################################ 
# hamming.window from Andreas Weingessel (library e1071)

hamming.w<-function (n)
{
if (n == 1)
  c <- 1
else {
  n <- n - 1
  c <- 0.54 - 0.46 * cos(2 * pi * (0:n)/n)
    }
return(c)
}


################################################################################
#                                HANNING.W                                        
################################################################################ 
# hanning.window from Andreas Weingessel (library e1071)

hanning.w<-function (n)
{
if (n == 1)
  c <- 1
else
  {
  n <- n - 1
  c <- 0.5 - 0.5 * cos(2 * pi * (0:n)/n)
  }
return(c)
}


################################################################################
#                         PEAKS, PEAKSIGN, CHECK.PCKS                                     
################################################################################
# Author: Martin Maechler, Date: 25 Nov 2005
# Martin Maechler <maechler@stat.math.ethz.ch>
# Peaksign: return (-1 / 0 / 1) if series[i] is ( trough / "normal" / peak )

peaks <- function(series, span = 3, do.pad = TRUE) {
    if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
    s1 <- 1:1 + (s <- span %/% 2)
    if(span == 1) return(rep.int(TRUE, length(series)))
    z <- embed(series, span)
    v <- apply(z[,s1] > z[, -s1, drop=FALSE], 1, all)
    if(do.pad) {
        pad <- rep.int(FALSE, s)
        c(pad, v, pad)
    } else v
}

peaksign <- function(series, span = 3, do.pad = TRUE)
{
    if((span <- as.integer(span)) %% 2 != 1 || span == 1)
        stop("'span' must be odd and >= 3")
    s1 <- 1:1 + (s <- span %/% 2)
    z <- embed(series, span)
    d <- z[,s1] - z[, -s1, drop=FALSE]
    ans <- rep.int(0:0, nrow(d))
    ans[apply(d > 0, 1, all)] <- as.integer(1)
    ans[apply(d < 0, 1, all)] <- as.integer(-1)
    if(do.pad) {
        pad <- rep.int(0:0, s)
        c(pad, ans, pad)
    } else ans
}

check.pks <- function(y, span = 3)
    stopifnot(identical(peaks( y, span), peaksign(y, span) ==  1),
              identical(peaks(-y, span), peaksign(y, span) == -1))


################################################################################
#                                REV.CM.COLORS                                       
################################################################################
# rev.cm.colors, reversion of cm.colors in grDevices package
# originally by R Development Core Team and contributors worldwide

rev.cm.colors<-
function (x)
{
    n<-x
    if ((n <- as.integer(n[1])) > 0) {
        even.n <- n%%2 == 0
        k <- n%/%2
        l1 <- k + 1 - even.n
        l2 <- n - k + even.n
        rev(c(if (l1 > 0) hsv(h = 6/12, s = seq(0.5, ifelse(even.n, 
            0.5/k, 0), length = l1), v = 1), if (l2 > 1) hsv(h = 10/12, 
            s = seq(0, 0.5, length = l2)[-1], v = 1)))
    }
    else character(0)
}



################################################################################
#                                REV.GRAY.COLORS                                       
################################################################################ 
rev.gray.colors<-
function (x, start=1, end=0, gamma = 1.7)
gray(seq(from = start^gamma, to = end^gamma, length = x)^(1/gamma))



################################################################################
#                                REV.HEAT.COLORS                                       
################################################################################
# rev.heat.colors, reversion of heat.colors in grDevices package
# originally by R Development Core Team and contributors worldwide 

rev.heat.colors<-
function (x) 
{
    n<-x
    if ((n <- as.integer(n[1])) > 0) {
        j <- n%/%4
        i <- n - j
        rev(c(rainbow(i, start = 0, end = 1/6), if (j > 0) hsv(h = 1/6, 
            s = seq(from = 1 - 1/(2 * j), to = 1/(2 * j), length = j), 
            v = 1)))
    }
    else character(0)
}



################################################################################
#                                REV.TERRAIN.COLORS                                       
################################################################################
# rev.terrain.colors, reversion of terrain.colors in grDevices package
# originally by R Development Core Team and contributors worldwide 

rev.terrain.colors<-
function (x)
{
    n<-x
    if ((n <- as.integer(n[1])) > 0) {
        k <- n%/%2
        h <- c(4/12, 2/12, 0/12)
        s <- c(1, 1, 0)
        v <- c(0.65, 0.9, 0.95)
        rev(c(
        hsv(h = seq(h[1], h[2], length = k),
        s = seq(s[1], s[2], length = k),
        v = seq(v[1], v[2], length = k)),
        
        hsv(h = seq(h[2], h[3], length = n - k + 1)[-1],
        s = seq(s[2], s[3], length = n - k + 1)[-1],
        v = seq(v[2], v[3], length = n - k + 1)[-1])
        ))
    }
    else character(0)
}



################################################################################
#                                REV.TOPO.COLORS                                       
################################################################################
# rev.topo.colors, reversion of topo.colors in grDevices package
# originally by R Development Core Team and contributors worldwide 

rev.topo.colors<-
function (x) 
{
    n<-x
    if ((n <- as.integer(n[1])) > 0) {
        j <- n%/%3
        k <- n%/%3
        i <- n - j - k
        rev(c(if (i > 0) hsv(h = seq(from = 43/60, to = 31/60, length = i)), 
            if (j > 0) hsv(h = seq(from = 23/60, to = 11/60, 
                length = j)), if (k > 0) hsv(h = seq(from = 10/60, 
                to = 6/60, length = k), s = seq(from = 1, to = 0.3, 
                length = k), v = 1)))
    }
    else character(0)
}



################################################################################
#                                SOSCILLO                                        
################################################################################

soscillo<-function
(
wave,
f,
env = FALSE,
smooth = 0,
xlab = "Times (s)",
ylab = "Amplitude",
colwave = "black",
colbg = "white",
coltitle = "black",
collab = "black",
colline = "black",
colaxis = "black",
coly0 = "grey47",
title = FALSE,
... 
)

{
n<-nrow(wave)
par(tcl=0.5, col.axis=colaxis, col=colline, col.lab=collab,las=0)

if (env == TRUE)
  {
    if (smooth)
    {
    z0<-abs(wave[0:n,])
    step<-seq(1,n-smooth,smooth)		# smooth windows
    z1<-numeric(length(step)) 
    for(i in step)
      {z1[which(step==i)]<-mean(z0[i:(i+smooth)])}
    data<-as.matrix(z1)
    n<-nrow(data)
    f<-f/smooth
    wave<-ts(z1,start=0,end=n/f,freq=f)
    }
    else
    wave<-ts(abs(wave[0:n,]),start=0,end=n/f,freq=f)
  }

else
  wave<-ts(wave[0:n,], start=0, end=n/f, freq=f)

plot(wave,
		col=colwave, type="l",
		xaxs="i", yaxs="i",
		xlab=xlab, ylab="",
		yaxt="n", bty="l",
		...)
axis(side=1, col=colline,labels=FALSE)
axis(side=2, at=max(wave,na.rm=TRUE), col=colline,labels=FALSE)
mtext(text=ylab,side=2,cex=0.85,line=2.75)
abline(h=0,col=coly0,lty=2)
}



################################################################################
#                                SSPECTRO                                        
################################################################################

sspectro <- function
(
wave,
f,
wl
)

{
n<-nrow(wave)
step<-seq(1,n-wl,wl)		# FT windows

z1<-matrix(data=numeric(wl*length(step)),wl,length(step))

for(i in step)
  {
  z1[,which(step==i)]<- Mod(fft(wave[i:(wl+i-1),]*hanning.w(wl)))
  }

z2<-z1[1:(wl/2),]
z3<-z2/max(z2)
return(z3)
}



################################################################################
#                                SPECTRO.COLORS                                        
################################################################################

spectro.colors<-
function (n)
{
if ((n <- as.integer(n[1])) > 0)
 {
 j <- n%/%3
 k <- n%/%3
 i <- n - j - k
 c(if (i > 0) hsv(h = seq(from = 31/60, to = 43/60, length = i), s = seq(0,1,length=i)),
   if (j > 0) hsv(h = seq(from = 21/60, to = 9/60, length = j), v = seq(0.5,0.8,length=j)),
   if (k > 0) hsv(h = seq(from = 8/60, to = 1/60, length = k), s = seq(from = 0.5, to = 1, length = k), v=1))
 }
else character(0)
}

