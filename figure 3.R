source("~/Documents/GitHub/immuno.var.vir.evo/cluster/code/analysis.setup.R")

colors<-c(colorRampPalette(c("green4","palegreen"))(8),"goldenrod",colorRampPalette(c("lightskyblue","royalblue2"))(8))

point.func<-function(file.path,pch=16,cex=2,wval,xval,yval,zval)
{
  
  w<<-function(sus.immunity.class) {wval*sus.immunity.class}
  x<<-function(sus.immunity.class) {xval*sus.immunity.class}
  y<<-function(inf.immunity.class) {yval*inf.immunity.class}
  z<<-function(sus.immunity.class) {zval*sus.immunity.class}
  
  p.immune<<-1
  birth.rate<<-.07
  death.rate<<-.05
  K<-1/(-death.rate/birth.rate+1)
  b1<<-3
  b2<<-.333
  d1<<-.1
  d2<<-0
  c1<<-1
  set.nonlinear.tradeoff(b1,b2,d1,d2,c1)
  
  data<-read.csv(file.path)
  es.virs<-data$ES.vir
  special.outcomes<-data$hyper.vir.
  
  set.immunity.dist.split()
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  if(is.na(special.outcomes[which(data$shape==-1)])) {v1<<-as.numeric(es.virs[which(data$shape==-1)])}
  if(isTRUE(special.outcomes[which(data$shape==-1)]=="selection for hypervirulence")) {v1<<-100}
  getCFRs()
  y<--2.25
  points(y,CFRs[22],col="black",cex=cex,pch=pch)
  
  for (alpha in alphas)
  {
    beta<-alpha  
    set.immunity.dist.beta(alpha,beta)
    initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
    get.startconds()
    if(is.na(special.outcomes[which.min(abs(data$shape-alpha))])) {v1<<-as.numeric(es.virs[which(data$shape==-2)])}
    if(isTRUE(special.outcomes[which.min(abs(data$shape-alpha))]=="selection for hypervirulence")) {v1<<-100}
    getCFRs()
    y<-log10(alpha)
    points(y,CFRs[22],col=colors[which(alphas==alpha)],cex=cex,pch=pch)
  }
  
  set.immunity.dist.single()
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  if(is.na(special.outcomes[which(data$shape==-2)])) {v1<<-as.numeric(es.virs[which(data$shape==-2)])}
  if(isTRUE(special.outcomes[which(data$shape==-2)]=="selection for hypervirulence")) {v1<<-100}
  getCFRs()
  y<-2.25
  points(y,CFRs[22],col="grey",cex=cex,pch=pch)
}

#layout(matrix(c(1,1,2,2,3,3,6,4,4,5,5,7),2,6,byrow = T))
par(mfrow=c(2,3))
par(mar=c(6.1,7.1,5.1,5.1))
#layout.mat<-matrix(c(1,1,2,2,3,3,6,4,4,5,5,7),2,6,byrow = T)
#layout(layout.mat)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(.1,.8),xlab=expression(log[10](alpha)),ylab="CFR",type="n",axes=F,cex.lab=2)
mtext("W=1 X=0 Y=0 Z=0",cex=2)
axis(1,cex.axis=1.75)
axis(2,at=seq(.1,.8,.1),cex.axis=1.75)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W1X0Y0Z0.csv",pch=16,cex=3.5,wval = 1, xval = 0, yval = 0, zval = 0)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(.1,.8),xlab=expression(log[10](alpha)),ylab="CFR",type="n",axes=F,cex.lab=2)
mtext("W=0 X=0 Y=1 Z=0",cex=2)
axis(1,cex.axis=1.75)
axis(2,at=seq(.1,.8,.1),cex.axis=1.75)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X0Y1Z0.csv",pch=16,cex=3.5,wval = 0, xval = 0, yval = 1, zval = 0)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(.1,.8),xlab=expression(log[10](alpha)),ylab="CFR",type="n",axes=F,cex.lab=2)
mtext("W=0 X=0.5 Y=0 Z=1",cex=2)
axis(1,cex.axis=1.75)
axis(2,at=seq(.1,.8,.1),cex.axis=1.75)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X0.5Y0Z1.csv",pch=16,cex=3.5,wval = 0, xval = 0.5, yval = 0, zval = 1)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(.1,.8),xlab=expression(log[10](alpha)),ylab="CFR",type="n",axes=F,cex.lab=2)
mtext("W=0 X=0 Y=0 Z=1",cex=2)
axis(1,cex.axis=1.75)
axis(2,at=seq(.1,.8,.1),cex.axis=1.75)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X0Y0Z1.csv",pch=16,cex=3.5,wval = 0, xval = 0, yval = 0, zval = 1)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(.1,.8),xlab=expression(log[10](alpha)),ylab="CFR",type="n",axes=F,cex.lab=2)
mtext("W=0 X=1 Y=0 Z=1",cex=2)
axis(1,cex.axis=1.75)
axis(2,at=seq(.1,.8,.1),cex.axis=1.75)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X1Y0Z1.csv",pch=16,cex=3.5,wval = 0, xval = 1, yval = 0, zval = 1)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(.1,.8),xlab=expression(log[10](alpha)),ylab="CFR",type="n",axes=F,cex.lab=2)
mtext("W=1 X=1 Y=1 Z=1",cex=2)
axis(1,cex.axis=1.75)
axis(2,at=seq(.1,.8,.1),cex.axis=1.75)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W1X1Y1Z1.csv",pch=16,cex=3.5,wval = 0, xval = 1, yval = 0, zval = 1)


