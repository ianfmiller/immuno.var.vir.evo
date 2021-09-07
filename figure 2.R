source("~/Documents/GitHub/immuno.var.vir.evo/cluster/code/analysis.setup.R")

colors<-c(colorRampPalette(c("green4","palegreen"))(8),"goldenrod",colorRampPalette(c("lightskyblue","royalblue2"))(8))

point.func<-function(file.path,pch=16,cex=2)
{
  data<-read.csv(file.path)
  es.virs<-data$ES.vir
  special.outcomes<-data$hyper.vir.

  
  set.immunity.dist.split()
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  #y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
  y<--2.25
  points(y,es.virs[which(data$shape==-1)],col="black",cex=cex,pch=pch)
  
  for (alpha in alphas)
  {
    beta<-alpha  
    set.immunity.dist.beta(alpha,beta)
    initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
    get.startconds()
    #y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
    y<-log10(alpha)
    points(y,es.virs[which.min(abs(data$shape-alpha))],col=colors[which(alphas==alpha)],pch=pch,cex=cex)
  }
  
  set.immunity.dist.single()
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  #y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
  y<-2.25
  points(y,es.virs[which(data$shape==-2)],col="darkgrey",pch=pch,cex=cex)
  
  for(index in which(special.outcomes=="selection for hypervirulence"))
  {
    if (data[index,"shape"]==-1)
    {
      set.immunity.dist.split()
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      #y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      y<--2.25
      points(y,.475,col="black",pch=pch,cex=cex)
    }
    
    if (!data[index,"shape"] %in% c(-1,-2))
    {
      alpha<-alphas[which.min(abs(data[index,"shape"]-alphas))]
      beta<-alpha  
      set.immunity.dist.beta(alpha,beta)
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      #y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      y<-log10(alpha)
      points(y,.475,col=colors[which(alphas==alpha)],pch=pch,cex=cex)
    }
    
    if (data[index,"shape"]==-2)
    {
      set.immunity.dist.single()
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      #y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      y<-2.25
      points(y,.475,col="grey",pch=pch,cex=cex)
    }
  }
}

#layout(matrix(c(1,1,2,2,3,3,6,4,4,5,5,7),2,6,byrow = T))
par(mfrow=c(2,3))
par(mar=c(6.1,7.1,5.1,5.1))
#layout.mat<-matrix(c(1,1,2,2,3,3,6,4,4,5,5,7),2,6,byrow = T)
#layout(layout.mat)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,.5),xlab=expression(log[10](alpha)),ylab="ES virulence",type="n",axes=F,cex.lab=2)
mtext("W=1 X=0 Y=0 Z=0",cex=2)
mtext("A",adj=-.25,cex=1.5)
axis(1,cex.axis=1.75)
axis(2,at=seq(0,.4,.1),cex.axis=1.75)
axis(2,at=.475,labels="hypervir.",las=2,tick = F,cex.axis=1.75)
abline(h=.45,lty=2,col="grey",lwd=2)
abline(h=.5,lty=2,col="grey",lwd=2)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W1X0Y0Z0.csv",pch=16,cex=3.5)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,.5),xlab=expression(log[10](alpha)),ylab="ES virulence",type="n",axes=F,cex.lab=2)
mtext("W=0 X=0 Y=1 Z=0",cex=2)
mtext("B",adj=-.25,cex=1.5)
axis(1,cex.axis=1.75)
axis(2,at=seq(0,.4,.1),cex.axis=1.75)
axis(2,at=.475,labels="hypervir.",las=2,tick = F,cex.axis=1.75)
abline(h=.45,lty=2,col="grey",lwd=2)
abline(h=.5,lty=2,col="grey",lwd=2)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X0Y1Z0.csv",pch=16,cex=3.5)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,.5),xlab=expression(log[10](alpha)),ylab="ES virulence",type="n",axes=F,cex.lab=2)
mtext("W=0 X=0.5 Y=0 Z=1",cex=2)
mtext("C",adj=-.25,cex=1.5)
axis(1,cex.axis=1.75)
axis(2,at=seq(0,.4,.1),cex.axis=1.75)
axis(2,at=.475,labels="hypervir.",las=2,tick = F,cex.axis=1.75)
abline(h=.45,lty=2,col="grey",lwd=2)
abline(h=.5,lty=2,col="grey",lwd=2)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X0.5Y0Z1.csv",pch=16,cex=3.5)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,.5),xlab=expression(log[10](alpha)),ylab="ES virulence",type="n",axes=F,cex.lab=2)
mtext("W=0 X=0 Y=0 Z=1",cex=2)
mtext("D",adj=-.25,cex=1.5)
axis(1,cex.axis=1.75)
axis(2,at=seq(0,.4,.1),cex.axis=1.75)
axis(2,at=.475,labels="hypervir.",las=2,tick = F,cex.axis=1.75)
abline(h=.45,lty=2,col="grey",lwd=2)
abline(h=.5,lty=2,col="grey",lwd=2)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X0Y0Z1.csv",pch=16,cex=3.5)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,.5),xlab=expression(log[10](alpha)),ylab="ES virulence",type="n",axes=F,cex.lab=2)
mtext("W=0 X=1 Y=0 Z=1",cex=2)
mtext("E",adj=-.25,cex=1.5)
axis(1,cex.axis=1.75)
axis(2,at=seq(0,.4,.1),cex.axis=1.75)
axis(2,at=.475,labels="hypervir.",las=2,tick = F,cex.axis=1.75)
abline(h=.45,lty=2,col="grey",lwd=2)
abline(h=.5,lty=2,col="grey",lwd=2)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X1Y0Z1.csv",pch=16,cex=3.5)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,.5),xlab=expression(log[10](alpha)),ylab="ES virulence",type="n",axes=F,cex.lab=2)
mtext("W=1 X=1 Y=1 Z=1",cex=2)
mtext("F",adj=-.25,cex=1.5)
axis(1,cex.axis=1.75)
axis(2,at=seq(0,.4,.1),cex.axis=1.75)
axis(2,at=.475,labels="hypervir.",las=2,tick = F,cex.axis=1.75)
abline(h=.45,lty=2,col="grey",lwd=2)
abline(h=.5,lty=2,col="grey",lwd=2)
axis(1,at=c(-2.25,2.25),labels=c("split","homog."),las=2,cex.axis=1.5,lwd=0,lwd.ticks = 1)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W1X1Y1Z1.csv",pch=16,cex=3.5)






