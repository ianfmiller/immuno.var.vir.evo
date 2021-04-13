source("~/Documents/GitHub/immuno.var.vir.evo/cluster/code/analysis.setup.R")

colors<-c(colorRampPalette(c("green4","palegreen"))(8),"goldenrod",colorRampPalette(c("lightskyblue","royalblue2"))(8))

point.func<-function(file.path,pch=16)
{
  data<-read.csv(file.path)
  es.virs<-data$ES.vir
  special.outcomes<-data$hyper.vir.

  
  set.immunity.dist.split()
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
  points(y,es.virs[which(data$shape==-1)],col="black",cex=1.5,pch=pch)
  
  for (alpha in alphas)
  {
    beta<-alpha  
    set.immunity.dist.beta(alpha,beta)
    initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
    get.startconds()
    y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
    points(y,es.virs[which.min(abs(data$shape-alpha))],col=colors[which(alphas==alpha)],pch=pch,cex=1.5)
  }
  
  set.immunity.dist.single()
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
  points(y,es.virs[which(data$shape==-2)],col="darkgrey",pch=pch,cex=1.5)
  
  for(index in which(special.outcomes=="selection for hypervirulence"))
  {
    if (data[index,"shape"]==-1)
    {
      set.immunity.dist.split()
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      points(y,.475,col="black",pch=pch,cex=1.5)
    }
    
    if (!data[index,"shape"] %in% c(-1,-2))
    {
      alpha<-alphas[which.min(abs(data[index,"shape"]-alphas))]
      beta<-alpha  
      set.immunity.dist.beta(alpha,beta)
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      points(y,.475,col=colors[which(alphas==alpha)],pch=pch,cex=1.5)
    }
    
    if (data[index,"shape"]==-2)
    {
      set.immunity.dist.single()
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      points(y,.475,col="grey",pch=pch,cex=1.5)
    }
  }
  
  for(index in which(special.outcomes=="global eradication"))
  {
    if (data[index,"shape"]==-1)
    {
      set.immunity.dist.split()
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      points(y,.025,col="black",pch=pch,cex=1.5)
    }
    
    if (!data[index,"shape"] %in% c(-1,-2))
    {
      alpha<-alphas[which.min(abs(data[index,"shape"]-alphas))]
      beta<-alpha  
      set.immunity.dist.beta(alpha,beta)
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      points(y,.025,col=colors[which(alphas==alpha)],pch=pch,cex=1.5)
    }
    
    if (data[index,"shape"]==-2)
    {
      set.immunity.dist.single()
      initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
      get.startconds()
      y<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
      points(y,.025,col="grey",pch=pch,cex=1.5)
    }
  }
}

par(mfrow=c(2,3))
par(mar=c(5.1,5.1,2.1,2.1))
#layout.mat<-matrix(c(1,1,2,2,3,3,6,4,4,5,5,7),2,6,byrow = T)
#layout(layout.mat)

plot(0,0,xlim=c(0,.25),ylim=c(0,.5),xlab="variance",ylab="ES virulence",type="n",axes=F)
mtext("W=1 X=0 Y=0 Z=0")
axis(1)
axis(2,at=seq(.1,.4,.1))
axis(2,at=c(.025,.475),labels=c("erad.","hypervir."),las=2,tick = F)
abline(h=.05,lty=2,col="grey",lwd=2)
abline(h=0,lty=2,col="grey",lwd=2)
abline(h=.45,lty=2,col="grey",lwd=2)
abline(h=.5,lty=2,col="grey",lwd=2)
point.func(file.path="~/Documents/GitHub/immuno.var.vir.evo/results/p1W0X0Y0Z1.csv",pch=16)


point.func("p0.9W1X0Y0Z0.csv",p.vacc=.9,pch=17)
point.func("p0.5W1X0Y0Z0.csv",p.vacc=.5,pch=15)






plot(0,0,xlim=c(0,.25),ylim=c(-.04,.24),xlab="variance",ylab="ES virulence",type="n",axes=F)
mtext("W=0 X=1 Y=0 Z=1")
axis(1)
axis(2,at=c(0,.1,.2))
axis(2,at=c(-.025,.225),labels=c("erad.","hypervir."),las=2,tick = F)
abline(h=-.04,lty=2,col="grey")
abline(h=-.01,lty=2,col="grey")
abline(h=.21,lty=2,col="grey")
abline(h=.24,lty=2,col="grey")
point.func("p1W0X1Y0Z1.csv",p.vacc=1,pch=16)
point.func("p0.9W0X1Y0Z1.csv",p.vacc=.9,pch=17)
point.func("p0.5W0X1Y0Z1.csv",p.vacc=.5,pch=15)

plot(0,0,xlim=c(0,.25),ylim=c(-.04,.24),xlab="variance",ylab="ES virulence",type="n",axes=F)
mtext("W=0 X=0 Y=1 Z=0")
axis(1)
axis(2,at=c(0,.1,.2))
axis(2,at=c(-.025,.225),labels=c("erad.","hypervir."),las=2,tick = F)
abline(h=-.04,lty=2,col="grey")
abline(h=-.01,lty=2,col="grey")
abline(h=.21,lty=2,col="grey")
abline(h=.24,lty=2,col="grey")
point.func("p1W0X0Y1Z0.csv",p.vacc=1,pch=16)
point.func("p0.9W0X0Y1Z0.csv",p.vacc=.9,pch=17)
point.func("p0.5W0X0Y1Z0.csv",p.vacc=.5,pch=15)


plot(0,0,xlim=c(0,.25),ylim=c(-.1,.6),xlab="variance",ylab="ES virulence",type="n",axes=F)
mtext("W=0 X=0 Y=0 Z=1")
axis(1)
axis(2,at=c(0,.1,.2,.3,.4,.5))
axis(2,at=c(-.0625,.5625),labels=c("erad.","hypervir."),las=2,tick = F)
abline(h=-.1,lty=2,col="grey")
abline(h=-.025,lty=2,col="grey")
abline(h=.525,lty=2,col="grey")
abline(h=.6,lty=2,col="grey")
point.func("p1W0X0Y0Z1.csv",p.vacc=1,pch=16)
point.func("p0.9W0X0Y0Z1.csv",p.vacc=.9,pch=17)
point.func("p0.5W0X0Y0Z1.csv",p.vacc=.5,pch=15)

plot(0,0,xlim=c(0,.25),ylim=c(-.1,.6),xlab="variance",ylab="ES virulence",type="n",axes=F)
mtext("W=0 X=0.5 Y=0 Z=1")
axis(1)
axis(2,at=c(0,.1,.2,.3,.4,.5))
axis(2,at=c(-.0625,.5625),labels=c("erad.","hypervir."),las=2,tick = F)
abline(h=-.1,lty=2,col="grey")
abline(h=-.025,lty=2,col="grey")
abline(h=.525,lty=2,col="grey")
abline(h=.6,lty=2,col="grey")
point.func("p1W0X0.5Y0Z1.csv",p.vacc=1,pch=16)
point.func("p0.9W0X0.5Y0Z1.csv",p.vacc=.9,pch=17)
point.func("p0.5W0X0.5Y0Z1.csv",p.vacc=.5,pch=15)

plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="")
points(0,.75,pch=16,cex=2)
points(0,.5,pch=15,cex=2)
points(0,.25,pch=17,cex=2)
text(.025,.75,"%100 vaccination",pos=4)
text(.025,.5,"%90 vaccination",pos=4)
text(.025,.25,"%50 vaccination",pos=4)

for (i in 1:19)
{
  points(.6,seq(0,1,length.out = 19)[i],col=c("black",colors,"grey")[i],pch=16,cex=1.5)
  if (i==19) {text(.65,seq(0,1,length.out = 19)[i],"fixed",pos=4)}
  if (i %in% c(6,10,14)) {{text(.65,seq(0,1,length.out = 19)[i],substitute(paste(alpha," = ",10^z),list(z=log10(alphas[i-1]))),pos=4)}}
  if (i==1) {text(.65,seq(0,1,length.out = 19)[i],"split",pos=4)}
}





