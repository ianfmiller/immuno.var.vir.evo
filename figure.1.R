### immunity distribution shapes, means, and variances
library(deSolve)
source("~/Documents/GitHub/immuno.var.vir.evo/cluster/code/analysis.setup.R")

plot.colors<-c(colorRampPalette(c("green4","palegreen"))(8),"goldenrod",colorRampPalette(c("lightskyblue","royalblue2"))(8))


alphas<-10^(seq(-2,2,.25))
x<-seq(0,1,.01)

layout.mat<-matrix(NA,6,10)
layout.mat[1,]<-rep(c(1,2,3,4,5),each=2)
layout.mat[2,]<-rep(c(6,7,8,9,10),each=2)
layout.mat[3,]<-rep(c(11,12,13,14,15),each=2)
layout.mat[4,]<-c(22,16,16,17,17,18,18,19,19,23)
layout.mat[5,]<-rep(c(20,21),each=5)
layout.mat[6,]<-rep(c(20,21),each=5)
layout(layout.mat)
par(mar=c(4,4,3,1))

set.immunity.dist.split()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
plot.startconds.2(ymax=.5,col="black",main="polarized")
mtext("A",3,adj=1,font=2)

for (alpha in alphas)
{
  beta<-alpha  
  set.immunity.dist.beta(alpha,beta)
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  plot.startconds.2(ymax=.5,col=plot.colors[which(alphas==alpha)],main=substitute(paste(theta," = ",10^z),list(z=log10(alpha))))
  mtext(c("B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R")[which(alphas==alpha)],3,adj=1,font=2)
}

set.immunity.dist.single()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
plot.startconds.2(ymax=1,col="darkgrey",main="homogeneous")
mtext("S",3,adj=1,font=2)


par(mar=c(7,5,3,5))
plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,1),xlab=expression(log[10](alpha)),ylab="mean",type="n",axes = F,cex.lab=2)
mtext("T",3,adj=1,font=2)
box()
axis(1,at=c(-2,-1,0,1,2),labels=c(-2,-1,0,1,2),cex.axis=1.5)
axis(1,at=c(-2.25,2.25),labels=c("polarized","homog."),las=2,cex.axis=1.5)
axis(2,cex.axis=1.5)
abline(v=-2.125,col="grey",lty=2,lwd=2.5)
abline(v=2.125,col="grey",lty=2,lwd=2.5)

set.immunity.dist.split()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
yy<-weighted.mean(immunity.categories,scaled.immunity.distribution(dummy.immunity.categories))
points(-2.25,yy,col="black",pch=16,cex=4)

for (alpha in alphas)
{
  beta<-alpha  
  set.immunity.dist.beta(alpha,beta)
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  yy<-weighted.mean(immunity.categories,scaled.immunity.distribution(dummy.immunity.categories))
  points(log10(alpha),yy,col=plot.colors[which(alphas==alpha)],pch=16,cex=4)
}

set.immunity.dist.single()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
yy<-weighted.mean(immunity.categories,scaled.immunity.distribution(dummy.immunity.categories))
points(2.25,yy,col="darkgrey",pch=16,cex=4)

par(mar=c(7,5,3,5))
plot(0,0,xlim=c(-2.25,2.25),ylim=c(-.05,.3),xlab=expression(log[10](alpha)),ylab="variance",type="n",axes = F,cex.lab=2)
mtext("U",3,adj=1,font=2)
box()
axis(1,at=c(-2,-1,0,1,2),labels=c(-2,-1,0,1,2),cex.axis=1.5)
axis(1,at=c(-2.25,2.25),labels=c("polarized","homog."),las=2,cex.axis=1.5)
axis(2,cex.axis=1.5)
abline(v=-2.125,col="grey",lty=2,lwd=2.5)
abline(v=2.125,col="grey",lty=2,lwd=2.5)

set.immunity.dist.split()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
yy<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
points(-2.25,yy,col="black",pch=16,cex=4)

for (alpha in alphas)
{
  beta<-alpha  
  set.immunity.dist.beta(alpha,beta)
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  yy<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
  points(log10(alpha),yy,col=plot.colors[which(alphas==alpha)],pch=16,cex=4)
}

set.immunity.dist.single()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
yy<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
points(2.25,yy,col="darkgrey",pch=16,cex=4)


