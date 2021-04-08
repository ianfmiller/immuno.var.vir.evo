### immunity distribution shapes, means, and variances
library(deSolve)
source("~/Documents/GitHub/immuno.var.vir.evo/analysis.setup.R")

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
par(mar=c(4,4,2,1))

set.immunity.dist.split()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
plot.startconds.2(ymax=.25,col="black",main="split")

for (alpha in alphas)
{
  beta<-alpha  
  set.immunity.dist.beta(alpha,beta)
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  z<-alpha
  plot.startconds.2(ymax=.25,col=plot.colors[which(alphas==alpha)],main=substitute(paste(alpha," = ",beta," = ",10^z)))
}

set.immunity.dist.single()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
plot.startconds.2(ymax=.25,col="darkgrey")

plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,1),xlab="log 10 alpha",ylab="mean",type="n",axes = F)
box()
axis(1,at=c(-2,-1,0,1,2),labels=c(-2,-1,0,1,2))
axis(1,at=c(-2.25,2.25),labels=c("split","fixed"),las=2)
axis(2)
abline(v=-2.125,col="grey",lty=2)
abline(v=2.125,col="grey",lty=2)

set.immunity.dist.split()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
yy<-weighted.mean(immunity.categories,scaled.immunity.distribution(dummy.immunity.categories))
points(-2.25,yy,col="black",pch=16,cex=1.5)

for (alpha in alphas)
{
  beta<-alpha  
  set.immunity.dist.beta(alpha,beta)
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  yy<-weighted.mean(immunity.categories,scaled.immunity.distribution(dummy.immunity.categories))
  points(log10(alpha),yy,col=plot.colors[which(alphas==alpha)],pch=16,cex=1.5)
}

set.immunity.dist.single()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
yy<-weighted.mean(immunity.categories,scaled.immunity.distribution(dummy.immunity.categories))
points(2.25,yy,col="darkgrey",pch=16,cex=1.5)

plot(0,0,xlim=c(-2.25,2.25),ylim=c(0,.25),xlab="log 10 alpha",ylab="var",type="n",axes = F)
box()
axis(1,at=c(-2,-1,0,1,2),labels=c(-2,-1,0,1,2))
axis(1,at=c(-2.25,2.25),labels=c("split","fixed"),las=2)
axis(2)
abline(v=-2.125,col="grey",lty=2)
abline(v=2.125,col="grey",lty=2)

set.immunity.dist.split()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
yy<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
points(-2.25,yy,col="black",pch=16,cex=1.5)

for (alpha in alphas)
{
  beta<-alpha  
  set.immunity.dist.beta(alpha,beta)
  initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
  get.startconds()
  yy<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
  points(log10(alpha),yy,col=plot.colors[which(alphas==alpha)],pch=16,cex=1.5)
}

set.immunity.dist.single()
initial.infection.structure<-function(x) {0.0*scaled.immunity.distribution(x)}
get.startconds()
yy<-var.func(scaled.immunity.distribution(dummy.immunity.categories),immunity.categories)
points(2.25,yy,col="darkgrey",pch=16,cex=1.5)


