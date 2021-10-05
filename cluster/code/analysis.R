### ESS analysis ###
## to iterate on cluster
library(deSolve)
library(parallel)
source("analysis.setup.R")

w<-0
x<-0
y<-0
z<-0
p.immune<-1
theta<-.01

### demographic parameters
birth.rate=.07
death.rate=.05

### trade-off parameters
K<-1/(-death.rate/birth.rate+1)
b1<<-3
b2<<-.333
d1<<-.1
d2<<-0
c1<<-1

virulence.steps<-c(seq(0,.1,.005),seq(.10005,.25,by=.0005),seq(.25,.5,.01))
set.nonlinear.tradeoff(b1,b2,d1,d2,c1)

source("writeSIResscDD.R")
system("R CMD SHLIB sirmodessDD.c")
dyn.load(paste("sirmodessDD", .Platform$dynlib.ext, sep = ""))

analysis.func<-function(i)
{
  v1<-i
  result<-data.frame("v1"=numeric(),"v2"=numeric(),"RE.res.val"=numeric(),"RE.inv.val"=numeric(),"xindex"=numeric(),"yindex"=numeric())
  #setup
  if (theta>0) {set.immunity.dist.beta(theta,theta)}
  if(theta==-1) {set.immunity.dist.split()}
  if(theta==-2) {set.immunity.dist.single()}
  
  #simulate epidemic & pull equilibirum parameters
  v2<-0
  initial.infection.structure<-function(x) {0.01*scaled.immunity.distribution(x)}
  get.startconds()
  parms <-c(birth.rate=birth.rate,death.rate=death.rate,K=K,v1=v1,v2=v2,b1=b1,b2=b2,d1=d1,d2=d2,c1=c1,w=w,x=x,y=y,z=z,startconds)
  out <- ode(startconds, seq(0,2000,.1), func = "derivs", parms = parms,dllname = "sirmodessDD",initfunc = "initmod", nout = 4*n.immunity.categories^2, outnames = paste("out",seq(0,4*n.immunity.categories^2-1,1),sep=""),verbose=F,method="lsoda")
  
  startconds.epi.equi<-as.numeric(out[dim(out)[1],2:(4*n.immunity.categories+1)])
  names(startconds.epi.equi)<-names(startconds)
  
  for (v2 in virulence.steps)
  {
    parms2 <-c(birth.rate=birth.rate,death.rate=death.rate,K=K,v1=v1,v2=v2,b1=b1,b2=b2,d1=d1,d2=d2,c1=c1,w=w,x=x,y=y,z=z,startconds) #use startdonds to set birth rates as this is carry over from first simulation. Use startconds.epi.equi for starting conditions.
    out <- ode(startconds.epi.equi, c(0,0), func = "derivs", parms = parms2,dllname = "sirmodessDD",initfunc = "initmod", nout = 4*n.immunity.categories^2, outnames = paste("out",seq(0,4*n.immunity.categories^2-1,1),sep=""),verbose=F,method="lsoda")
    get.matricies(output=out)
    RE.res<-getR0(Fmat.res,Vmat.res,output=F)
    RE.inv<-getR0(Fmat.inv,Vmat.inv,output=F)
    result<-rbind(result,cbind(v1,v2,RE.res,RE.inv))
  }
  result
}

raw.analysis.out<-mcmapply(analysis.func,virulence.steps,SIMPLIFY = F,mc.cores = 4)
analysis.out<-data.frame("RE.res.val"=numeric(),"RE.inv.val"=numeric(),"xindex"=numeric(),"yindex"=numeric())

for (i in 1:length(raw.analysis.out))
{
  analysis.out<-rbind(analysis.out,raw.analysis.out[[i]])
}

RE.inv.mat<-matrix(NA,length(virulence.steps),length(virulence.steps)) #rows=res strategy #columns=invader strategy
colnames(RE.inv.mat)<-virulence.steps
row.names(RE.inv.mat)<-virulence.steps

RE.res.mat<-matrix(NA,length(virulence.steps),length(virulence.steps)) #rows=res strategy #columns=invader strategy
colnames(RE.res.mat)<-virulence.steps
row.names(RE.res.mat)<-virulence.steps

for (i in 1:dim(analysis.out)[1])
{
  xindex<-which(virulence.steps==analysis.out[i,"v1"])
  yindex<-which(virulence.steps==analysis.out[i,"v2"])
  RE.res.mat[xindex,yindex]<-analysis.out[i,"RE.res"]
  RE.inv.mat[xindex,yindex]<-analysis.out[i,"RE.inv"]
}

out<-c(theta,ess.analysis(RE.inv.mat))
names(out)<-c("shape","PIP","ES vir","repeller?","repeller vir","hyper vir?","upper erad?","upper erad vir","lower erad?","lower erad vir","mid.erad.lower?","mid.erad.lower.vir","mid.erad.upper?","mid.erad.upper.vir")

# save file if working on cluster
if(dir.exists("~/immuno.var.vir.evo/output"))
{
  setwd("~/immuno.var.vir.evo/output")
  saveRDS(out,file=paste("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,"theta",theta,".RDS",sep=""))
}
