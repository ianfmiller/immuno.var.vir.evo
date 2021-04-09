## data compilation
w<-0
x<-0.5
y<-0
z<-0
p.immune<-1
alpha<--1

out.data<-c()
setwd("~/immuno.var.vir.evo/output")

for (i in 1:19)
{
  data<-readRDS(paste("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,"alpha",alpha,".RDS",sep=""))
  out.data<-rbind(out.data,data)
}

write.csv(out.data,file=paste("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,".csv",sep=""))