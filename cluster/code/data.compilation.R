## data compilation
w<-0
x<-0.5
y<-0
z<-0
p.vacc<-1
alpha<--1

out.data<-c()
      
for (i in 1:19)
{
  setwd(paste("~/hetero/",paste("p",p.vacc,"W",w,"X",x,"Y",y,"Z",z,sep=""),"/dir.",i,sep=""))
  data<-readRDS(dir()[grep(".Rds",dir())])
  out.data<-rbind(out.data,data)
}

setwd(paste("~/hetero/",paste("p",p.vacc,"W",w,"X",x,"Y",y,"Z",z,sep=""),sep=""))
write.csv(out.data,file=paste("p",p.vacc,"W",w,"X",x,"Y",y,"Z",z,".csv",sep=""))