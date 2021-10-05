rm(list=ls())

# parameter sets to explore
params.mat<-data.frame("p.immune"=c(1),"w"=c(1),"x"=c(1),"y"=c(1),"z"=c(1))
thetas<-c(-2,-1,10^(seq(-2,2,.25)))

for(k in 1:dim(params.mat)[1])
{
  p.immune<-params.mat[k,"p.immune"]
  w<-params.mat[k,"w"]
  x<-params.mat[k,"x"]
  y<-params.mat[k,"y"]
  z<-params.mat[k,"z"]
  
  source("~/immuno.var.vir.evo/code/folder.setup.R")
  setwd(new.dir)
  
  cmd<-paste0("jid1=$(sbatch ",paste0("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,".q);"),
              " sbatch --dependency=afterany:${jid1##* } ",paste0("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,".data.q"))
  system(cmd)
}