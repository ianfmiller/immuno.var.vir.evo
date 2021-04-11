rm(list=ls())

# parameter sets to explore
params.mat<-data.frame("p.immune"=c(1,1,1,1,1),"w"=c(1,0,0,0,0),"x"=c(0,1,0,0,.5),"y"=c(0,0,1,0,0),"z"=c(0,1,0,1,1))
alphas<-10^(seq(-2,2,.25))

for(i in 1:dim(params.mat)[1])
{
  p.immune<-params.mat[i,"p.immune"]
  w<-params.mat[i,"w"]
  x<-params.mat[i,"w"]
  y<-params.mat[i,"w"]
  z<-params.mat[i,"w"]
  
  source("~/immuno.var.vir.evo/code/folder.setup.R")
  setwd(new.dir)
  
  cmd<-paste0("jid1=$(sbatch ",paste0("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,".q;"),
              " jid2=$(sbatch --dependency=afterany:${jid1##* } ",paste0("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,".data.q"))
  system(cmd)
}