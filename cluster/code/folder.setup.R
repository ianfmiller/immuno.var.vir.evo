#setting up folders--to be sourced in SLURM

setwd("~/immuno.var.vir.evo")
new.dir<-paste("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,sep="")
if(!dir.exists(new.dir)) {dir.create(new.dir)}
setwd(new.dir)
new.dir<-getwd()
if (length(dir()) > 0) {for (name in dir()) {file.remove(name)}}

batch.lines<-readLines("~/immuno.var.vir.evo/code/batch.script.q")
batch.lines[6]<-paste("#SBATCH -J",paste('"p',p.immune,"W",w,"X",x,"Y",y,"Z",z,'"',sep=""))
batch.lines[14]<-paste("cd ~/immuno.var.vir.evo/",paste("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,sep=""),"/dir.$INDEX",sep="")
writeLines(batch.lines,paste("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,".q",sep=""))

data.batch.lines<-readLines("~/immuno.var.vir.evo/code/data.batch.script.q")
data.batch.lines[6]<-paste("#SBATCH -J",paste('"p',p.immune,"W",w,"X",x,"Y",y,"Z",z,'"',sep=""))
data.batch.lines[14]<-paste("cd ~/immuno.var.vir.evo/",paste("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,sep=""),"/dir.$INDEX",sep="")
writeLines(data.batch.lines,paste("p",p.immune,"W",w,"X",x,"Y",y,"Z",z,".data.q",sep=""))

data.compilation.lines<-readLines("~/immuno.var.vir.evo/code/data.compilation.R")
data.compilation.lines[grep("w<-",data.compilation.lines)]<-paste("w<-",w,sep="")
data.compilation.lines[grep("x<-",data.compilation.lines)]<-paste("x<-",x,sep="")
data.compilation.lines[grep("y<-",data.compilation.lines)]<-paste("y<-",y,sep="")
data.compilation.lines[grep("z<-",data.compilation.lines)]<-paste("z<-",z,sep="")
data.compilation.lines[grep("p.immune<-",data.compilation.lines)]<-paste("p.immune<-",p.immune,sep="")
data.compilation.lines[grep("alpha<-",data.compilation.lines)]<-paste("alpha<-",alphas[i],sep="")
writeLines(data.compilation.lines,"data.compilation.R")

for (i in 1:length(alphas))
{
  setwd(new.dir)
  dir.create(paste("dir.",i,sep=""))
  
  setwd(paste(new.dir,"/dir.",i,sep=""))
  
  analysis.lines<-readLines("~/immuno.var.vir.evo/code/analysis.R")
  
  analysis.lines[grep("w<-",analysis.lines)]<-paste("w<-",w,sep="")
  analysis.lines[grep("x<-",analysis.lines)]<-paste("x<-",x,sep="")
  analysis.lines[grep("y<-",analysis.lines)]<-paste("y<-",y,sep="")
  analysis.lines[grep("z<-",analysis.lines)]<-paste("z<-",z,sep="")
  analysis.lines[grep("p.immune<-",analysis.lines)]<-paste("p.immune<-",p.immune,sep="")
  analysis.lines[grep("alpha<-",analysis.lines)]<-paste("alpha<-",alphas[i],sep="")
  writeLines(analysis.lines,"analysis.R")
  
  file.copy("~/immuno.var.vir.evo/code/analysis.setup.R",paste(new.dir,"/dir.",i,sep=""))
  file.copy("~/immuno.var.vir.evo/code/writeSIResscDD.R",paste(new.dir,"/dir.",i,sep=""))
}
