#NEEDS TO BE REDONE
#setting up folders--to be sourced in SLURM
p.vacc<-.5
w<-0
x<-1
y<-0
z<-1

alphas<-10^(seq(-2,2,.25))

setwd("/Users/ifmiller/Documents/Princeton/immunological heterogeneity modeling/heterogeneity/cluster/")
new.dir<-paste("p",p.vacc,"W",w,"X",x,"Y",y,"Z",z,sep="")
if(!dir.exists(new.dir)) {dir.create(new.dir)}
setwd(new.dir)
new.dir<-getwd()
if (length(dir()) > 0) {for (name in dir()) {file.remove(name)}}

setwd("/Users/ifmiller/Documents/Princeton/immunological heterogeneity modeling/heterogeneity/cluster")
batch.lines<-readLines("batch.script.q")
batch.lines[6]<-paste("#SBATCH -J",paste('"p',p.vacc,"W",w,"X",x,"Y",y,"Z",z,'"',sep=""))
batch.lines[14]<-paste("cd ~/hetero/",paste("p",p.vacc,"W",w,"X",x,"Y",y,"Z",z,sep=""),"/dir.$INDEX",sep="")
setwd(new.dir)
writeLines(batch.lines,paste("p",p.vacc,"W",w,"X",x,"Y",y,"Z",z,".q",sep=""))

setwd("/Users/ifmiller/Documents/Princeton/immunological heterogeneity modeling/heterogeneity/cluster")
data.batch.lines<-readLines("data.batch.script.q")
data.batch.lines[6]<-paste("#SBATCH -J",paste('"p',p.vacc,"W",w,"X",x,"Y",y,"Z",z,'data"',sep=""))
data.batch.lines[14]<-paste("cd ~/hetero/",paste("p",p.vacc,"W",w,"X",x,"Y",y,"Z",z,sep=""),sep="")
setwd(new.dir)
writeLines(data.batch.lines,paste("p",p.vacc,"W",w,"X",x,"Y",y,"Z",z,".data.q",sep=""))

data.comp.lines<-readLines("/Users/ifmiller/Documents/Princeton/immunological heterogeneity modeling/heterogeneity/cluster/data.compilation.R")

data.comp.lines[2]<-paste("w<-",w,sep="")
data.comp.lines[3]<-paste("x<-",x,sep="")
data.comp.lines[4]<-paste("y<-",y,sep="")
data.comp.lines[5]<-paste("z<-",z,sep="")
data.comp.lines[6]<-paste("p.vacc<-",p.vacc,sep="")
data.comp.lines[7]<-paste("alpha<-",alphas[i],sep="")

writeLines(data.comp.lines,"data.compilation.R")

alphas<-c(-1,10^seq(-2,2,.25),-2)

for (i in 1:length(alphas))
{
  setwd(new.dir)
  dir.create(paste("dir.",i,sep=""))
  
  setwd(paste(new.dir,"/dir.",i,sep=""))
  
  analysis.lines<-readLines("/Users/ifmiller/Documents/Princeton/immunological heterogeneity modeling/heterogeneity/analysis.R")
  
  analysis.lines[6]<-paste("w<-",w,sep="")
  analysis.lines[7]<-paste("x<-",x,sep="")
  analysis.lines[8]<-paste("y<-",y,sep="")
  analysis.lines[9]<-paste("z<-",z,sep="")
  analysis.lines[10]<-paste("p.vacc<-",p.vacc,sep="")
  analysis.lines[11]<-paste("alpha<-",alphas[i],sep="")
  writeLines(analysis.lines,"analysis.R")
  
  file.copy("/Users/ifmiller/Documents/Princeton/immunological heterogeneity modeling/heterogeneity/analysis.functions.R",paste(new.dir,"/dir.",i,sep=""))
  file.copy("/Users/ifmiller/Documents/Princeton/immunological heterogeneity modeling/heterogeneity/writeSIResscDD.R",paste(new.dir,"/dir.",i,sep=""))
}
