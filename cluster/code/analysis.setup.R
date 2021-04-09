# setup and functions for ESS analysis ###

## setup

### set immunity category midpoints for calculating frequencies
dummy.immunity.categories<-seq(.005,.995,length.out = 21)

### set values of immunity
immunity.categories<-seq(0,1,length.out = 21)

### set misc params
n.immunity.categories<-length(immunity.categories) #### number of immunity categories
alphas<-10^(seq(-2,2,.25)) #### distribution parameters


## functions

### get relative frequency of immune categories within total pop
scaled.immunity.distribution<-function(x) {immunity.distribution(x)/sum(mapply(immunity.distribution,dummy.immunity.categories))}

### set initial infection structure, assume 1% infected
initial.infection.structure<-function(x) {0.01*scaled.immunity.distribution(x)}

### set initial recovered structure, assume no recovered individuals
initial.recovered.structure<-function(x) {0*scaled.immunity.distribution(x)}

### get starting conditions for evo analyses
get.startconds<-function()
{
  startconds<-rep(NA,4*n.immunity.categories)
  
  for(n in 1:n.immunity.categories) {names(startconds)[n]<-paste("S",immunity.categories[n],sep="")}
  for(n in 1:n.immunity.categories) {names(startconds)[n.immunity.categories+n]<-paste("I",immunity.categories[n],sep="")}
  for(n in 1:n.immunity.categories) {names(startconds)[2*n.immunity.categories+n]<-paste("Ii",immunity.categories[n],sep="")}
  for(n in 1:n.immunity.categories) {names(startconds)[3*n.immunity.categories+n]<-paste("R",immunity.categories[n],sep="")}
  
  for (i in 1:n.immunity.categories) #Si +Ii
  {
    startconds[i]<-scaled.immunity.distribution(dummy.immunity.categories[i])-initial.infection.structure(dummy.immunity.categories[i])-initial.recovered.structure(dummy.immunity.categories[i])
    startconds[n.immunity.categories+i]<-initial.infection.structure(dummy.immunity.categories[i])
    startconds[2*n.immunity.categories+i]<-0
    startconds[3*n.immunity.categories+i]<-initial.recovered.structure(dummy.immunity.categories[i])
  }
  
  startconds<<-startconds
  
}

### get starting conditions for epi analyses/calculating IFRs
get.startconds.epi.effects<-function()
{
  get.startconds()
  startconds<<-startconds[-c((2*n.immunity.categories+1):(3*n.immunity.categories))]
  startconds<<-c(startconds,rep(0,times=2*n.immunity.categories))
  names(startconds)[(3*n.immunity.categories+1):(5*n.immunity.categories)]<<-c(paste("Rtot",immunity.categories,sep=""),paste("Dtot",immunity.categories,sep=""))
}

### plot starting conditions
plot.startconds<-function(ymax=1,main=paste("rate=",rate,sep=" "),border=NULL)
{
  bin.size<-immunity.categories[2]-immunity.categories[1]
  plot(0,0,col="white",xlim=c(-.05,1.05),ylim=c(0,ymax),axes = T,xlab="immunity",ylab="frequency",main=main)
  for(i in 1:n.immunity.categories)
  {
    rect(immunity.categories[i]-bin.size*.5,0,immunity.categories[i]+bin.size*.5,startconds[i],col = "black",lty=0,border=border)
    rect(immunity.categories[i]-bin.size*.5,startconds[i],immunity.categories[i]+bin.size*.5,startconds[i]+startconds[n.immunity.categories+i],col = "green",lty=0,border=border)
    rect(immunity.categories[i]-bin.size*.5,startconds[i]+startconds[n.immunity.categories+i],immunity.categories[i]+bin.size*.5,startconds[i]+startconds[n.immunity.categories+i]+startconds[2*n.immunity.categories+i],col = "blue",lty=0,border=border)
  }
}

### define effects of immunity

set.immunity.effects<-function(w,x,y,z)
{
  w<<-function(sus.immunity.class) {w}
  x<<-function(sus.immunity.class) {x}
  y<<-function(inf.immunity.class) {y}
  z<<-function(sus.immunity.class) {z}
}

### set trade-off shape

set.nonlinear.tradeoff<-function(b1,b2,d1,d2,c1) #creates nonlinear trade-off between transmission and virulence
{  
  beta.func<<-function(virulence,sus.immunity.class,inf.immunity.class) {(1-w(sus.immunity.class))*(1-y(inf.immunity.class))*(b1*((1-x(sus.immunity.class))*virulence)^b2)}
  gamma.func<<-function(virulence,inf.immunity.class) {d1*((1-x(inf.immunity.class))*virulence)^d2}
  death.rate.disease.func<<-function(virulence,immunity.class) {(1-x(immunity.class))*(1-z(immunity.class))*virulence*c1}
  b1<<-b1
  b2<<-b2
  d1<<-d1
  d2<<-d2
  c1<<-c1
}

### plot trade-off

plot.tradeoffs<-function() #plots association between viruence, disease associated mortality, recovery, and transmission
{
  par(mfrow=c(2,2))
  plot(virulence.steps,mapply(beta.func,virulence=virulence.steps,sus.immunity.class=0,inf.immunity.class=0),type="l",xlab="virulence",ylab="beta")
  plot(virulence.steps,mapply(death.rate.disease.func,virulence=virulence.steps,immunity.class=0),type="l",xlab="virulence",ylab="death.rate.disease")
  plot(virulence.steps,mapply(gamma.func,virulence=virulence.steps,inf.immunity.class=0),type="l",xlab="virulence",ylab="gamma")
  plot(mapply(beta.func,virulence=virulence.steps,sus.immunity.class=0,inf.immunity.class=0),mapply(death.rate.disease.func,virulence=virulence.steps,immunity.class=0),type="l",xlab="beta",ylab="death.rate.disease")      
}

### calculates F and V matricies (used in next-gen R0 calculation) from output of SIR model
get.matricies<-function(output=out)
{
  Fmat.res<<-matrix(output[1,(4*n.immunity.categories+2):(4*n.immunity.categories+1+n.immunity.categories^2)],n.immunity.categories,n.immunity.categories,byrow = T)
  Vmat.res<<-matrix(output[1,(4*n.immunity.categories+2+n.immunity.categories^2):(4*n.immunity.categories+1+2*n.immunity.categories^2)],n.immunity.categories,n.immunity.categories,byrow=T)
  Fmat.inv<<-matrix(output[1,(4*n.immunity.categories+2+2*n.immunity.categories^2):(4*n.immunity.categories+1+3*n.immunity.categories^2)],n.immunity.categories,n.immunity.categories,byrow=T)
  Vmat.inv<<-matrix(output[1,(4*n.immunity.categories+2+3*n.immunity.categories^2):(4*n.immunity.categories+1+4*n.immunity.categories^2)],n.immunity.categories,n.immunity.categories,byrow=T)
}

### calculates R0 from F and V matricies
getR0<-function(Fmat,Vmat,output=T)
{values<-eigen(Fmat %*% solve(Vmat))$values
R0<<-0
for(i in 1:length(values))
{
  if(Im(values)[i]==0) {R0<<-max(R0,Re(values)[i])}
}
if(output==T) {print(R0)}
return(R0)
}

### gets disease associated mortality rates for each immunity category, as well as average across all categories
#### additional values returned: CFRvacc=CFR for individuals with any immunity; CFRnonvacc=CFR for individuals with no immunity

getCFRs<-function()
{
  CFRs<<-c()
  
  for (k in 1:n.immunity.categories) 
  {
    new.cfr<-death.rate.disease.func(v1,immunity.categories[k])/(death.rate.disease.func(v1,immunity.categories[k])+death.rate+gamma.func(v1,immunity.categories[k]))
    CFRs<-c(CFRs,new.cfr)
  }
  
  CFRall<-sum(CFRs*startconds[1:n.immunity.categories])
  
  if (p.immune<1)
  {CFRvacc<-(sum(CFRs[2:n.immunity.categories]*startconds[2:n.immunity.categories])+ifelse(startconds[1]>0,CFRs[1]*startconds[1]*(1-(1-p.immune)/startconds[1]),0))/p.immune}
  else {CFRvacc<-CFRall}
  
  CFRnonvacc<-(CFRs[1]*startconds[1]*((1-p.immune)/startconds[1]))/(1-p.immune) #equal to CFRs[1]
  
  CFRs<<-c(CFRs,CFRall,CFRvacc,CFRnonvacc)
  names(CFRs)[]<<-c(immunity.categories,"all","vacc","nonvacc")
}

### gets prevalence
getPrevs<-function()
{
  total.prevs<<-matrix(NA,dim(out)[1],2)
  total.prevs[,1]<<-times
  for (t in times) {
    total.prevs[which(total.prevs[,1]==t),2]<<-sum(as.numeric(out[which(times==t),(n.immunity.categories+2):(2*n.immunity.categories+1)]))
  }
  max.prev<<-max(total.prevs[,2])
  final.prev<<-total.prevs[dim(total.prevs)[1],2]
}

### Scans vector from start to finish to find 1st local maximum
find.peak<-function(x)
{
  peak.index<-NA
  for (z in 4:(length(x)-3)) 
  {
    if (x[z] >= x[z-1] && x[z] >= x[z-2] && x[z] > x[z-3] && x[z] >= x[z+1] && x[z] >= x[z+2] && x[z] > x[z+3]) {peak.index<-z;break}
  }
  return(peak.index)
}

### Used to find virulence strategy with fitness equal to that of local optimum 
find.break.even.point<-function(x,peak.index)
{
  break.even.point<-NA
  for (z in (peak.index+1):(length(x)))
  {
    if (x[z] >= x[peak.index]) {break.even.point<-z; break}
  }
  return(break.even.point)
}

### plots output of SIR model
plot.simulation<-function(ylim=c(0,2))
{  
  #par(mfrow=c(1,3))
  plot(times,rowSums(out[,2:(n.immunity.categories+1)]),type="l",col="black",xlab="time",ylab="freq",ylim=ylim)
  points(times,rowSums(out[,(n.immunity.categories+2):(2*n.immunity.categories+1)]),type="l",col="green")
  if (sum(out[1,(2*n.immunity.categories+2):(3*n.immunity.categories+1)]) > 0) {points(times,rowSums(out[,(2*n.immunity.categories+2):(3*n.immunity.categories+1)]),type="l",col="darkgreen")}
  points(times,rowSums(out[,(3*n.immunity.categories+2):(4*n.immunity.categories+1)]),type="l",col="red")
  
  #plot(rowSums(out[,2:(n.immunity.categories+1)]),rowSums(out[,(n.immunity.categories+2):(2*n.immunity.categories+1)]),xlab="S",ylab="I",type="l")
}

### introduces invader strain into population at prevalence p. Individuals infected with invader strain come from both infected and susceptible categories. 
introduce.invader.strain<-function(p)
{
  startconds.inv<-rep(NA,times=n.immunity.categories)
  names(startconds.inv)<-paste("i",names(startconds)[(n.immunity.categories+1):(2*n.immunity.categories)],sep="")
  startconds.new<-startconds
  for(l in 1:length(startconds.inv))
  {
    startconds.inv[l]<-p*startconds[l] #makes x% of sus individuals infected with invader strain
    startconds.new[l]<-(1-p)*startconds.new[l]
  }
  startconds.new[(2*n.immunity.categories+1):(3*n.immunity.categories)]<-startconds.inv
  startconds<<-startconds.new
}

#### sets immunity to be beta distributed
#### when p.immune parameter is> 0, a fraction 1-p.immune of population has 0 immunity and a fraction p.immune has immunity defined by beta function
set.immunity.dist.beta<-function(alpha,beta,p.immune=1) {immunity.distribution<<-function(x)
{
  weights<-c()
  for (i in dummy.immunity.categories) {weights<-c(weights,dbeta(i,alpha,beta))}
  weight.adj<-sum(weights)
  
  fun<-function(xx) {if (xx==min(dummy.immunity.categories)) {(1-p.immune)+p.immune*dbeta(xx,alpha,beta)/weight.adj} else {p.immune*dbeta(xx,alpha,beta)/weight.adj}}
  mapply(fun,x)
}}

### sets immunity to be binary
#### when p.immune parameter is> 0, a fraction 1-p.immune of population has 0 immunity and a fraction p.immune has immunity defined by binary funciton

set.immunity.dist.split<-function(p.immune=1) {immunity.distribution<<-function(x)
{
  
  fun<-function(xx){
    dens<-0
    if (xx==min(dummy.immunity.categories)) {dens<-.5*p.immune+(1-p.immune)}
    if (xx==max(dummy.immunity.categories)) {dens<-.5*p.immune}
    dens}
  mapply(fun,x)
  
}}

### sets immunity to have no variation
#### when p.immune parameter is > 0, a fraction 1-p.immune of population has 0 immunity and a fraction p.immune has immunity

set.immunity.dist.single<-function(p.immune=1) {immunity.distribution<<-function(x)
{
  
  fun<-function(xx){
    dens<-0
    if (xx==min(dummy.immunity.categories)) {dens<-1-p.immune}
    if (xx==median(dummy.immunity.categories)) {dens<-p.immune}
    dens}
  mapply(fun,x)
  
}}

### modified version of function for plotting immunity distributions
plot.startconds.2<-function(ymax=1,col="black",main)
{
  bin.size<-immunity.categories[2]-immunity.categories[1]
  plot(0,0,col="white",xlim=c(-.05,1.05),ylim=c(0,ymax),axes = F,xlab="",ylab="",bty="n")
  mtext("immunity",side=1,cex=1,line=2.5)
  mtext("frequency",side=2,cex=1,line=2.5)
  axis(1,cex.axis=1.1)
  axis(2,cex.axis=1.1)
  mtext(main,side = 3,cex=1.1,line = .25)
  for(i in 1:n.immunity.categories)
  {
    rect(immunity.categories[i]-bin.size*.5,0,immunity.categories[i]+bin.size*.5,startconds[i],col = col,lty=0,border = NA)
    rect(immunity.categories[i]-bin.size*.5,startconds[i],immunity.categories[i]+bin.size*.5,startconds[i]+startconds[n.immunity.categories+i],col = "green",lty=0,border = NA)
    rect(immunity.categories[i]-bin.size*.5,startconds[i]+startconds[n.immunity.categories+i],immunity.categories[i]+bin.size*.5,startconds[i]+startconds[n.immunity.categories+i]+startconds[2*n.immunity.categories+i],col = "blue",lty=0,border = NA)
  }
}

### analysis of ESS
ess.analysis<-function(mat)
{
  output<-rep(NA,length.out=3)
  mod.mat<-matrix(NA,dim(mat)[1],dim(mat)[2])
  colnames(mod.mat)<-colnames(mat)
  rownames(mod.mat)<-rownames(mat)
  for(j in 1:dim(mod.mat)[1]) {
    for(k in 1:dim(mod.mat)[1]) {
      mod.mat[j,k] <- 1*(mat[j,k]>1)
    }}
  
  row.max<-find.peak(rowSums(mod.mat))
  row.min<-find.peak(-1*rowSums(mod.mat))
  
  col.max<-find.peak(colSums(mod.mat))
  col.min<-find.peak(-1*colSums(mod.mat))
  
  if (is.numeric(col.min)&&is.numeric(row.max)) {if(col.min-row.max<=1) {output[1]<-"ESS"; output[2]<-virulence.steps[col.min]}}
  
  
  if (isFALSE(any(diff(colSums(mod.mat))<0)) && any(mod.mat>0)) 
  {
    output[4]<-"selection for hypervirulence"
  }
  
  return(output)
}

var.func<-function(pdf,values)
{
  out<-0
  for (k in 1:length(pdf))
  {
    out<-out+pdf[k]*(values[k]-mean(values))^2
  }
  out
}



