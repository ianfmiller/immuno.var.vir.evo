susceptible.eqns<-c()

for (i in 1:n.immunity.categories) {susceptible.eqns<-paste(susceptible.eqns,paste("startcondsdot[",i-1,"] = birthrate * N * (1 - N/K) * ((init",i-1,"+init",i+(n.immunity.categories-1),"+init",i+(2*n.immunity.categories-1),"+init",i+(3*n.immunity.categories-1),")/Ninit) - deathrate * startconds[",i-1,"] - startconds[",i-1,"] * (",paste("betafunc(cat[",i-1,"],cat[",seq(0,n.immunity.categories-1),"],v1)*startconds[",seq(n.immunity.categories,(2*n.immunity.categories-1)),"]",sep="",collapse=" + "),") - startconds[",i-1,"] * (",paste("betafunc(cat[",i-1,"],cat[",seq(0,n.immunity.categories-1),"],v2)*startconds[",seq((n.immunity.categories*2),(n.immunity.categories*3-1)),"]",sep="",collapse=" + "),");",sep=""),sep='\n')}

inf.eqns.res<-c()
for (i in (n.immunity.categories+1):(2*n.immunity.categories)) {inf.eqns.res<-paste(inf.eqns.res,paste("startcondsdot[",i-1,"] = startconds[",i-1,"] * -1*(deathrate + drdfunc(cat[",i-1-n.immunity.categories,"],v1) + gammafunc(cat[",i-1-n.immunity.categories,"],v1)) + startconds[",i-1-n.immunity.categories,"] * (",paste("betafunc(cat[",i-n.immunity.categories-1,"],cat[",seq(0,n.immunity.categories-1),"],v1)*startconds[",seq((n.immunity.categories),(n.immunity.categories*2-1)),"]",sep="",collapse=" + "),");",sep=""),sep='\n')}

inf.eqns.inv<-c()
for (i in (2*n.immunity.categories+1):(3*n.immunity.categories)) {inf.eqns.inv<-paste(inf.eqns.inv,paste("startcondsdot[",i-1,"] = startconds[",i-1,"] * -1*(deathrate + drdfunc(cat[",i-1-2*n.immunity.categories,"],v2) + gammafunc(cat[",i-1-2*n.immunity.categories,"],v2)) + startconds[",i-1-2*n.immunity.categories,"] * (",paste("betafunc(cat[",i-2*n.immunity.categories-1,"],cat[",seq(0,(n.immunity.categories-1)),"],v2)*startconds[",seq((n.immunity.categories*2),(n.immunity.categories*3-1)),"]",sep="",collapse=" + "),");",sep=""),sep='\n')}

rec.eqns<-c()
for(i in (3*n.immunity.categories+1):(4*n.immunity.categories)) {rec.eqns<-paste(rec.eqns,paste("startcondsdot[",i-1,"] = startconds[",i-2*n.immunity.categories-1,"] * gammafunc(cat[",i-1-3*n.immunity.categories,"],v1) + startconds[",i-n.immunity.categories-1,"] * gammafunc(cat[",i-1-3*n.immunity.categories,"],v2) - startconds[",i-1,"] * deathrate;",sep=""),sep='\n')}

Fmat.res.eqns<-c()
for (i in 1:n.immunity.categories) {
  for (j in 1:n.immunity.categories)
  {
    Fmat.res.eqns<-paste(Fmat.res.eqns,paste("out[",(i-1)*n.immunity.categories+j-1,"] = startconds[",i-1,"]*betafunc(cat[",i-1,"],cat[",j-1,"],v1);",sep=""),sep='\n')
  }
}

Vmat.res.eqns<-c()
for (i in 1:n.immunity.categories) {
  for (j in 1:n.immunity.categories)
  {
    Vmat.res.eqns<-paste(Vmat.res.eqns,ifelse(i==j,paste("out[",(i+n.immunity.categories-1)*n.immunity.categories+j-1,"] = deathrate + gammafunc(cat[",i-1,"],v1) + drdfunc(cat[",i-1,"],v1);",sep=""),paste("out[",(i+n.immunity.categories-1)*n.immunity.categories+j-1,"] = 0;",sep="")),sep='\n')
  }
}

Fmat.inv.eqns<-c()
for (i in 1:n.immunity.categories) {
  for (j in 1:n.immunity.categories)
  {
    Fmat.inv.eqns<-paste(Fmat.inv.eqns,paste("out[",(i+2*n.immunity.categories-1)*n.immunity.categories+j-1,"] = startconds[",i-1,"]*betafunc(cat[",i-1,"],cat[",j-1,"],v2);",sep=""),sep='\n')
  }
}

Vmat.inv.eqns<-c()
for (i in 1:n.immunity.categories) {
  for (j in 1:n.immunity.categories)
  {
    Vmat.inv.eqns<-paste(Vmat.inv.eqns,ifelse(i==j,paste("out[",(i+(3*n.immunity.categories-1))*n.immunity.categories+j-1,"] = deathrate + gammafunc(cat[",i-1,"],v2) + drdfunc(cat[",i-1,"],v2);",sep=""),paste("out[",(i+(3*n.immunity.categories-1))*n.immunity.categories+j-1,"] = 0;",sep="")),sep='\n')
  }
}

cat.string<-c(immunity.categories[1])
for (i in 2:n.immunity.categories) {cat.string<-paste(cat.string,immunity.categories[i],sep=",")}

Salls<-c("startconds[0]")
for(i in 2:n.immunity.categories) {Salls<-paste(Salls, paste("startconds[",i-1,"]",sep=""),sep=" + ")}

Ialls<-paste("startconds[",n.immunity.categories,"]",sep="")
for(i in (n.immunity.categories+1):(2*n.immunity.categories)) {Ialls<-paste(Ialls, paste("startconds[",i-1,"]",sep=""),sep=" + ")}

Ralls<-paste("startconds[",2*n.immunity.categories,"]",sep="")
for(i in (2*n.immunity.categories+2):(3*n.immunity.categories)) {Ralls<-paste(Ralls, paste("startconds[",i-1,"]",sep=""),sep=" + ")}

init.parms<-c()
for(i in 1:(n.immunity.categories*4))
{init.parms<-paste(init.parms,paste("# define init",i-1," parms[",13+i,"]",sep=""),sep="\n")}

Ninit<-c("init0")
for(i in 2:(4*n.immunity.categories)) {Ninit<-paste(Ninit,paste(" + ","init",i-1,sep=""))}

text<-{
  paste(
    "/* file sirmodess.c */
#include <R.h>
#include <math.h>
static double parms[",14+n.immunity.categories*4,"];
#define birthrate parms[0]
#define deathrate parms[1]
#define K parms[2]
#define v1 parms[3]
#define v2 parms[4]
#define b1 parms[5]
#define b2 parms[6]
#define d1 parms[7]
#define d2 parms[8]
#define c1 parms[9]
#define w parms[10]
#define x parms[11]
#define y parms[12]
#define z parms[13]

",init.parms,
    
    "
    
double cat[",n.immunity.categories,"] = {",cat.string,"};
    
/* wfunc */
double wfunc (double immunityclass) {
double result;
result = immunityclass*w;",
    "
return result; 
    }
    
    "
    ,
    "/* xfunc */
double xfunc (double immunityclass) {
double result;
result = immunityclass*x;",
    "
return result; 
    }
    
    "
    ,
    "/* yfunc */
double yfunc (double immunityclass) {
double result;
result = immunityclass*y;",
    "
return result; 
    }
    
    "
    ,
    "/* zfunc */
double zfunc (double immunityclass) {
double result;
result = immunityclass*z;",
    "
return result; 
}
    
"
,
"double gammafunc (double immunityclass, double virulence) {
double result;
result = d1*pow((1-xfunc(immunityclass))*virulence,d2);
return result; 
}

"
,
"double betafunc (double susimmunityclass, double infimmunityclass, double virulence) {
double result;
result = (1-wfunc(susimmunityclass))*(1-yfunc(infimmunityclass))*(b1*pow((1-xfunc(infimmunityclass))*virulence,b2));
return result; 
}

"
,
"double drdfunc (double immunityclass, double virulence) {
double result;
result = (1-zfunc(immunityclass))*virulence*c1; 
return result; 
}

"
,
"/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
  int N=",14+n.immunity.categories*4,";
  odeparms(&N, parms);
}

",
'/* ODEs */
void derivs (int *neq, double *t, double *startconds,double *startcondsdot, double *out, int *ip)
{
  if (ip[0] <1) error("nout should be at least 1");
  ',
  "
double Sall;
double Iall;
double Rall;
double N;
double Ninit;
  
Sall = ", paste(Salls,";"), 
  "
Iall = ",paste(Ialls,";"),  
  "
Rall = ",paste(Ralls,";"),
  "
N = Sall + Iall + Rall;",
  "
Ninit = ",paste(Ninit,";"),
  "
  
  /* susceptible classes */
  ",
  susceptible.eqns,
  "
  
  /* infected classes-resident */
  ",
  inf.eqns.res,
  "
  
  /* infected classes-invader */
  ",
  inf.eqns.inv,
  "
  
  /* recovered classes */
  ",
  rec.eqns,
  "
  
  /* output for Fmat.res */
  ",
  
  Fmat.res.eqns,
  "
  
  /* output for Vmat.res */
  ",
  Vmat.res.eqns,
  "
  
  /* output for Fmat.inv */
  ",
  
  Fmat.inv.eqns,
  "
  
  /* output for Vmat.inv */
  ",
  Vmat.inv.eqns,
  
  '
  
}

/* END file mymod.c */
'
,
sep="")  
}
write(text,file="sirmodessDD.c")
