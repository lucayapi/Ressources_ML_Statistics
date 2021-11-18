#######################################################################
########   Code Chapitre 1 - Séries temporelles avec R
########   Y. Aragon

# Mise à jour effectuée le 17/08/2016

# packages utilisés dans le chapitre
require(caschrono)
# require(fSeries)  # ce package a été supprimé sur le CRAN,
# mais cela ne doit pas affecter les codes ci-dessous (ed 1)

# 1.1 Exemples de séries temporelles
# Exemple 1.1 (Populations)
data(popfr)
op=par(mfrow=c(2,1),mar=c(2.5,2.5,2,2),mgp=c(1.5,.5,0),
  oma=c(0,0,0,0),cex.main=.8,cex.lab=.7,cex.axis=.7)
plot.ts(popfr,xlab="année", ylab="population",
   main ="Population française, 1846-1951")
plot.ts(uspop,xlab="année", ylab="population",
  main ="Population des Etats-Unis, 1790-1970")
par(op)

# Exemple 1.2 (Morts par accident)
data(m30)
op = par(mfrow=c(2, 1),mar=c(2,3,2,2),oma=c(1,1,0,0),        
 mgp=c(2,.4,0),cex.main=.75,cex.axis=.8,cex=.8,cex.lab=.8)
plot.ts(m30, xlab="année", ylab="nombre de morts",las=1, 
 main ="Accidents de la route en France, mortalité mensuelle")
polygon(c(1995,2005,2005,1995),c(340,340,910,910),lty=2)
debut = c(1995,1); fin = c(2004,12) # zoom
plot.ts(window(m30,start=debut,end=fin),
  xlab="année",ylab="nombre de morts", las=1,
  main ="Années 1995 à 2004")
par(op)



# Exemple 1.3 (Champagne)
aa=scan(system.file("/import/champagne_2001.txt",package="caschrono"))
champa = aa/1000; # en milliers de bouteilles
ytr=cbind(champa,log(champa))                  
colnames(ytr)=c("champagne","log(champagne)") 
ytr.ts = ts(ytr, start=c(2001,1), frequency=12)
op = par(oma=rep(0,4),mgp=c(2.5,.7,0),mar=c(4.5,3,2,2),
 cex.main=.7,cex.axis=.8,cex=.9,cex.lab=.9)
plot.ts(ytr.ts,xlab="temps",main="",oma.multi=c(4.5,3,.2, 0),
  mar.multi=c(0, 4, 0, .5),las=0)  
par(op)


# Exemple 1.5 (Lac Huron)
temps = time(LakeHuron)
reglin = lm(LakeHuron~temps)
resi.mco = residuals(reglin)
ychap = fitted(reglin)
v12=c(var(LakeHuron),var(resi.mco))

op=par(mar=c(4,3,2,2),mgp=c(2,.4,0),oma=c(0,0,0,0))  
plot.ts(LakeHuron, las=1, ylab="niveau", xlab="année")
abline(coef= coef(reglin))
s = c(12,48,59,80); 
segments(temps[s],ychap[s],temps[s],LakeHuron[s],lty=1,lwd=2)
y.gra=1.0009*par()$usr[3]
text(temps[1],y.gra,labels="obs.",cex=.8)
text(temps[s],rep(y.gra,length(s)),labels=s,cex=.8)
par(op)


op=par(mfrow=c(2,1),mar=c(3.5,3.5,2,1),mgp=c(2,.5,0),
 oma=c(0,0,0,0),cex.main=.8,cex.axis=.8,cex=.8,cex.lab=.8)
plot(as.vector(temps),resi.mco,type="l",xlab="année",ylab="résidu")
abline(h=0)
zero = rep(0, length(s))
segments(temps[s], zero, temps[s], resi.mco[s], lty=3,lwd=1 )
y.gra=0.9*par()$usr[3]
text(temps[1],y.gra,labels="obs.",cex=.8)
text(temps[s],rep(y.gra,length(s)),labels=s,cex=.8)
n =length(resi.mco)
plot(resi.mco[-n],resi.mco[-1],xlab="résidu en t-1",asp=1,
     ylab="résidu en t")
par(op)


require("forecast")
mod.lac=Arima(LakeHuron,order=c(1,0,0),xreg=temps,method="ML")
ychap = fitted(mod.lac)
resi.inno = residuals(mod.lac)
v23= c(var(resi.mco),var(resi.inno))

ychap = fitted(mod.lac)
op= par(mfrow=c(2,1),mar=c(3.5,3.5,2,1),mgp=c(2,.5,0),
 oma=c(0,0,0,0),cex.main=.8,cex.axis=.8,cex=.8,cex.lab=.8)
plot.ts(cbind(LakeHuron,ychap),lty=1:2,type="l",las=1,
 plot.type="single",ylab="niveau",xlab="année")
leg.txt = c("Série observée","Prédiction (MCG + erreur AR(1))")
legend(x=1876, y= 577, leg.txt, lty=1:2 , cex=.8)
n=length(resi.inno)
plot(resi.inno[-n],resi.inno[-1],xlab="résidu en t-1",asp=1,
     ylab="résidu en t")
title("Lag plot des résidus ARMAX")
par(op)



# 1.2 Graphiques pour les séries temporelles
# 1.2.1 Chronogramme
# 1.2.2 Lag plot
set.seed(2761)
innov1 = rnorm(290,sd=4.18)
y = arima.sim(list(order = c(12,0,1), ma=-.7, ar=c(rep(0,11),.9)),
  innov =innov1, n.start =50, n = 240) +50
y.ts = ts(y,frequency=12, start=c(1920,1))
ytr=cbind(y.ts,nottem)
colnames(ytr)=c("série simulée","température")

op=par(oma=rep(0,4))
plot.ts(ytr,main="",xlab="temps",cex=1.2)
par(op)

lag.plot(rev(nottem),12,layout=c(4,3),do.lines=FALSE,diag.col="red",
   main = "Températures mensuelles moyennes à Nottigham,
   "1920-1939", col.main = "blue")


lag.plot(rev(y),12,layout=c(4,3),ask=NA,
 do.lines=FALSE, diag.col="red",
 main=expression(paste(SARMA(0,1)(1,0)[12], " simulé")))


# 1.2.3 Month plot
op= par(mfrow=c(2,1),mar=c(3.5,3.5,2,1),mgp=c(2,.5,0), oma=c(0,0,0,0))
monthplot(nottem,  ylab="Température",main="",
 cex.main=1)
par(mar=c(3,3,1,1))
monthplot(y.ts,xlab="saison",ylab="y",main="",cex.main=1)
par(op)

# 1.3 Tendance, saisonnalité, résidus

# 1.4 Etapes et objectifs de l'analyse d'une série temporelle

# 1.5 Opérateur retard

