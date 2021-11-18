######################################################
######   Codes du chapitre 8
######   Y. Aragon
######################################################

# Mise à jour effectuée le 17/08/2016


# packages utilisés dans le chapitre
require(caschrono)
require(fBasics)
require(urca)
require(polynom)
require(dse)
require(xtable)


# 8.1 Préparation des données
aa=read.table(file= system.file("/import/trafquoti.txt",package="caschrono"),
header=FALSE,quote="",sep="", colClasses=c("numeric","character"),
col.names =c("trafic","date"))
str(aa)
summary(aa)

date.1 = as.Date(aa$date)
date.1[1:10]
date.2=seq(from=as.Date("1993-01-01"),to=as.Date("2007-10-31"),by="day")
c(length(date.1),length(date.2))

# agrégation par an et par mois
an =  substr(aa$date,1,4)
mois =substr(aa$date,6,7)
trafan =  aggregate(aa$traf, list( An =  an), sum)
str(trafan)
trafan = ts(trafan$x/1000, start= 1993, frequency=1)
trafan.1 = window(trafan, end= 2006)

mois.an= as.numeric(paste(an,mois,sep=""))
trafmens=aggregate(aa$traf, list(Mois.An = mois.an), sum)
str(trafmens)
trafmensu=ts(trafmens$x/1000,start= c(1993,1),frequency= 12)



# 8.2 Exploration
# 8.2.1 Décomposition de la série en tendance, saisonnalité et erreur
dec.m = decompose(trafmensu)
op=par(oma=rep(0,4))
plot(dec.m)
abline(v=2001.75)
par(op)



# 8.2.2 Month plot
op=par(oma=rep(0,4))
monthplot(trafmensu, main="Month plot du trafic")
par(op)


window(trafmensu, start=c(2000,1), end=c(2002,12))


# 8.2.3 Lag plot
trafav=window(trafmensu, start=c(1996,1), end=c(2001,8))
traf.apr01=window(trafmensu, start=c(2001,10), end=c(2007,10))
time.apr01=time(traf.apr01)

op=par(oma=rep(0,4))
lag.plot(rev(trafav),set.lags=1:12,main="Avant le 9/11",
 asp=1,diag =TRUE, diag.col="red",type="p",do.lines=FALSE)
par(op)



op=par(oma=rep(0,4))
lag.plot(rev(traf.apr01), set.lags=1:12,main="Après le 9/11",
 asp=1,diag=TRUE,diag.col="red",type="p",do.lines=FALSE)
par(op)



# 8.3 Modélisation avant septembre 2001
aa=dagoTest(trafav)
aa@test$p.value[1]

temps = time(trafav)
(mod1 = lm(trafav~temps))
resid= ts(residuals(mod1),start=c(1996,1),frequency=12)

op=par(mfrow=c(3,1),oma=rep(0,4))
par(mar=c(0,4,1,0))
plot(trafav,xlab="", ylab="trafic",xaxt="n")
par(mar=c(4,4,0,0))
plot(resid, xlab="temps",ylab="résidu MCO")
abline(h=0)
par(mar=c(4,4,2,0))
acf(resid, lag.max=50, xlab="retard",main="")
par(op)


# 8.3.1 Modélisation manuelle
op=par(oma=rep(0,4))
xy.acfb(diff(trafav,12),40,numer=FALSE)
par(op)


require("forecast")
(mod2=Arima(trafav,seasonal=list(order=c(1,1,0),period=12),
  include.drift=TRUE))

op=par(oma=rep(0,4))
xy.acfb(mod2$residuals,numer=FALSE)
par(op)


(mod3=Arima(trafav,order=c(2,0,0),seasonal=list(order=c(1,1,0),period=12),
 method="ML",include.drift=TRUE))

op=par(oma=rep(0,4))
xy.acfb(residuals(mod3),numer=FALSE)
par(op)


ret= seq(6,30,6)
t(Box.test.2(residuals(mod3),ret,type="Ljung-Box",fitdf=4,decim=4))

mod4=Arima(trafav,order=c(2,0,0),
  seasonal=list(order= c(1,1,0),period=12),
  fixed=c(0,NA,NA,NA),method="ML",include.drift=TRUE)
summary(mod4)

ret6= seq(6,30,6)
t(Box.test.2(residuals(mod4),ret6,type="Ljung-Box",fitdf=3,decim=4))



# 8.3.2 Modélisation automatique
best.av = auto.arima(trafav, D = 1)
summary(best.av)

t(Box.test.2(residuals(best.av),ret6,type="Ljung-Box",decim=4,fitdf =4))
t_stat(best.av)

ec80 = mod4$sigma2^.5 * qnorm(0.90)
vajust = fitted(mod4)
matri=as.ts(cbind(trafav,vajust-ec80,vajust+ec80),start=c(1996,1),
   frequency=12)

op=par(oma=rep(0,4))
plot(matri, plot.type="single", lty=c(1,2,2), xlab="temps",
 ylab="trafic", main="", cex.main=0.8 )
legend( par("usr")[1], par("usr")[4], c("Valeur observée","Bande de prédiction"), lwd=1, lty=c(1,2) )
par(op)


indi = (trafav - (vajust-ec80)) > 0 & ( vajust+ec80 - trafav) > 0
prop = 100*sum(indi)/length(indi)

ur1 = ur.df(diff(trafav,12),lags=1,type="drift")
summary(ur1)

require("tseries")
kpss.test(diff(trafav,12))


# 8.4 Impact sur le volume de traffic
# 8.4.1 Prévision ponctuelle
prev2002 = forecast(mod4, h=16, level=80)

apredire=window(trafmensu,start=c(2001,9),end=c(2002,12))
matri2 = cbind(apredire, prev2002$lower, prev2002$upper,prev2002$mean)

plot(matri2,plot.type="single",lty=c(1,2,2,3),
    xlab="temps",xaxt="n",ylab="trafic")
axis(1,at=c(2001.7,2002,2002.9),
  labels=c("sept-2001","jan-2002","dec-2002"),
  cex.axis=.8)
legend("topleft",c("Valeur obs.","Bande de préd. à 80%","Préd. moyenne"),lwd=1, lty=c(1,2,3),cex=.8 )



pr2002 = prev2002$mean[5:16]

sigpred = (prev2002$upper - prev2002$mean)/qnorm(.9)
ordres.q = pnorm((apredire - prev2002$mean)/sigpred)
round(100*ordres.q,digits=2)


# 8.4.2 Simulation de trajectoires
poly.ar = polynomial(coef=c(1,0,-mod4$coef[2]))*polynomial(coef=c(1,rep(0,11),-mod4$coef[3]))

cte = predict(poly.ar,1)*12*mod4$coef[4]

fac.rsaiso = polynomial(coef = c(1,rep(0,11), -1))
coef.yg = fac.rsaiso*poly.ar

poly.ma=1

traf.ini= window(trafmensu, start=c(1999,7), end=c(2001,8))
y0= traf.ini; ly0 = length(y0)
nsim=10000 ; hori= 16
ysim= matrix(0, ncol=nsim, nrow=hori)
set.seed(347) # choix d'une graine 
bb=matrix(rnorm(nsim*hori,sd=mod4$sigma2^.5),nrow=hori,ncol=nsim)
AR =  array(as.vector(coef.yg),c(length(coef.yg),1,1))
Cte = array(cte,c(1,1,1))
entree = as.matrix(rep(1,hori))
BM = array(as.vector(poly.ma),c(length(poly.ma),1,1))
mod.sim = ARMA(A=AR, B= BM, C=Cte)

for (sim in 1:nsim)
{ysim[,sim] = simulate(mod.sim, y0=rev(y0), input=entree,
     sampleT=hori, noise=as.matrix(bb[,sim]))$output }

moy.mois = apply(ysim[5:16,],1,mean)
var.mois = apply(ysim[5:16,],1,var)
compar=round(rbind(moy.mois,prev2002$mean[5:16],
             var.mois^.5, sigpred[5:16]),digits=2)

aa= seq(as.Date("2000/1/1"), by="month", length.out=12)
id.mois= months(aa, abbreviate = TRUE)
rownames(compar) = c("moy. emp.","moy. théo.",
            "é.-types emp.", "é.-types théo.")
colnames(compar) = id.mois
paires = seq(2,12,by=2)

xtable(compar[,paires], caption="Moyennes et écarts-types des prévisions - simulations et calcul théorique.", digits=2,label="comparsimest")



# Distribution de la perte en 2002
traf.02= window(trafmensu, start=c(2002,1), end=c(2002,12))
perte.sim=ysim[5:16,]-matrix(traf.02,nrow=12,ncol=nsim)
perte.an = apply(perte.sim,2,sum)
perte.an.pct = 100*(perte.an / sum(traf.02))
perte= rbind(perte.an, perte.an.pct)
q.perte=t(apply(perte,1,quantile,probs=c(0,.10,.25,.5,.75,.90,1)))

xtable(q.perte, caption="Quantiles de la distribution des pertes en valeur (milliers de passagers) et en pourcentage pour 2002.", digits=2,label="quant.pert")


# 8.5 Etude après le 9/11 - lissage exponentiel
op=par(oma=rep(0,4))
plot(traf.apr01, xaxt="n",xlab="", ylab="trafic")
axis(1,at=time.apr01[c(4,16,28,40,52,64)],labels=c("2002","2003",
     "2004","2005","2006","2007"), cex.axis=.8)
par(op)


traf0405 = window(trafmensu,start=c(2004,1),end=c(2005,12))

traf06_1_6=window(trafmensu,start=c(2006,1),end=c(2006,6))

es.1 =  ets(traf0405,model="ANA")
summary(es.1)


# Simulation d'un modèle de lissage exponentiel
saiso=matrix(c(rep(0,11),1), nrow=1)
saiso=rbind(saiso, cbind(diag(x=1,ncol=11,nrow=11),
            matrix(rep(0,11),ncol=1)))
# blocs ligne de la matrice de transition
f1 = matrix(c(1,rep(0,12)), nrow=1)
f2 = cbind(matrix(rep(0,12),ncol=1),saiso)
fmat = rbind(f1,f2)
gmat=NULL
kmat = matrix(c(es.1$fit$par[1:2],rep(0,11)), ncol=1)
hmat = matrix(c(1,rep(0,11),1),nrow=1)
mod.es=SS(F=fmat,G=gmat,K=kmat,H=hmat,z0=es.1$states[nrow(es.1$states)-1,])

nsim=10000 ; resulsim = matrix(0,ncol= nsim, nrow=6)
set.seed(975)
for(i in 1:nsim)
{resulsim[,i]=simulate(mod.es,
  noise=list(w=matrix(rnorm(6,sd=es.1$sigma2^.5),ncol=1),
          w0=rnorm(1,sd=es.1$sigma2^.5)))$output}


q2006 = apply(resulsim,1, quantile, probs= c(0.05,.10,.90,.95))
mat.rep = t(rbind(q2006,as.numeric(traf06_1_6)))
op=par(oma=rep(0,4))
matplot(1:6,mat.rep,type="l",lty=c(1,2,2,3,3),ylab="trafic",
   xlab="temps",col="black",lwd=1.5)
title(main="Bandes de prédiction basées sur les quantiles",
  sub="Premier semestre 2006")
leg.txt=c("Série",expression(q[0.05]),expression(q[0.10]),
  expression(q[0.90]),expression(q[0.95]))
legend(1,540, leg.txt, lty=c(1,2,2,3,3))
par(op)


# 8.6 Estimation d'un SARIMA dans R - Vérification
xmat = as.matrix(1:length(trafav))
(mod4x=Arima(trafav,order=c(2,0,0),
  seasonal=list(order=c(1,1,0),period=12),
  fixed = c(0,NA,NA,NA), xreg= xmat))
