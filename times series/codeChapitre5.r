################################################
#####  Codes du chapitre 5
##### Yves Aragon
#################################################

# Mise à jour effectuée le 17/08/2016

# packages utilisés dans le chapitre
require(caschrono)
require(urca)


# Introduction
# 5.1 Séries intégrées - Modèles ARIMA et SARIMA 
# 5.2 Construction d'un modèle SARIMA
data(nottem)

set.seed(2761)
innov1 = rnorm(290,sd=4.18)
y = arima.sim(list(order = c(12,0,1), ma=-.7, ar=c(rep(0,11),.9)),
  innov =innov1, n.start =50, n = 240) +50
y.ts = ts(y,frequency=12, start=c(1920,1))
ytr=cbind(y.ts,nottem)
colnames(ytr)=c("série simulée","température")

op=par(oma=rep(0,4))
plot2acf(nottem,y, main=c("ACF nottem","ACF SAR"))
par(op)


# 5.3 Non stationnarité stochastique ou déterministe
# 5.3.1 Test de non stationnarité : introduction et pratique
# Exemple 5.1
data(Raotbl3)
attach(Raotbl3,warn.conflicts=FALSE)

op=par(mar=c(3,3,0,1))
plot(lc, type = "l", xaxt = "n", xlab="temps", cex=.8)
axis(1, at = seq(from = 2, to = 98, by=12),
     labels = as.character(seq(from = 1967, to = 1991, by=3)))
par(op)


lc.df1 = ur.df(y=lc, lags=3, type="trend")
summary(lc.df1)

# Exemple 5.2
data(UKpppuip)
i1=ts(UKpppuip$i1,start=c(1971,1),frequency=4)
op=par(oma=rep(0,4),mar=c(4,4,0,1))
attach(UKpppuip,warn.conflicts=FALSE)
plot(i1,type="l",xlab="temps",ylab="taux")
par(op)

i1.df0 = ur.df(y=i1, lags=6, type="drift")
summary(i1.df0)

i1.df1 = ur.df(y=i1, lags=1, type="drift")
summary(i1.df1)


# 5.3.2 Test de stationnarité à une tendance déterministe près
# Exemple 5.3
set.seed(231)
x = rnorm(1000)
y = cumsum(x)
xy = ts(cbind(x,y))
colnames(xy)=c("BB","M. aléat.")

plot.ts(xy, xlab="temps",main="",
  oma.multi=c(0,0,.2,0),mar.multi=c(0,4,0,.5),cex=.8)

summary(ur.kpss(x, type = "mu"))

summary( ur.kpss(x, type = "tau"))

summary(ur.kpss(y, type = "mu"))
summary(ur.kpss(y, type = "tau"))

# Exemple 5.4
summary(ur.kpss(i1, type = "mu"))

summary( ur.kpss(i1, type = "tau"))

require("forecast")
(m1 = Arima(i1, order=c(1,1,0),include.drift=TRUE))

(m2 = Arima(i1, order=c(2,0,0),include.mean=TRUE))
ret= c(3,6,9,12)
t(Box.test.2(residuals(m2),nlag=ret,type="Ljung-Box",fitdf =3))
t_stat(m2)


# 5.4 Significativité illusoire en régression
# Exemple 5.5
data(indbourse)
nikkei =indbourse[,1]
nasdaq =indbourse[,4]

fac0 = as.numeric(nikkei[3]/nasdaq[3])
nas1= fac0*nasdaq

# valan = substr(dimnames(nikkei@.Data)[[1]],1,4)
valan = substr(index(nikkei),1,4)
# lesdates = dimnames(nikkei@.Data)[[1]]
lesdates = index(nikkei)
anunique = unique(valan)
man1 = rep(NA,length(anunique))
for(i  in 1:length(anunique))
{sousdate = lesdates[valan == anunique[i]]
  man1[i] = which(lesdates == sousdate[1]) }

op = par(mfrow=c(2,1), mar=c(3.5,4,1,1),oma=c(0,0,0,0),
 mgp=c(2,.5,0),cex.lab=1.2, cex.axis=1.2)
plot(nikkei,nas1,xlab="nikkei",ylab="nasdaq",pch="+")
par(mar=c(3.5,4,0,1))
# plot.ts(cbind(nikkei,nas1),plot.type="single",
# ylab="nikkei, nasdaq",xlab="temps",xaxt="n",lty=1:2)
plot(cbind(nikkei,nas1),plot.type="single",
 ylab="nikkei, nasdaq",xlab="temps",xaxt="n",lty=1:2)
legend("bottomleft",legend=c("nikkei","nasdaq"),lty=1:2)
axis(1,at=man1, labels= substr(lesdates[man1],1,4))
par(op)


# mod2 = lm(nikkei@.Data~nasdaq@.Data)
mod2 = lm(na.approx(nikkei)~na.approx(nasdaq))
aa=summary(mod2)
aa

op <- par(mfrow=c(2,1),  mar=c(4,4,1,1), oma=c(0,0,0,0))
# plot.ts(residuals(mod2), xlab="temps",ylab="résidu MCO")
plot(residuals(mod2), xlab="temps",ylab="résidu MCO")
abline(h=0)
acf(residuals(mod2),main="",xlab="retard")
par(op)
