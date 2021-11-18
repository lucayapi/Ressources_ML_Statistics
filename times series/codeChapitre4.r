###########################################################
######   Codes du chapitre 4
######   Y. Aragon
##########################################################

# Mise à jour effectuée le 17/08/2016

# packages utilisés dans le chapitre
require(caschrono)
require(FitARMA)
require(xtable)
require(polynom)
require(nlme)
require(dynlm)

# 4.1 Stationnarité 
# 4.1.1 Fonction d'autocorrélation d'une série stationnaire
# 4.1.2 Bruit blanc
# Exemple 4.1
set.seed(123)
y1 = arima.sim(n=100,list(ar=-.7), sd=sqrt(4))
y2 = arima.sim(n=100,list(ar=c(rep(0,11),-.7)), sd=sqrt(4))
ret= c(3, 6, 9, 12)
a1=Box.test.2(y1, nlag=ret,type="Ljung-Box",decim=2)
a2=Box.test.2(y2, nlag=ret, type="Ljung-Box",decim=2)
a12 = cbind(a1,a2[,2])
colnames(a12)= c("Retard","p-val. y1","p-val. y2")
a12

# 4.2 Série linéaire 
# 4.3 Fonctions d'autocorrélation
# 4.3.1 Fonction d'autocorrélation d'un AR 
# 4.3.2 Fonction d'autocorrélation d'un MA
# Exemple 4.3
require(FitARMA)
phi0=c(-.7, .2) ; theta0=0.6
g = TacvfARMA(phi = phi0, theta=-theta0, lag.max=5)
g/g[1]

ARMAacf(ar=phi0, ma = theta0, lag.max = 5)

# Exemple 4.4
set.seed(219)
y3 = arima.sim(n=200, list(ma=c(-.3,.6)),sd = sqrt(1.5))

mod0=arima(y3,order=c(0, 0, 2))

acf(residuals(mod0),xlab="retard",main="")


# 4.4 Prévision 
# 4.4.1 Principe 
# 4.4.2 Fonction d'autocorrélation partielle
# Exemple 4.5
g = TacvfARMA(phi=c(-.7,.2),lag.max=4)
(a=PacfDL(g/g[1], LinearPredictor = TRUE))
g[1] - t(as.matrix(a$ARCoefficients))%*%as.matrix(g[-1])

# 4.4.3 Prévision d'un modèle autorégressif
# 4.4.4 Prévision d'un MA(q) 
# 4.5 Estimation
# 4.5.1 Exemples
# Exemple 4.1
require("forecast")
(my1 = Arima(y1, order=c(1,0,0),include.mean = FALSE))

ret=c(3,6,9,12)
aa=Box.test.2(residuals(my1),nlag=ret,type="Ljung-Box",decim=4,fitdf=1)
colnames(aa)= c("Retard","p-val.")
t(aa)

t_stat(my1)

# Exemple 4.2
(my2a = Arima(y2, order=c(12,0,0),include.mean=FALSE))

(my2a = Arima(y2, order=c(12,0,0),include.mean=FALSE,
 fixed=c(rep(0,11),NA)))

(my2a=Arima(y2,include.mean=FALSE,
 seasonal=list(order=c(1,0,0),period=12)))



# 4.5.2 Modèle ARMA saisonnier (modèle SARMA)
# Exemple 4.7
set.seed(7392)
(autopol = polynomial(c(1,0.8))*polynomial(c(1,0,0,0,-0.7)))
yd=arima.sim(n=200,list(ar=-autopol[-1],ma=c(0,0.6)),sd=sqrt(1.5))
yd=yd+4
acf.th=ARMAacf(ar=-autopol[-1],ma=c(0,0.6),lag.max=20,pacf=FALSE)
pacf.th=ARMAacf(ar=-autopol[-1],ma=c(0,0.6),lag.max=20,pacf=TRUE)

plotacfthemp(yd,ar=-autopol[-1], ma=c(0,0.6), lag.max=20)


op=par(oma=rep(0,4))
lag.plot(rev(yd),9,layout=c(3,3),ask=NA,do.lines=FALSE,main="SARMA yd",diag.col="red")
par(op)

(msarma=Arima(yd,order=c(1,0,2),
seasonal=list(order=c(1,0,0),period=4)))

(msarma=Arima(yd,order=c(1,0,2),
  seasonal=list(order=c(1,0,0),period=4),
 fixed=c(NA,0,NA,NA,NA)))

# 4.5.3 Modèle ARMAX
# Exemple 4.8
temps = time(LakeHuron)
mod1.lac=Arima(LakeHuron,order=c(1,0,0),xreg=temps,method="ML")
(cf.arima = mod1.lac$coef)
ychap.arima = fitted(mod1.lac)
resi.arima = residuals(mod1.lac)
nret = c(3,6,9,12)
bl1=Box.test.2(resi.arima,nlag=nret,type="Ljung-Box",decim=4,fitdf=2)

tu = 1:length(temps)
mod2.lac=gls(LakeHuron~temps,correlation=corAR1(form=~tu),method="ML")
(cf.gls = mod2.lac$coef)
ychap.gls = fitted(mod2.lac)
resi.gls = residuals(mod2.lac)
bl2=Box.test.2(resi.gls,nlag=nret,type="Ljung-Box",decim=4,fitdf=2)

mod3.lac = dynlm(LakeHuron~ 1+temps+L(LakeHuron, 1))
cf.dynlm= mod3.lac$coef
ychap.dynlm = fitted(mod3.lac)
resi.dynlm = residuals(mod3.lac)
bl3=Box.test.2(resi.dynlm,nlag=nret,type="Ljung-Box",decim=4,fitdf=2)

c2= cf.arima[1] ; c1 = cf.arima[3] *(1- cf.arima[1])
c0 = cf.arima[2]*(1 - cf.arima[1])+ cf.arima[3]*(1875 -1 -cf.arima[1]*1875+ 2* cf.arima[1])
coefx= round(cbind(as.matrix(cf.dynlm),c(c0,c1,c2)),digits=5)
colnames(coefx) = c("dynlm","Arima")
rownames(coefx)= c("Intercept","temps","phi")

bl123 = cbind(bl1,bl2[,2],bl3[,2])
colnames(bl123 )= c("retard","Arima","gls","dynlm")
y123=cbind(ychap.arima,ychap.gls,ychap.dynlm)
r123=cbind(round((1:length(ychap.arima)),digits=0),resi.arima, resi.gls,resi.dynlm )[1:5,]
colnames(r123)=c("t","resi.arima","resi.gls","resi.dynlm")

require(xtable)
xtable(coefx, caption="Estimation du modèle avec variable dépendante retardée, (1) directe par   \texttt{dynlm()} et (2)
via   \texttt{Arima()}.",label="compar4", digits=4)

xtable(r123, caption="Résidus, suivant les fonctions utilisées.",label="compar3", digits=4)

xtable(bl123, caption="Test de blancheur des résidus, suivant les fonctions utilisées.",label="compar2", digits=4)


# 4.6 Construction d'un ARMA ou d'un SARMA
# 4.6.1 Identifcation d'un ARMA
# Exemple 4.8
set.seed(4123)
n2 = 210
yc = arima.sim(n = 200, list( ar = -0.8, ma= c(-0.3, 0.6)), sd = sqrt(1.5))
yc = yc-10

armaselect(yc, nbmod=5)

# 4.6.2 La méthode MINIC
# 4.7 Exercices 