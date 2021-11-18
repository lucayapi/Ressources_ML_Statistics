##############################################
#####  Codes du chapitre 11
#####  Y. Aragon
##############################################

# Mise à jour effectuée le 17/08/2016


# packages utilisés dans le chapitre
require(caschrono)
require(fBasics)
require(xtable)
require(car)

# 11.1 Analyse exploratoire 

lait2=read.table(system.file("/import/collecteLait.txt",package="caschrono"),
 header=FALSE, sep=";",colClasses=c("character",rep("numeric",3)),dec=".",
col.names=c("mois","an","evol","coll.v","cum.v","coll.m","cum.m"))
lait= ts(lait2$coll.v/1000,start=c(1979,1), frequency = 12)
head(lait,3)

# Tendance
op=par(oma=rep(0,4))
decomp=decompose(lait)
plot(decomp)
abline(v=1984)
par(op)


num=which.min(lait)
t.lait = time(lait)
cat("temps de la collecte minimale : ",t.lait[num],"\n")

# Month plot
op=par(oma=rep(0,4))
monthplot(lait, xlab="mois") 
par(op)




# Year plot
ans = 1979:2009
freq=12
y.m=as.matrix(window(lait,start=c(ans[1],1),end=c(ans[1],freq)))
for(i in ans[-1])
{
y.m=cbind(y.m, as.matrix(window(lait,start=c(i,1),end=c(i,freq))))
}
q2=t(apply(y.m[,-(1:5)],1,quantile,c(.25,.75)))
ypl=cbind(y.m[,1:5],q2)
colnames(ypl) = c(unlist(lapply(ans[1:5],toString)),"q.25","q.75")

op=par(oma=rep(0,4))
matplot(1:(freq+1),rbind(ypl,NA),xlab="mois",ylab="séries",
     lty= c(rep(1,5),rep(2,2)), type="o", lwd=2,
     col=c(gray(0:4 / 8),rep("black",2)),pch=c(21:25,NA,NA))
x.t = freq
y.t = ypl[freq,]
y.t = ypl[freq,]*c(.97,.98,1.011,.99,1.011,1.009,1.014)
text(x.t,y.t,colnames(ypl),pos=4,cex=.8,col=c(gray(0:4 / 8),rep("black",2)))
par(op)



# Lag plot
op=par(oma=rep(0,4))
lag.plot(rev(lait),set=c(1:12),pch="+",col="black")
par(op)



# Normalité
aa= dagoTest(lait)
a1 = aa@test$p.value[1]
a2 = dagoTest(log(lait))@test$p.value[1]
a3 = dagoTest(lait^.5)@test$p.value[1]
aa= as.matrix(c(a1,a2,a3))
rownames(aa) = c("Série brute","log","racine")
colnames(aa) ="p-value"

# Tableau 11.1
xtable(aa, caption="Lait - Test de normalité de D'Agostino"  , label="laitnor", digits=6)

(ptr=powerTransform(lait)$lambda)
lait.tr =(lait^ptr - 1)/ptr
dagoTest(lait.tr)@test$p.value[1]


# 11.2 Modélisation avant 1984
log.lait=log(lait)
lait.avant=window(log.lait,end=c(1983,12))

op=par(oma=rep(0,4))
xy.acfb(lait.avant,lag.max=26,numer=FALSE)
par(op)




require("forecast")
(mod0=auto.arima(lait.avant))

ret=c(6,12,18,24,30)
t(Box.test.2(residuals(mod0),ret,type ="Ljung-Box",decim=2,fitdf=3))
t_stat(mod0)

(mod1b=Arima(lait.avant, order=c(1,0,5),seasonal=list(order=c(1,0,0)),
  fixed=c(NA,rep(0,4),rep(NA,3)),method="CSS"))
t(Box.test.2(residuals(mod1b),ret,type ="Ljung-Box",decim=2,fitdf=4))
t_stat(mod1b)

aa= dagoTest(residuals(mod1b))

(mod1bm=Arima(lait.avant, order=c(1,0,5),seasonal=list(order=c(1,0,0)),
   fixed=c(NA,rep(0,4),rep(NA,3))))
t(Box.test.2(residuals(mod1bm),ret,type ="Ljung-Box",decim=2,fitdf=4))
t_stat(mod1bm)

cmod1=mod1bm$coef[mod1bm$mask==TRUE]

demi.b = qnorm(.9)*mod1bm$sigma2^.5
b.inf =  fitted(mod1bm) - demi.b
b.sup = fitted(mod1bm) + demi.b
dans = (b.inf < lait.avant) & (lait.avant < b.sup)

op=par(oma=rep(0,4))
plot.ts(cbind(lait.avant,b.inf,b.sup),plot.type="single", ylab="Log collecte",
xlab="temps",lty=c(1,2,2))
title(main="Log collecte de lait et bande de confiance à 80%")
par(op)



pred56 = forecast(mod1bm,h=12, level=80)
et.pred = (pred56$upper-pred56$mean)/qnorm(.9)


# 11.3 Essai de modélisation après l'introduction des quota
apres95 = window(log.lait,  start = c(1995,1))
(mod2=Arima(apres95,order=c(1,0,0),seasonal=list(order=c(1,0,0))))
t(Box.test.2(residuals(mod2),ret,type="Ljung-Box",decim=2,fitdf=2))

(mod2=Arima(apres95,order=c(1,0,0),seasonal=list(order=c(1,1,0)),
 include.drift=FALSE))
t(Box.test.2(residuals(mod2),ret,type="Ljung-Box",decim=2,fitdf=2))
t_stat(mod2,decim=2)


# 11.4 Modélisation SARIMA de toute la série 
mod2.tot = Arima(log.lait,order=c(1,0,0),seasonal=list(order=c(1,1,0)),
  include.drift=FALSE)
summary(mod2.tot)
t(Box.test.2(residuals(mod2.tot),ret,type ="Ljung-Box",decim=2,fitdf=2))
t_stat(mod2.tot,decim=2)



mod3.tot=Arima(log.lait,order=c(1,0,0),
  seasonal=list(order=c(2,1,0)),include.drift=FALSE)
summary(mod3.tot)
t(Box.test.2(residuals(mod3.tot),ret,type ="Ljung-Box",decim=2,fitdf=3))
t_stat(mod3.tot,decim=2)

acf2y(residuals(mod3.tot),numer=FALSE)



# 11.5 Modélisation ARMAX de la collecte 
# 11.5.1 Modélisation MCO 
f = t(as.matrix(1:6))/12
num.temps = as.matrix(1:length(lait))
temps = time(log.lait)
ind.av = ifelse(as.numeric(temps)<1984,1,0)        
ind.apr = 1-ind.av
xmat0 = cbind(cos(2*pi*num.temps%*%f), sin(2*pi*temps%*%f))[,-12]
xmat0 = as.data.frame(xmat0)
xmat0b = cbind(ind.av,ind.apr,xmat0)
colnames(xmat0b) = c("ind.av","ind.apr","cos1","cos2","cos3","cos4","cos5",
                       "cos6","sin1","sin2","sin3","sin4","sin5")

mod1 = lm(log.lait~ind.av+ind.apr+cos1+ cos2+ cos3+cos4+cos5+cos6+sin1+sin2+sin3+sin4+sin5 -1, data=xmat0b)
summary(mod1)
resid.mod1= ts(residuals(mod1),start=c(1979,1),frequency=12)

# 11.5.2 Identification des résidus de l'ajustement MCO 
op=par(oma=rep(0,4))
acf2y(resid.mod1,numer=FALSE)
par(op)




(modar1=Arima(resid.mod1,order=c(1,0,0),seasonal=list(order=c(2,0,0))))
t(Box.test.2(residuals(modar1),ret,fitdf=3,decim=4))

op=par(oma=rep(0,4))
acf2y(residuals(modar1),numer=FALSE)
par(op)



modar1b=Arima(resid.mod1,order=c(9,0,0),seasonal=list(order =c(2,0,0)))
t(Box.test.2(residuals(modar1b),ret,fitdf=11,decim=4))
t_stat(modar1b)

(modar1c=Arima(resid.mod1,order=c(9,0,0),seasonal=list(order=c(2,0,0)),
 fixed =c(NA,rep(0,6),rep(NA,4)),include.mean=FALSE ))
t(Box.test.2(residuals(modar1c),ret,fitdf=5,decim=4))



# 11.5.3 Modélisation simultanée de la moyenne et de l'erreur 
(modar1X=Arima(log.lait,order=c(9,0,0),seasonal=list(order=c(2,0,0)),
include.mean=FALSE,xreg=xmat0b,fixed=c(NA,rep(0,6),rep(NA,4),rep(NA,13))))
t(Box.test.2(residuals(modar1X),ret,fitdf=18,decim=4))

t_stat(modar1X)

xmat0c = xmat0b[,c(1:3,12) ]
(modar2X=Arima(log.lait,order=c(9,0,0),seasonal=list(order=c(2,0,0)),
 include.mean=FALSE,xreg=xmat0c,fixed=c(NA,rep(0,6),rep(NA,8))))
t(Box.test.2(residuals(modar2X),ret,fitdf=9,decim=4))
t_stat(modar2X)

aa= modar2X$var.coef
(aa0 = aa[6:7,6:7])
difa = as.matrix(c(1,-1))
(vardif = t(difa)%*%aa0%*%difa )
(t.stat= (t(difa)%*%modar2X$coef[6:7])/vardif^.5)

(modar3X=Arima(log.lait,order=c(9,0,0),seasonal=list(order=c(2,0,0)),
 xreg=xmat0c[,-(1:2)],fixed=c(NA,rep(0,6),rep(NA,7))))
t(Box.test.2(residuals(modar3X),ret,fitdf=8,decim=4))
round(t_stat(modar3X),digits=2)
resid.x = residuals(modar3X)

aa= dagoTest(resid.x)

cf3x = round(modar3X$coef[modar3X$mask],digits=2)
noms = names(cf3x)