########################################################
#######  Codes du chapitre 5
#######  Y. Aragon
########################################################

# Mise à jour effectuée le 17/08/2016


# packages utilisés dans le chapitre
require(caschrono)
require(xtable)
require(dse)
require(polynom)



# 10.1 Identification de la série des résidus obtenus par MCO 
data(khct)
plot.ts(khct,xlab="temps",main="",cex.lab=.9,cex.axis=.8,
 oma.multi=c(4.5,3,.2,0),mar.multi=c(0, 4, 0, .5),las=0)



khct.df<-as.data.frame(window(cbind(khct,time(khct),
(time(khct)-1977)^2), end=c(1983,12)))
colnames(khct.df) <- c("kwh","htdd","cldd","t1","t1.2")

mod2 = lm(sqrt(kwh) ~ htdd + cldd + t1 + t1.2, data=khct.df)
u = ts(residuals(mod2), start=c(1970,1), frequency=12)

acf2y(u,numer=FALSE)



require("forecast")
(modar1=Arima(u,order=c(3,0,1),seasonal=list(order=c(1,0,1)),
  include.mean=FALSE))

llag = seq(6,30,6)
t(Box.test.2(residuals(modar1),llag,type ="Ljung-Box",decim=2,fitdf=6))

t_stat(modar1)

(modar2=Arima(u,order=c(2,0,1),seasonal=list(order= c(1,0,1)),
  include.mean= FALSE))
t_stat(modar2)

(modar3=Arima(u,order=c(1,0,1),seasonal=list(order=c(1,0,1)),
 include.mean=FALSE))

(mod.auto=auto.arima(u,max.p=4,max.q=4,max.P=1,approximation=FALSE))

t(Box.test.2(residuals(mod.auto),llag,type="Ljung-Box",decim=2,fitdf=4))


# 10.2 Estimation du modèle ARMAX 

kwh1rc = window(sqrt(khct[,"kwh"]), end=c(1983,12))
xreg1 = khct.df[ ,c("htdd","cldd","t1","t1.2")]

mdarx1=Arima(kwh1rc,order=c(1,0,1),seasonal=list(order=c(1,0,1)),
 xreg=xreg1)

xreg2 = xreg1[,-4]
(mdarx2=Arima(kwh1rc,order=c(1,0,1),seasonal=list(order=c(1,0,1)),
  xreg=xreg2))

t(Box.test.2(residuals(mdarx2),seq(6,30,6),type="Ljung-Box",
 decim=2,fitdf=8))

t_stat(mdarx2)

(mdarx3c=Arima(kwh1rc,order=c(1,0,0),seasonal=list(order=c(1,0,1)),
  xreg=xreg2))

t(Box.test.2(residuals(mdarx3c),seq(6,30,6),type="Ljung-Box",
 decim=2,fitdf=7))

u.3c=kwh1rc-as.matrix(xreg2)%*%as.matrix(mdarx3c$coef[5:7])-mdarx3c$coef[4]

op=par(oma=rep(0,4))
plot(u.3c,type="l")
abline(h=0)
par(op)


tt = (khct.df$t1 - 1973)^2
mmo = lm(u.3c ~tt)
summary(mmo)

(mod2s = lm(sqrt(kwh) ~ htdd + cldd + t1, data=khct.df))



# 10.3 Estimation d'un modèle à erreur non stationnaire
(mdarx1 = Arima(kwh1rc,order=c(1,0,1),seasonal=list(order=c(1,0,1)),
 xreg = xreg1))

(modarimax1=Arima(kwh1rc,order=c(1,0,1),seasonal=list(order=c(0,1,1)),
  xreg=xreg1))

(modarimax1b=Arima(kwh1rc,order=c(1,0,1),seasonal=list(order=c(0,1,1)),
  xreg = xreg1, include.drift=TRUE) )

xreg1b = xreg1[,-3]
(mx2b=Arima(kwh1rc,order=c(1,0,1),seasonal=list(order=c(0,1,1)),
  xreg = xreg1b, include.drift= TRUE))
t(Box.test.2(residuals(mx2b),seq(6,30,6),type ="Ljung-Box",decim=2,fitdf=7))

t_stat(mx2b)

(modarimax2c=Arima(kwh1rc,order=c(1,0,0),seasonal=list(order=c(0,1,1)),
  xreg=xreg1b,include.drift=TRUE) )
t(Box.test.2(residuals(modarimax2c),seq(6,30,6),type ="Ljung-Box",decim=2,fitdf=6))

t_stat(modarimax2c)


# 10.4 Prévision de l'année 1984
# Prévision par le modèle MCO.
khct.df.84<-as.data.frame(window(cbind(khct,time(khct),
(time(khct)-1977)^2), start = c(1984,1)))
colnames(khct.df.84) <- c("kwh","htdd","cldd","t1","t1.2")
p2=predict(mod2,khct.df.84,interval="prediction",level=0.80,se.fit=TRUE)

# Prévision par le modèle ARMAX.
xreg.p=khct.df.84[,2:4]
prev.3c=forecast(mdarx3c,h=12,level=c(80,95),fan=FALSE,xreg=xreg.p)

str(prev.3c)

etyp.pred  = (prev.3c$upper[,1]-prev.3c$mean)/qnorm(0.9)
etyp.pred2 = (prev.3c$upper[,2]-prev.3c$mean)/qnorm(0.975)

# Prévision par le modèle ARIMAX.
xreg.2c=khct.df.84[,c(2,3,5)]
prev.2c=forecast(modarimax2c, h=12,level=c(80,95),fan=FALSE,xreg=xreg.2c)

kwh2rc <- window(sqrt(khct[,"kwh"]), start = c(1984,1))
op=par(oma=rep(0,4))
aa= seq(as.Date("2000/1/1"), by="month", length.out=12)
id.mois= months(aa, abbreviate = TRUE)
plot(ts(cbind( kwh2rc, prev.3c$lower[,1],prev.3c$upper[,1],p2$fit[,2:3]),frequency=12,start=c(1984,1)),
 plot.type = "single", lty=c(1,2,2,3,3),  ylab=expression(sqrt(kwh)),cex.main=.8, xlab="1984",
 xaxt = "n")
axis(1, at=seq(1984, 1984.917, length.out=12),labels=id.mois)
title(main="Année 1984 - Bandes de prédiction à 80% : modélisations  ARMAX et MCO", cex.main=0.9)
legend( par("usr")[1], par("usr")[4],c("Valeur  observée","Prédiction ARMAX","Prédiction MCO"),
 lwd=1, lty=c(1,2,3))
par(op)



# ARMAX
un80=sum( (kwh2rc < prev.3c$lower[,1])|(kwh2rc > prev.3c$upper[,1]))
un95=sum( (kwh2rc <  prev.3c$lower[,2] ) | (kwh2rc > prev.3c$upper[,2]))
cat("taux de non appartenance 95 (ARMAX)= ",sum(un95)/12,"\n")
cat("taux de non appartenance 80 (ARMAX)= ",sum(un80)/12,"\n")
pp= c(sum(un80),sum(un95))/12
# ARIMAX
un80i=sum((kwh2rc < prev.2c$lower[,1])|(kwh2rc > prev.2c$upper[,1]))
un95i=sum((kwh2rc < prev.2c$lower[,2])|(kwh2rc > prev.2c$upper[,2]))
ppi= c(sum(un80i),sum(un95i))/12
cat("taux de non appartenance 80 (ARIMAX)= ",sum(un80i)/12,"\n")
cat("taux de non appartenance 95 (ARIMAX)= ",sum(un95i)/12,"\n")

# Comparaison des prédictions
# MCO
p0=predict(mod2, khct.df.84, se.fit=TRUE)
# ARMAX
prev.3c=forecast(mdarx3c,h=12,level=c(80,95),fan=FALSE,xreg=xreg.p)
# ARIMAX
prev.2c=forecast(modarimax2c,h=12,level=c(80,95),fan=FALSE,xreg=xreg.2c)
# EQM
b.arimax = cumsum((kwh2rc - prev.2c$mean)^2)/ 1:12
b.armax = cumsum((kwh2rc - prev.3c$mean)^2)/ 1:12
b.mco = cumsum((kwh2rc - p0$fit)^2)/ 1:12
aaa= t(cbind(b.mco,b.armax, b.arimax))
rownames(aaa) = c("MCO","ARMAX","ARIMAX")
colnames(aaa) = id.mois

# Tableau 10.1
xtable(aaa[,c(2,4,6,8,10,12)],caption="Erreurs quadratiques de prévision pour 1984.", label="eqm84")


# 10.5 Prédiction sur la série non transformée

ret.u= rev(u.3c)[1:13]
# 12 valeurs retardées de z_t
ret.z = rev(residuals(mdarx3c))[1:12]
#  modèle du bruit u_t
coef0=mdarx3c$coef

ret.u= rev(u.3c)[1:13] # 13 valeurs retournées
ret.z = rev(residuals(mdarx3c))[1:12] # 12 valeurs retournées de z_t
A.u=polynomial(c(1,-coef0[1]))*polynomial(c(1,rep(0,11),-coef0[2]))
A.arr=array(A.u ,c(length(A.u),1,1))
B.arr=array(c(1,coef0[3]),c(2,1,1))
mod.u=ARMA(A=A.arr, B=B.arr)

pred.moy=mdarx3c$coef[4]+
    as.matrix(xreg.p)%*%as.matrix(mdarx3c$coef[5:7])
nsim=10000
pred.y=matrix(NA,ncol=nsim,nrow=12)
set.seed=539
wsim=matrix(rnorm(12*nsim, sd=mdarx3c$sigma2^.5),ncol=nsim,nrow=12)
for ( i in 1:nsim)
{pred.y[,i]=pred.moy + simulate(mod.u, y0=ret.u,
      noise=list(w=as.matrix(wsim[,i]),w0=ret.z),sampleT=12)$output}
# retour aux données initiales
pred.kwh=pred.y^2
# quantiles
quant=apply(pred.kwh,1,quantile,probs=c(.05,.10,.90,.95,.5))

op=par(oma=rep(0,4))
plot(ts(cbind(as.numeric(khct.df.84[,"kwh"]), t(quant)),frequency=12,start=c(1984,1)),
       plot.type="single",lty=c(1,2,3,3,2,1), ylab= "kwh",cex.main=.8, xlab="1984",
       xaxt = "n",col=c(rep("black",5),"grey"))
axis(1, at=seq(1984, 1984.917, length.out=12),labels=id.mois)
title(main="Année 1984 - Intervalle de prédiction interquantiles à 90% et 80% et série", cex.main=0.9)
legend(par("usr")[1], par("usr")[4],c("Série","Intervalle à 90%","Intervalle à 80%","Médiane"),
         lwd=1, lty=c(1,2,3,1),col=c(rep("black",5),"grey"))
par(op)
