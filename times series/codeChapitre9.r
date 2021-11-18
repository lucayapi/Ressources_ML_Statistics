#########################################
######  Codes du chapitre 9
######  Y. Aragon
#######################################

# Mise à jour effectuée le 17/08/2016


# packages utilisés dans le chapitre
require(caschrono)
require(TSA)
require(xtable)


# 9.1 Exploration
data(nottem)
op=par(oma=rep(0,4))
plot.ts(nottem, xlab="temps",ylab="température")
par(op)



# 9.2 Modélisation 
nott1 = window(nottem, end=c(1936,12))
nott2 = window(nottem, start=c(1937,1))

op=par(oma=rep(0,4))
plot2acf(nott1, diff(nott1,12), lag.max=40, main=c("nott1",expression(paste("(1-",B^{12},") nott1",sep=""))))
par(op)


# 9.2.1 Modèle SARIMA 
require("forecast")
fitm = Arima(nott1,order=c(1,0,0), list(order=c(2,1,0), period=12))
summary(fitm)

ret6 = seq(6,30,6)
t(Box.test.2(residuals(fitm),ret6,type="Ljung-Box",decim=2,fitdf=3))

ccfm = round(fitm$coef,digits=3)


# 9.2.2 Régression sur fonctions trigonométriques
f = t(as.matrix(1:6))/12
temps = as.matrix(1:length(nottem))
xmat0 = cbind(cos(2*pi*temps%*%f), sin(2*pi*temps%*%f))[,-12]
xmat0 = as.data.frame(xmat0)
colnames(xmat0) = c("cos_1","cos_2","cos_3","cos_4","cos_5","cos_6",
                    "sin_1","sin_2","sin_3","sin_4","sin_5")

xmat1 = xmat0[1:length(nott1),]
xmat2 = xmat0[(length(nott1)+1):length(nottem),]
attach(xmat1,warn.conflicts = FALSE)
# on régresse  sur les colonnes de la matrice
mod1 = lm(nott1~cos_1+ cos_2+ cos_3+cos_4+cos_5+cos_6+sin_1+
                  sin_2+sin_3+sin_4+sin_5, data=xmat1)
summary(mod1)

mod2 = lm(nott1~cos_1+ sin_1+ sin_2+sin_4, data=xmat1)
summary(mod2)

var.mod2= summary(mod2)$sigma^2

md2sm= round(mod2$coefficients,digits=4)
s2simp=sum((fitted(mod2)-mean(nott1))^2)

attach(xmat1,warn.conflicts = FALSE)
xmat1a = cbind(cos_1, sin_1,sin_2,sin_4)
mod2b = Arima(nott1,order=c(0,0,0), xreg =  xmat1a)

attach(xmat2,warn.conflicts = FALSE)
xmat2a = cbind(cos_1, sin_1,sin_2,sin_4)
pred.mco2 = forecast(mod2b,  xreg=xmat2a)

op=par(oma=rep(0,4))
acf2y(residuals(mod2),lag.max=30,numer=FALSE)
par(op)


modar1b=Arima(residuals(mod2),order=c(1,0,0),
 list(order=c(2,0,0),period=12),include.mean=FALSE)
summary(modar1b)

t(Box.test.2(residuals(modar1b),ret6,type="Ljung-Box",decim=2,fitdf=3))

attach(xmat1,warn.conflicts = FALSE)
mod3b=Arima(nott1,order=c(1,0,0),list(order=c(2,0,0),period=12),
              xreg =cbind(cos_1,sin_1,sin_2,sin_4))
summary(mod3b)

m3bcf = round(coef(mod3b),digits=4)

t(Box.test.2(residuals(mod3b),ret6,type="Ljung-Box",decim=2,fitdf=8))


sinit = round(sum((nott1-mean(nott1))^2),digits=0)
smco = round(var(residuals(mod2b)),digits=4)
sarmax= round(var(residuals(mod3b)),digits=4)



# 9.3 Prévision
pred.sarima = forecast(fitm, h=36)
tnupp = ts(pred.sarima$upper[,1], start= c(1937,1), frequency=12)
tnlow = ts(pred.sarima$lower[,1], start= c(1937,1), frequency=12)

attach(xmat2,warn.conflicts = FALSE)
pred.armax=forecast(mod3b,h=36,xreg=cbind(cos_1,sin_1,sin_2,sin_4))

pred.mco=predict(mod2,xmat2,se.fit=TRUE,level=0.80,
  interval="prediction")



# 9.4 Comparaison 
mxupp = ts(pred.armax$upper[,1], start= c(1937,1), frequency=12)
mxlow = ts(pred.armax$lower[,1], start= c(1937,1), frequency=12)
op=par(oma=rep(0,4))
plot.ts(cbind(nott2, mxupp, mxlow, tnupp,tnlow), plot.type = "single",axes=FALSE,
  col=c(1,1,1,1,1), lty=c(1,3,3,5,5),xlab= "temps", ylab="température",
  main="Température, prévisions et réalisations",cex.main=.8, cex= .7,cex.lab=.7)
aa= time(nott2)[36]
# examiner le résultat de ce qui précède avant exécution
axis(1,at = c(1937:1939,aa), lab= c(as.character(1937:1939),""))
axis(2, pretty(nott2),las=1)
# examiner le résultat de ce qui précède ...
topleft= par()$usr[c(1,4)]  # coordonnées du coin supérieur gauche du graphique
legend(1937.8,topleft[2], c("observé", "ARMAX 80%", "ARMAX 80%", "SARIMA 80%","SARIMA 80%"),
       lty = c(1,3,3,5,5),  merge = TRUE, cex=.7)
par(op)


aa=cbind(as.matrix(nott2),pred.mco2$mean,
         pred.armax$mean,pred.sarima$mean)
eq = (aa[,1]%*%matrix(1,nrow=1,ncol=3) - aa[,c(2,3,4)])^2
colnames(eq) = c("mco", "armax","sarima")
eqm = apply(eq, 2, "cumsum")/ 1:36
colnames(eqm) = c("mco", "armax","sarima")

op=par(oma=rep(0,4))
par(fig= c(0,1,0.5,1), mar=c(0,5,3,2), mgp=c(2.5,1,0))
matplot(1:36,eqm, type="l",lty=1:3, col="black",xlab="",xaxt="n",lwd=2)
legend(x=25, y=3,c("mco","armax","sarima"),lty=1:3)
par(new=T,fig= c(0,1,0,0.5),mar=c(5,5,0,2), mgp=c(2.5,1,0))
matplot(1:36,eq,type="l",lty=1:3, col="black",xlab="horizon",lwd=2)
legend(x=25, y=25, c("mco","armax","sarima"), lty=1:3)
par(op)



# 9.5 Analyse spectrale 
op=par(oma=rep(0,4))
periodogram(nottem)
par(op)


aa= periodogram(nottem,plot=FALSE)
aa$freq[which.max(as.vector(aa$spec))]

ab =  order(-aa$spec)[1:5]
frq.spe= rbind(aa$freq[ab],aa$spec[ab])
rownames(frq.spe)= c("Fréquence","Périodogramme")

xtable(frq.spe,caption="Fréquences et périodogramme de plus grande énergie.", 
 label="tabl.freq")