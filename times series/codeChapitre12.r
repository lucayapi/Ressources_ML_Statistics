###############################################
########  Code du chapitre 12
########  Y. Aragon
###############################################

# Mise à jour effectuée le 17/08/2016

# packages utilisés dans le chapitre
require(caschrono)
require(fBasics)
require(xtable)
require(fGarch)
require(FinTS)

# Introduction
# 12.1 Notions de base 
data(csdl)
aa = returns(csdl, percentage = TRUE)
aab = aa[complete.cases(aa) == TRUE,]
# r.csdl = its(aab, as.POSIXct(row.names(aab)))
r.csdl = zoo(aab, as.POSIXct(row.names(aab)))
aa= dagoTest(r.csdl[,"Danone"], title = "Danone", description = NULL)
res.aa = cbind(aa@test$statistic,aa@test$p.value)

colnames(res.aa)= c("Stat. de test","p-value")

# Tableau 12.1
xtable(res.aa, caption="Rendement Danone -Test de normalité de D'Agostino"  , label="danonor", digits=4)

density.plot=function(x,legende=FALSE,...)
{
H<-hist(x, sub = NULL, ylab = "densité", freq=FALSE, ...)
abline(v=0,lwd=2)
rug(x, ticksize = 0.01)
xmin=par()$usr[1];xmax=par()$usr[2]
tab<-seq(xmin,xmax,0.1)
lines(tab, dnorm(tab, mean(x),sd(x)), col="red", lty=2, lwd=2)
lines(density(x),lwd=2, col="orange")
if(legende)
 {lg0 =c("estimation n.p. de la densité","estimation d'une   gaussienne")
 legend("topright",legend=lg0,lty=c(1,2),lwd=2,
 col=c("orange","red"),cex=0.9)
 }
}
density.plot(r.csdl[,3], xlab = "rendement", ylim=c(0,0.3),
main="Danone", nclass=12, xlim=c(-10,10), legende=TRUE )


# 12.2 Modèles d'hétéroscédasticité conditionnelle
# 12.3 ARCH(1) et test d'hétéroscédasticité conditionnelle
# 12.3.1 Simulation d'un ARCH(1) 
# Exemple 12.1
spec.1=garchSpec(model=list(mu=5,omega=0.1,alpha=0.9,beta=0),
   rseed=397)
archsim.1 = garchSim(extended = TRUE,spec.1, n = 300, n.start=10)
head(archsim.1,2)


# 12.3.2 Moments d'un ARCH(1)
# 12.3.3 Tests d'hétéroscédasticité conditionnelle 
ArchTest(archsim.1[,1], lag=12)
ArchTest(archsim.1[,3], lag=12)

# 12.4 Estimation et diagnostic d'ajustement d'un GARCH 
# Exemple 12.1 (Simulation et estimation d'un GARCH(2,1))
spec=garchSpec(model=list(mu=2,omega=0.09,alpha =c(0.15, 0.3),
  beta = 0.4), rseed=9647)
var.margi = 0.09/(1 - 0.15 - 0.3-0.4)
y = garchSim(spec, n = 420, extended = TRUE)
y1 = y[21:420,1]

op<-par(oma=rep(0,4),cex.lab=1.2)
plot.ts(y1, xlab="temps")
par(op)



m1 = mean(y1)
(q5 = quantile(abs(y1-m1), probs=c(.975,.98,.985,.99,.995)))
extrem985 = which(abs(y1-m1) > q5[3])
cat("nombre de points : ",length(extrem985),"\n")
y1b = y1
y1b[extrem985]  = q5[3]* sign(y1[extrem985]-m1)

op<-par(oma=rep(0,4),cex.lab=1.2)
plot.ts(y1b, xlab="temps")
title("Série corrigée des points extrêmes")
par(op)



mod1=garchFit(~garch(2,1),data=y1b,trace=FALSE,include.mean=TRUE)
summary(mod1)

var.marg.est<-function(mod)
{param.estim = mod@fit$par
 std.estim = mod@fit$se.coef
 k<-which(names(param.estim)=="omega")
 value = param.estim[k]/(1-sum(param.estim[(k+1):length(param.estim)]))
cat("variance marginale : ",value,"\n")
}
var.marg.est(mod1)

op<-par(oma=rep(0,4),cex.lab=1.2)
mu= coef(mod1)["mu"]
sigma.t=mod1@sigma.t
ymat=cbind(y1b, mu-2*sigma.t,mu+2*sigma.t)
matplot(1:nrow(ymat),ymat,type="l",col="black",ylab="",xlab="temps",
           lwd=1.5,lty=c(1,2,2))
title("GARCH(2,1) simulé +/- 2 écarts-types conditionnels",cex=0.8)
par(op)



prop.obs= mean((y1b > mu-2*sigma.t) &
               (y1b < mu+2*sigma.t))
prop.att = 1 - 2*(1 - pnorm(2))

mod1fx=garchFit(~garch(1,0),data=y1b,trace=FALSE,include.mean=TRUE)
summary(mod1fx)


# 12.5 Prévision
p.cond=predict(mod1, mse="cond", n.ahead=20, plot=FALSE)
p.cond[c(1,2,19,20),]

b.inf= p.cond$meanForecast -1.96*p.cond$standardDeviation
b.sup= p.cond$meanForecast +1.96*p.cond$standardDeviation
matpr= cbind(y[401:420,1],b.inf,b.sup)
matplot(1:20,matpr,type="l", lty=c(1,2,2))



# 12.6 Modèles à erreur conditionnellement hétéroscédastique 
# 12.7 Etude de séries du marché parisien autour de la crise de 2007-2008 

# 12.7.1 Exploration 
op<-par(oma=rep(0,4),cex.lab=1.2,mar=c(4,4,4,2))
max.csdl = apply(csdl,2,which.max)
#datelor = csdl@dates[max.csdl["L_Oreal"]]
#datemax = csdl@dates[max.csdl["Cac40"]]
datelor = index(csdl)[max.csdl["L_Oreal"]]
datemax = index(csdl)[max.csdl["Cac40"]]
#zz<-zoo(as.data.frame(csdl),csdl@dates)
zz <- csdl
opar <- par(mai = c(.8, .8, .2, .8))
plot(zz[,"Cac40"], type = "l", xlab = "année", ylab = "indice Cac 40")
par(new = TRUE)
plot(zz[,2:4], type = "l", ann = FALSE, yaxt = "n", lty=2,
col = rainbow(3), plot.type="single")
abline(v=datemax)
abline(v=datelor)
axis(side = 4)
mtext(side=4,"indices Société générale, Danone, L'Oréal",line=2.5)
legend(x = "topright", bty = "n", lty = c(1,2,2,2),
col = c("black", rainbow(3)),
legend = paste(colnames(zz), c("(échelle de gauche)", rep("(échelle de droite)",3))))
usr <- par("usr")
par(op)


# 12.7.2 Etude des rendements 
#rendav.06 = rangeIts(r.csdl, end= "2007-06-01")
#rendapr.06 = rangeIts(r.csdl, start= "2007-06-02")
rendav.06 = window(r.csdl, end= "2007-06-01")
rendapr.06 = window(r.csdl, start= "2007-06-02")
sk.av06   =  apply(rendav.06,2,skewness, na.rm=TRUE)
kurt.av06 =  apply(rendav.06,2,kurtosis, na.rm=TRUE)
sk.apr06  =  apply(rendapr.06,2,skewness, na.rm=TRUE)
kurt.apr06=  apply(rendapr.06,2,kurtosis, na.rm=TRUE)
sk06 =  rbind(sk.av06,  sk.apr06, kurt.av06 , kurt.apr06)[,2:3]
colnames(sk06) =  c("Socgen", "Danone")
rownames(sk06) = c("asym.av","asym.apr","aplat.av","aplat.apr")


# Tableau 12.3
sk = cbind(sk06)
colnames(sk) = c("Socgen", "Danone")
xtable(sk,caption="Asymétrie et aplatissement des rendements",label="anorm3")


# 12.7.3 Hétéroscédasticité conditionnelle des rendements
aa.av= apply(rendav.06,2, ArchTest, lag=20)

# 12.8 Etude du rendement de L'Oréal
#r.lor <- rangeIts(r.csdl[,"L_Oreal"], start= "2007-12-28")
r.lor <- window(r.csdl[,"L_Oreal"], start= "2007-12-28")
r.lor.0 = r.lor[1:(length(r.lor)-51)]
r.lor.1 = r.lor[(length(r.lor)-50):length(r.lor)]

# 12.8.1 Estimation
# Modélisation du rendement
op<-par(oma=rep(0.5,4))
# xy.acfb(r.lor.0, numer=FALSE)
xy.acfb(as.timeSeries(r.lor.0), numer=FALSE)
par(op)



require("forecast")
#mod.r.lor=Arima(r.lor.0@.Data,order=c(0,0,4),include.mean=FALSE,
# fixed=c(NA,NA,0,NA))
mod.r.lor=Arima(as.numeric(r.lor.0),order=c(0,0,4),include.mean=FALSE,
 fixed=c(NA,NA,0,NA))
summary(mod.r.lor)
rr.lor = residuals(mod.r.lor)
ret=c(6,12,18,24)
t(Box.test.2(rr.lor,ret,type="Box-Pierce",fitdf=3,decim=4))
t_stat(mod.r.lor, decim=4)

ArchTest(rr.lor, lag=20)

op<-par(oma=rep(0.5,4))
density.plot(rr.lor, main="", xlab = "", ylim=c(0,0.2),
nclass=12, xlim=c(-12,12))
par(op)




# Modélisation de l'hétéroscédasticité du rendement
moda=garchFit(~garch(1,1),data=rr.lor,trace=FALSE,
   include.mean=TRUE,na.action=na.pass)
summary(moda)
var.marg.est(moda)


# Modélisation simultanée du rendement et de son hétéroscédasticité
#mod.lor=garchFit(formula=~arma(0,4)+garch(1,1),data=r.lor.0@.Data,
# trace=FALSE,include.mean=FALSE)
mod.lor=garchFit(formula=~arma(0,4)+garch(1,1),data=r.lor.0,
 trace=FALSE,include.mean=FALSE)
summary(mod.lor)

var.marg.est(mod.lor)

param.estim=round(coef(mod.lor),digits=4)

dans = (r.lor.0 < qnorm(.9)*mod.lor@sigma.t)& (-qnorm(.9)*mod.lor@sigma.t < r.lor.0 )
prop.lor=sum(dans)/length(dans)

op<-par(oma=rep(0.5,4))
n1 = length(r.lor.0)-99; n2 = length(r.lor.0);n1.n2=n1:n2;
mat.est=cbind(r.lor.0[n1.n2],qnorm(.9)*mod.lor@sigma.t[n1.n2],
              -qnorm(.9)*mod.lor@sigma.t[n1.n2])
matplot(n1:n2,mat.est,type="l",col="black",lty =c(1,2,2),
  xlab="100 dernières observations",ylab="rendement",xaxt="n")
# axis(1,at=axTicks(1),labels=r.lor.0@dates[axTicks(1)])
axis(1,at=axTicks(1),labels=index(r.lor.0)[axTicks(1)])
par(op)




# 12.8.2 Prédiction du rendement
npred = 50
pred.lor=predict(mod.lor,n.ahead=npred,plot=FALSE,nx=0)
dem.garch = qnorm(0.9)*pred.lor$standardDeviation
binf.garch =  pred.lor$meanForecast-dem.garch
bsup.garch =  pred.lor$meanForecast+dem.garch

pred.arima=predict(mod.r.lor,n.ahead=npred)
str(pred.arima)
demi = qnorm(0.9)*pred.arima$se
binf.arima =  pred.arima$pred-demi
bsup.arima =  pred.arima$pred+demi

op<-par(oma=rep(0.5,4))
mat.p = cbind(bsup.garch,
 binf.garch,binf.arima,bsup.arima,r.lor.1[1:npred])
matplot(1:npred,mat.p,type="l", col="black", lty=c(1,1,2,2,3),
  lwd=2,xlab="prévision à l'horizon 50", ylab="rendement")
leg.txt = c("ARMA-GARCH","ARMA","réalisation")
legend(14,3, leg.txt, lty=c(1,2,3))
par(op)


dans1 = (pred.arima$pred-demi < r.lor.1[1:npred])& ( r.lor.1[1:npred] < pred.arima$pred+demi)
plor1 = sum(dans1)/length(dans1)
dans2 = (binf.garch < r.lor.1[1:npred])& ( r.lor.1[1:npred] < bsup.garch)
plor2 = sum(dans2)/length(dans2)


# 12.9 Etude du rendement de Danone
# r.dan = rangeIts(r.csdl[,"Danone"], start= "2007-12-28")
r.dan = window(r.csdl[,"Danone"], start= "2007-12-28")
str(r.dan)
r.dan.0 = r.dan[1:(length(r.dan)-50)]
r.dan.1 = r.dan[(length(r.dan)-49):length(r.dan)]


# 12.9.1 Modélisation
# Modélisation du rendement
# xy.acfb(r.dan.0,numer=FALSE)
xy.acfb(as.timeSeries(r.dan.0),numer=FALSE)


#mod3=Arima(r.dan.0@.Data, order=c(2,0,0),include.mean= FALSE)
mod3=Arima(as.numeric(r.dan.0), order=c(2,0,0),include.mean= FALSE)
summary(mod3)
rr.dan = residuals(mod3)
t(Box.test.2(rr.dan,c(6,12,18),type="Box-Pierce",decim=8,fitdf=2))
t_stat(mod3)

ArchTest(rr.dan, lag=20)

# Modélisation de l' hétéroscédasticité du rendement
modarch=garchFit(~garch(1,1),data=rr.dan,trace=FALSE,include.mean=FALSE,
 na.action=na.pass)
summary(modarch)

var.marg.est(modarch)

# Modélisation simultanée du rendement et de son hétéroscédasticité
#mod.dan=garchFit(~arma(2,0)+garch(1,1),data=r.dan.0@.Data,trace=FALSE,
# include.mean=FALSE,na.action=na.pass)
mod.dan=garchFit(~arma(2,0)+garch(1,1),data=r.dan.0,trace=FALSE,
 include.mean=FALSE,na.action=na.pass)
summary(mod.dan)

coef.dan = round(mod.dan@fit$par,digits=4)

dans = (r.dan.0 < qnorm(.9)*mod.dan@sigma.t)& (-qnorm(.9)*mod.dan@sigma.t < r.dan.0 )
prop.dan.0=sum(dans)/length(dans)


# 12.9.2 Prédiction du rendement
npred = 50
pred.dan=predict(mod.dan,n.ahead=npred,plot=FALSE,nx=0)
dem.garch = qnorm(0.9)*pred.dan$standardDeviation
binf.garch =  pred.dan$meanForecast-dem.garch
bsup.garch =  pred.dan$meanForecast+dem.garch

pred.dan.ar =predict(mod3,n.ahead=npred,se.fit = TRUE)
dem.ar = qnorm(0.9)*pred.dan.ar$se
binf.ar =  pred.dan.ar$pred-dem.ar
bsup.ar = pred.dan.ar$pred+dem.ar

dans1 = (binf.ar  < r.lor.1[1:npred])& ( r.lor.1[1:npred] < bsup.ar)
pdn1 = sum(dans1)/length(dans1) # ARMA
dans2 = (binf.garch < r.lor.1[1:npred])& ( r.lor.1[1:npred] < bsup.garch)
pdn2 = sum(dans2)/length(dans2) # ARMA-GARCH

op<-par(oma=rep(0.5,4))
mat.p = cbind(bsup.garch,
 binf.garch,binf.ar,bsup.ar,r.dan.1[1:npred])
matplot(1:npred,mat.p,type="l", col="black", lty=c(1,1,2,2,3),
  lwd=2,xlab="horizon 50", ylab="rendement")
leg.txt = c("ARMA-GARCH","AR","réalisation")
legend(14,3, leg.txt, lty=c(1,2,3))
par(op)




