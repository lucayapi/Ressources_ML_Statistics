############################################################
######   Codes utilisés dans le chapitre 3
######   Y. Aragon
############################################################
  
  # Mise à jour effectuée le 17/08/2016
  
# packages utilisés dans le chapitre
require(caschrono)
require(fBasics)
require(xtable)

# 3.1 Principe de la régression linéaire
# 3.2 Signifcativité de la régression 
# 3.3 Comparaison de modèles et critères d'information 
# 3.4 Intervalle de confiance (IC)
# 3.5 Prédiction
# 3.6 Exemple : consommation d'électricité 

data(khct)

plot.ts(khct,xlab='temps',main="",cex.lab=.9,cex.axis=.8,
 oma.multi=c(4.5,3,.2,0),mar.multi=c(0, 4, 0, .5),las=0)

fig3-1

ytr= cbind(khct[,"kwh"], log(khct[,"kwh"]), (khct[,"kwh"])^.5)
colnames(ytr)= c("kwh", "log(kwh)","kwh^.5")
my.panel<- function(x, y,..., pf = parent.frame()) {
 lines(x, y, ...)
 abline(lm(y~x),lty=2)
}
plot(zoo(ytr),ylab=c("kwh","log(kwh)",expression(sqrt(kwh))), 
 main="", xlab='temps', panel=my.panel,cex=1.2)


khct.df<-as.data.frame(window(cbind(khct,time(khct)),
 end=c(1983,12)))
colnames(khct.df) <- c("kwh","htdd","cldd","t1")
mod1 = lm(sqrt(kwh) ~ htdd +cldd+ t1, data=khct.df)

# Examen des résidus du modèle 1
khct.df$t1.2 = (khct.df$t1 - 1977)^2
mod2 = lm(sqrt(kwh) ~ htdd + cldd + t1 + t1.2, data=khct.df)

resi.ts=ts(cbind(residuals(mod1),residuals(mod2)),start=c(1970,1),frequency=12)
my.panel2<- function(x, y,..., pf = parent.frame()) {
lines(x, y, ...)
abline(h=0,lty=3,lwd=2)
}

plot(zoo(resi.ts),ylab=c("Ajust. linéaire","Ajust. quadratique"), main="Examen des résidus",
xlab='temps', panel=my.panel2,cex.lab=1.4,cex.main=1.4,cex=1.4)


op=par(oma=rep(0,4))
lag.plot(rev(residuals(mod2)),12,layout=c(4,3),main="Résidus du modèle 2", diag.col="red",cex.main=.7)
par(op)


# Examen de la normalité des résidus du modèle 2.

op=par(mar=c(4,4,3,2),oma=rep(0,4))
plot(mod2, which = 2)
par(op)


a.dag = dagoTest(residuals(mod2))
aa2 = cbind(a.dag@test$statistic ,a.dag@test$p.value)
colnames(aa2)=c( "Stat.","p-value")
rownames(aa2)=c("Omnibus D'Agos.",
       "Skewness D'Agos.","Kurtosis D'Agos.")

xtable(aa2, caption="Tests de normalité des résidus.", digits=4,label="nor.resrcb")

# Examen de l'estimation du modèle 2.
summary(mod2)


# Signicativité de la régression
resum2b = summary(mod2)
vcov2b = vcov(mod2)

# Prévision de la consommation en 1984
khct.df.84<-as.data.frame(window(cbind(khct,time(khct),
(time(khct)-1977)^2), start = c(1984,1)))
colnames(khct.df.84) <- c("kwh","htdd","cldd","t1","t1.2")
p0=predict(mod2, khct.df.84, se.fit=TRUE)
p1=predict( mod2,khct.df.84, interval="confidence",level=0.80,se.fit=TRUE)
p2=predict( mod2,khct.df.84, interval="prediction",level=0.80,se.fit=TRUE)
cbind(p1$fit[1:6,2:3], p2$fit[1:6,2:3])

(etyp.pmco = (p2$fit[,3] - p0$fit)/ qnorm(.9))

aa= seq(as.Date("2000/1/1"), by="month", length.out=12)
id.mois= months(aa, abbreviate = TRUE)

kwh2rc<-window(sqrt(khct[,"kwh"]), start = c(1984,1))
temps2<-time(kwh2rc)
aa= seq(as.Date("2000/1/1"), by="month", length.out=12)
id.mois= months(aa, abbreviate = TRUE)

op=par(oma=rep(0,4),mar=c(4,4,3,2))
plot(kwh2rc, ylim= c(12,15.5), xlab='temps',
     ylab=expression( sqrt(kwh)), xaxt="n")
axis(1, at=seq(1984, 1984.917, length.out=12), labels= id.mois)

# intervalle de confiance
points( p1$fit[,2] ~ temps2, type='l', lty="dotted")
points( p1$fit[,3] ~ temps2,  type='l', lty="dotted")

# intervalle de prédiction
points( p2$fit[,2] ~ temps2, type='l', lty="dashed" )
points( p2$fit[,3] ~ temps2,  type='l', lty="dashed" )
title(main="Année 1984 - Bandes de conf. et de préd. à 80%",
     cex.main=0.8)
legend("topleft", c("Bande de confiance", "Bande de prédiction"), lwd=1, lty=c("dotted", "dashed"),cex=.8 )
par(op)


sum((kwh2rc < p2$fit[,2]) | (kwh2rc > p2$fit[,3])) /length(kwh2rc)