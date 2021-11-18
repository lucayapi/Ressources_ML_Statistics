###########################################
###### Codes du chapitre 6
######  Y. Aragon
###########################################

# Mise à jour effectuée le 17/08/2016


# packages utilisés dans le chapitre
require(caschrono)
require(expsmooth)


# Introduction 
# 6.1 Lissage exponentiel 
# 6.1.1 Lissage exponentiel simple
# Exemple 6.1
require(expsmooth)
ets0 = ets(fmsales,model="ANN")
summary(ets0)

(aaa=predict(ets0, 4))

e_t = fmsales[62] - ets0$fitted[62]
ets0$fitted[62] + ets0$par[1]*e_t

op=par(mar=c(2,2,2,1)+0.1)
plot(fmsales,xlab="temps",ylab="Ventes",
  main=expression(paste("Série ", plain( fmsales))))
par(op)


op=par(mar=c(2,2,2,1)+0.1)
plot(aaa,xlab="temps",ylab="Ventes",
  main=expression(paste("Prédiction de ",plain( fmsales)," à l'horizon 4")))
par(op)



# 6.1.2 Lissage exponentiel double
# Exemple 6.2
ses.2 = ets(fmsales, model="AAN")
summary(ses.2)

predict(ses.2,h=4)


# 6.1.3 Méthode de Holt-Winters et modèle de lissage correspondant 
