####################################################
#######  Codes du chapitre 7
#######  Y. Aragon
####################################################

# Mise à jour effectuée le 17/08/2016


# packages utilisés dans le chapitre
require(caschrono)
require(dse)
require(polynom)
require(TSA)


# Introduction 
# 7.1 Principe
# 7.2 Simulation de séries temporelles
# 7.2.1 Principe 

# 7.2.2 Illustration numérique
# Exemple 7.1
# A
set.seed(951)
n2 = 250
z0 = rnorm(n2+2,sd=sqrt(1.5))
require("Hmisc")
vt = z0- 0.3*Lag(z0,1) + 0.6*Lag(z0,2)
str(vt)
vt = vt[-(1:2)]

moy=-.5
cc= moy*(1+0.8)
wt = vt+cc
y.n=filter(wt[-1],c(-0.8),method="recursive")
y.n=y.n[-(1:50)]

# B
yc.n=moy+arima.sim(n=200,list(ar=-0.8,ma=c(-0.3,0.6)),
 innov=z0[53:252],n.start=50, start.innov= z0[1:50])

set.seed(281)
yd.n=moy+arima.sim(n=200,list(ar=-0.8,ma=c(-0.3, 0.6)),
          sd = sqrt(1.5),n.start=50)


# C
require(dse)
vecauto= c(1,0.8)
(AR =array( vecauto,c(length(vecauto),1,1)))
vecma = c(1, - 0.3,0.6)
MA=array(vecma,c(length(vecma),1,1))
MA
Cmat = 1
mod1 = ARMA(A=AR  ,B=MA, C=Cmat)


# Exemple 7.2
y.int = diffinv(y.n)

y2.int=moy*(0:199)+arima.sim(n=199,list(order=c(1,1,2),
       ar=-0.8,ma=c(-0.3, 0.6)),sd=sqrt(1.5),n.start=50)


# Exemple 7.3
autopol = polynomial(c(1,0.8))*polynomial(c(1,0,0,0,-0.7))
# vérification (facultative) de la stationnarité
Mod(polyroot(autopol))
# simulation
ys=4+arima.sim(n=200,list(ar=-autopol[-1],ma=c(0, 0.6)),
                sd=sqrt(1.5))

bt=arima.sim(n=230, list (ar=c(0,0,0,+.7)), sd=sqrt(1.5))
yt=4+arima.sim(n=200, list (ar=-0.8, ma = c(0, 0.6)),innov=bt[31:230], n.start=30, start.innov=bt[1:30])


# 7.3 Construction de séries autorégressives
# Exemple 7.4
(autopol=polynomial(c(1,-1/1.4))*polynomial(c(1,1/1.2))*
         polynomial(c(1,-1/2)))



# 7.4 Construction de séries subissant une intervention
# 7.4.1 Réponses typiques à une intervention
# Exemple 7.5
imp.fun= function(temps,t0)
{a=rep(0,length(temps)); a[t0]=1 ;a }
t0 = 3
dates=1:20; ldat=length(dates)
y.imp = imp.fun(dates,t0)
delta = .8; omega = 4 ; vecauto= c(1,-delta)
AR =array( vecauto,c(length(vecauto),1,1))
vecma = 1; MA=array(vecma,c(length(vecma),1,1))
mod1 = ARMA(A=AR  ,B=MA)
y1 = as.vector(simulate(mod1,y0=0,noise=as.matrix(omega*y.imp),
        sampleT=ldat)$output)
y1[1:6]

ech.fun= function(temps,t0)
{a=rep(0,length(temps)); a[temps >= t0]=1 ;a }
y.ech = ech.fun(dates,t0)
y2 = as.vector(simulate(mod1,y0=0,noise=as.matrix(omega*y.ech),
          sampleT=ldat)$output)
y2[1:6]


# 7.4.2 Simulation d'une intervention
temps = 1:100 ; t0 = 51
echel =   rep(1,length(temps))*(temps >= t0)
delta = .6
omega = -2
vecauto= c(1,-delta)
(AR =array( vecauto,c(length(vecauto),1,1)))
vecma = -2
(MA=array(vecma,c(length(vecma),1,1)))
mod2 = ARMA(A=AR  ,B=MA)
moy0 = as.matrix(rep(3.1,100))
moy = moy0+simulate(mod2, y0=0, 
  noise=as.matrix(echel),sampleT=100)$output

set.seed(329)
vecauto1 = c(1,-.8)
mod3 = ARMA(A=array( vecauto,c(length(vecauto),1,1)), B=1 )
bruit1 = simulate(mod3, y0=0,  sd=.5,sampleT=100)
bruit2 = simulate(mod3, y0=0,  sd=1,sampleT=100)
bruit3 = simulate(mod3, y0=0,  sd=1.5,sampleT=100)

# datation de la série
ts.temps= ts( temps, start=c(60,3), frequency=4)
str(ts.temps)
signal1 = moy+bruit1$output
signal2 = moy+bruit2$output
signal3 = moy+bruit3$output
ser2 = cbind(moy,signal1,signal2,signal3)
colnames(ser2)= c("Moyenne","sigma=0.5","sigma=1","sigma=1.5")
sign.br = ts(ser2, start=c(60,3), frequency=4)

plot(sign.br, xlab="temps",  main="Moyenne et simulation - écart-type = .5, 1, 1.5")


# 7.4.3 Estimation d'une intervention
require(TSA)
temps= 1:length(signal1)
echel1 = rep(1,length(temps))* (temps >= 51)
(pirate.m1=arimax(signal1,order=c(1,0,0),xtransf=data.frame(echel1),
        transfer=list(c(1,0)),method="ML"))


