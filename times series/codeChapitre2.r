########################################################
#######   Codes du Chapitre 2
#######   Y. Aragon
########################################################

# Mise à jour effectuée le 17/08/2016

# packages utilisés dans le chapitre
require(caschrono)
require(chron)

# 2.1 Dates et heures
# 2.1.1 Considérations générales 
# 2.1.2 Dates et heures en R
(z = Sys.time())
str(z)
format(z, "%a %b %d %X %Y %Z")

unlist(z)

ll = as.POSIXlt(z, "CET")
str(ll)
unlist(ll)

as.POSIXct(1268736919, origin="1970-01-01", tz="GMT")


# Exemple 2.1
cat(str("2008-05-16"), "", str(as.Date("2008-05-16")), "\n")

x= "1970-01-01"
as.numeric(as.Date(x))
st = as.Date("1999-09-17")
en = as.Date("2000-1-7")
(ll = seq(en, st, by="-1 month"))
str(ll)
xs=strptime(x, "%Y-%m-%d")
str(xs)
unlist(xs)
str(as.Date(x))


# Exemple 2.2 (Lecture de dates dans un chier importé)
don.mois1=read.csv2(file=system.file("/import/Tel_extrait.csv", 
          package="caschrono"),skip=0,stringsAsFactors=FALSE)

str(don.mois1,width=60, strict.width="cut")

# Récupération de la date
date= as.Date(don.mois1$Date.app, "%d/%m/%Y")
str(date)

# Récupération simultanée de la date et de l'heure par chron()
x=chron(dates=don.mois1$Date.app,times=don.mois1$H.deb.app,
           format = c(dates = "d/m/Y", times = "h:m:s"))
x[1:2]


# Récupération simultanée de la date et de l'heure comme heure POSIX
dh = paste(don.mois1$Date.app,don.mois1$H.deb.app)
xp = strptime(dh, "%d/%m/%Y %H:%M:%S")
xp[1:2]


# Extraction des composantes
options(digits=12)
dates(x)[1:4]
times(x)[1:4]
hours(xp); minutes(xp); seconds(xp)
quarters(xp);quarters(x)
months(xp)[1:3]; days(xp)[1:3]; weekdays(xp)[1:3]
years(xp)

# Extraction naïve des éléments de date
heure0 = don.mois1$H.deb.app
heure = substr(heure0,1,2)
minute = substr(heure0,4,5)
seconde = substr(heure0,7,8)
an = substr(don.mois1$Date.app,7,10)
an.num = as.numeric(an)


# 2.2 Les structures de séries temporelles dans R

# 2.2.1 La fonction ts()
# Exemple 2.3 (Série temporelle mensuelle)
set.seed(493)
x = 2 + round(rnorm(15),3)
x.ts = ts(x,start=c(2005,2),frequency=12)
x.ts

str(x.ts)

options(digits = 12)
time(x.ts)
frequency(x.ts)
cycle(x.ts)


# Exemple 2.4 (Série multidimensionnelle)
set.seed(8515)
y = 10 + 0.5 * round(rnorm(30), 3)
xy.ts = ts(cbind(x,y),start=c(2005,2),frequency =12)
str(xy.ts)
xy.ts[1:4,]

# Exemple 2.5 (Série quotidienne)
set.seed(571)
z = round(runif(900), digits = 4) + 5
z.ts = ts(z, start = c(2005, 64), frequency = 365)
time(z.ts)[1:8]
num = 31 + 28 + 4
(frac = num/365)

date.z = seq(as.Date("2005-03-05"), length=900, by="day")
date.z[1:4]

# Exemple 2.6 (Données boursières)
str(EuStockMarkets)

CAC40 = EuStockMarkets[, "CAC"]
str(CAC40)
frequency(CAC40)

options(digits=12)
time(CAC40)[c(1:3,1856:1860)]



# 2.2.2 Récupération de données boursières et séries de classe its
deb = "2006-01-01"; fin = "2009-06-30"
instru = c("^FCHI", "GLE.PA", "BN.PA", "OR.PA")

# Erratum (le package its n'étant plus maintenu
# le code suivant risque de ne plus fonctionner
# csdl= priceIts(instrument=instru,start=deb ,end=fin,
#                  quote="Close")
# colnames(csdl@.Data) = c("Cac40","Socgen","Danone","L_Oréal")

# importer les données en utilisant le package tseries
require("tseries")
deb = "2006-01-01"; fin = "2009-06-30"
instru = c("^FCHI", "GLE.PA", "BN.PA", "OR.PA")
cac40 <- get.hist.quote(instrument = instru[1], start=deb, end=fin, quote="Close")
socgen <- get.hist.quote(instrument = instru[2], start=deb, end=fin, quote="Close")
BN.PA <- get.hist.quote(instrument = instru[3], start=deb, end=fin, quote="Close")
OR.PA <- get.hist.quote(instrument = instru[4], start=deb, end=fin, quote="Close")
csdl <- cbind(cac40 ,socgen, BN.PA, OR.PA)

# Télécharger ensuite le package timeSeries
require("timeSeries")

aa = returns(csdl, percentage = TRUE)
r.csdl = aa[complete.cases(aa) == TRUE,]

cacsoc=union(cac40,socgen)
cs=intersect(cac40,socgen)

deb2 = deb = "2008-02-01"

# la fonction rangeits ne fonctionne plus.
#wxx = rangeIts(cac40, start=deb2, end=fin)

# il faudra à présent utiliser la commande suivante :
wxx <- window(cac40, start=deb2, end=fin)


# 2.3 Retard, différence, union, intersection, sous-série - série de classe ts
aa=scan(system.file("/import/champagne_2001.txt",package="caschrono"))

champa.ts = ts(aa/1000, start=c(2001,1), frequency=12)

chL1 = lag(champa.ts, -1)
chL12 = lag(champa.ts, -12)
dd = ts.union(champa.ts, chL1, chL12)
dd[c(1:2, 12:13), ]

window(dd, start = c(2002, 10), end = c(2003, 2))

dd0 = ts.intersect(champa.ts, chL1, chL12)
dd0[1:3,]


# 2.4 Traitement des manquants
# Repérage des manquants
(xx = c(0:3,99,7,99))
which(xx == 7)   
is.na(xx) = (xx == 99)
xx 
which(is.na(xx) == TRUE)

DF = data.frame(x=c(1,2,3),y=c(0, 10, NA),z=c(NA,1,2))
(DF1 = na.omit(DF))
m=as.matrix(DF)
(m1 = na.omit(m))
(imq = na.action(m1))

# Repérage d'une valeur manquante sur les cours de bourse
manq = !complete.cases(csdl)
i.manq = which(manq == TRUE)
# le package its n'existant plus, il faut utiliser une autre syntaxe
# (date.manq = csdl@dates[i.manq][1:3])
(date.manq = index(csdl)[i.manq][1:3])

# Repérage d'une valeur exceptionnelle.
data(essil)
op=par(mar=c(2.5,2.5,2,2),mgp=c(1.5,.5,0),
 oma=c(0,0,0,0),cex.main=.8,cex.lab=.7,cex.axis=.7)
plot(essil, ylab="cours")
par(op)

fig2.1

# paragraphe "Repérage dune valeur exceptionnelle" 
# remplacer les commandes suivantes :
# i0 = which(essil@.Data > 60)
i0 = which(essil > 60)

# essil@dates[i0]
index(essil)[i0]

# essil@.Data[i0] = 50
essil[i0] <- 50

essil[(i0-2):(i0+2),]

# essil@.Data[i0] = NA
essil[i0] <- NA

# z = zoo(essil@.Data,essil@dates)

# z.corr = na.approx(z)
z.corr <- na.approx(essil)
z.corr[(i0-2):(i0+2),]