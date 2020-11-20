####### OA experiment - Mudflats - La Rochelle #######
####### script by : C. Hubas                   #######
####### Data : C. Mazière, M. Bodo et al       #######


#### Packages

library(ade4)
library(ggplot2)
library(factoextra)
library(cowplot)
library(pls)
library(RVAideMemoire)
library(vegan)

#### palettes
my.palette <- colorRampPalette(c("#DE4064","#FFAA5A","#FEEB94","#06D6A0","#4B7696"))

#### Upload pigment data
all.pig.data <- read.table("Pourcentages.txt",h=T)
all.pig.data$treatment<-factor(gsub(pattern = "PH","pH",all.pig.data$treatment))
short.pig.data <- read.table("POURCENT.txt",h=T)
all.pig.data$treatment

##########################################################
#### Pigment analysis : BCA
###########################################################
all.pig.data2 <- all.pig.data[raw.pig.data$treatment!="Ctr.bis",]
short.pig.data2 <- short.pig.data[raw.pig.data$treatment!="Ctr.bis",]
short.pig.data2$treatment<-all.pig.data2$treatment

all.pig.list<-split(all.pig.data2,all.pig.data2$time)
short.pig.list<-split(short.pig.data2,all.pig.data2$time)

names(short.pig.list$P1)

data=all.pig.data2
f="treatment"
n=37

make.subplot.all.pig<-function(data){
  res.pca <- dudi.pca(data[,1:37],scannf=F,nf=2)
  factor<-factor(data[,which(names(data)=="treatment")])
  res.bca <- bca(res.pca, factor, scannf = FALSE, nf = 2)
  plot(res.bca$ls,
       pch=21,
       col="white",
       bg=my.palette(length(levels(factor)))[factor],
       cex=3)
  ordicluster(res.bca$ls,hclust(vegdist(data[,1:37],"bray"),"ward.D2",))
}

make.subplot.short.pig<-function(data){
  res.pca <- dudi.pca(data[,1:20],scannf=F,nf=2)
  factor<-factor(data[,which(names(data)=="treatment")])
  res.bca <- bca(res.pca, factor, scannf = FALSE, nf = 2)
  plot(res.bca$ls,
       pch=21,
       col="white",
       bg=my.palette(length(levels(factor)))[factor],
       cex=3)
  ordicluster(res.bca$ls,hclust(vegdist(data[,1:20],"bray"),"ward.D2",))
}

par(mfrow=c(2,4),mar=c(0,0,0,0))
lapply(all.pig.list,make.subplot.all.pig)

par(mfrow=c(2,4),mar=c(0,0,0,0))
lapply(all.pig.list,make.subplot.short.pig)


gdat<-pig.data[raw.pig.data$treatment!="Ctr.bis",]
Treatments<-raw.pig.data$treatment[raw.pig.data$treatment!="Ctr.bis"]
Treatments<-gsub(pattern = "PH","pH",Treatments)
ggplot(gdat,aes(x = raw.pig.data$time[raw.pig.data$treatment!="Ctr.bis"],
                y= Ca.iso,
                col=Treatments)) +
    geom_boxplot()+
    theme_bw()


##########################################################
#### Unsupervised analysis : multiple factor analysis (MFA)
###########################################################

#### Upload data for MFA
all.data <- read.csv("Table_physicochem_mesocosm - Feuille 1.csv",dec=",",h=T)
no.na.data <- na.omit(all.data)
names(no.na.data)


# extraction des k tables
enviro <- no.na.data[,5:14]
pigments <- no.na.data[,15:51]
eps <- no.na.data[,52:55]
blocs <- c(dim(enviro)[2],
           dim(pigments)[2],
           dim(eps)[2]
)

# jeu de données
pca.data <- data.frame(
  enviro,
  pigments,
  eps)

# changement des poids
res.mfa <- dudi.pca(pca.data,col.w = rep(c(1/dudi.pca(enviro,scannf=F,nf=2)$eig[1],
                                           rep(1/dudi.pca(pigments,scannf=F,nf=2)$eig[1]),
                                           rep(1/dudi.pca(eps,scannf=F,nf=2)$eig[1])),blocs),
                    scannf=F,nf=2)

plot(res.mfa$li,main="MFA")
s.class(res.mfa$li,
        as.factor(paste(no.na.data$week,no.na.data$treatment)),
        col=my.palette(24),add.plot=T)



##########################################################
#### Supervised analysis : analyse discriminate PPLS-DA
##########################################################

# transformation du facteur en variable indicatrice
var.ind <- dummy(paste(no.na.data$week,no.na.data$treatment,sep=""))
fac <- as.factor(paste(no.na.data$week,no.na.data$treatment,sep=""))

# transformation racine carrée
DATA.M <- as.matrix(decostand(pca.data,"hellinger"))

# modèle Canonical Powered Partial least squared (PLS) regression
PLSDA<-cppls(var.ind ~ DATA.M, ncomp=20) 
summary (PLSDA)
attributes(PLSDA)

par(mfrow=c(2,1))
MVA.plot(PLSDA, "scores", fac=fac, cex = 0.8, col=my.palette(length(levels(fac))))
MVA.plot(PLSDA, "corr", thres = 0, fac=rep(c("enviro","pigments","eps"),blocs), arrows = FALSE, cex = 0.7,intcircle = 0.5,points=F,col=rainbow(3))

# cross model validation (2CV) with different PLS analyses = Classification error rate
valid<-MVA.cmv(DATA.M, fac ,repet=50,ncomp=6,model="PPLS-DA",crit.inn="NMC")

# Test différences entre groupes (si non significatif, absolument test 2 a 2)
testglobal <-MVA.test(DATA.M, fac, model="PPLS-DA",cmv=TRUE,ncomp=6,kout=7,kinn=6)

# permet de voir si groupes sont différents 2 a 2 (étape très longue si grand nombre de variables)#
test2a2 <-pairwise.MVA.test(DATA.M, fac ,model="PPLS-DA",cmv=TRUE,ncomp=6,kout=5,kinn=4)#

##########################################################
#### Supervised analysis 2 : between-class analysis
##########################################################

# BCA interaction
between.interaction <- bca(res.mfa,
                           as.factor(paste(no.na.data$week)),
                           scannf=F,nf=2)
between.interaction$ratio
randtest(between.interaction)
plot(between.interaction)

##########################################################
#### Full analysis : P8
##########################################################

# extraction des k tables
enviro.P8 <- no.na.data[no.na.data$week=="P8",5:14]
pigments.P8 <- no.na.data[no.na.data$week=="P8",15:51]
eps.P8 <- no.na.data[no.na.data$week=="P8",52:55]
blocs.P8 <- c(dim(enviro.P8)[2],
           dim(pigments.P8)[2],
           dim(eps.P8)[2]
)

#  standardisation par unité d'inertie
pca.data.P8 <- data.frame(
  enviro.P8,
  pigments.P8,
  eps.P8)

# MFA
res.mfa.P8 <- dudi.pca(pca.data.P8,col.w = rep(c(1/dudi.pca(enviro.P8,scannf=F,nf=2)$eig[1],
                                                 rep(1/dudi.pca(pigments.P8,scannf=F,nf=2)$eig[1]),
                                                 rep(1/dudi.pca(eps.P8,scannf=F,nf=2)$eig[1])),blocs.P8),
                       scannf=F,nf=2)

plot(res.mfa.P8$li,main="MFA")
s.class(res.mfa.P8$li,
        as.factor(paste(no.na.data$week[no.na.data$week=="P8"],
                        no.na.data$treatmen[no.na.data$week=="P8"])),
        col=my.palette(4),add.plot=T)


# transformation du facteur en variable indicatrice
var.ind.P8 <- dummy(paste(no.na.data$week,no.na.data$treatment,sep="")[no.na.data$week=="P8"])
fac.P8 <- as.factor(paste(no.na.data$week,no.na.data$treatment,sep="")[no.na.data$week=="P8"])

# transformation racine carrée
DATA.M.P8 <- as.matrix(decostand(mfa.data.P8,"hellinger"))

# modèle Canonical Powered Partial least squared (PLS) regression
PLSDA.P8<-cppls(var.ind.P8 ~ DATA.M.P8) 

par(mfrow=c(2,1))
MVA.plot(PLSDA.P8, "scores", fac=fac.P8, cex = 0.8, col=my.palette(length(levels(fac.P8))))
MVA.plot(PLSDA.P8, "corr", thres = 0, fac=rep(c("enviro","pigments","eps"),blocs), arrows = FALSE, cex = 0.7,intcircle = 0.5,points=F,col=rainbow(3))

# BCA treatment
between.treatment.P8 <- bca(res.mfa.P8,fac.P8,scannf=F,nf=2)
between.treatment.P8$ratio

plot(between.treatment.P8$ls,col=my.palette(length(levels(fac.P8)))[fac.P8])
s.class(between.treatment.P8$ls,fac.P8,
        add.plot=T,
        col=my.palette(length(levels(fac.P8))),
        axesell = F,
        cellipse=0
        )
plot(between.treatment.P8$co,type="n")
arrows(x0 = 0,
       y0=0,
       x1=between.treatment.P8$co[,1],
       y1=between.treatment.P8$co[,2],
       col=alpha(rainbow(3)[as.factor(rep(c("enviro","pigments","eps"),blocs))],0.5),
       length=0.2
       )
text(x=between.treatment.P8$co[,1],
     y=between.treatment.P8$co[,2],
     rownames(between.treatment.P8$co),
     col=rainbow(3)[as.factor(rep(c("enviro","pigments","eps"),blocs))],
     cex=0.7
     )
abline(h=0,col="grey")
abline(v=0,col="grey")

randtest(between.treatment.P8)
plot(between.treatment.P8)

