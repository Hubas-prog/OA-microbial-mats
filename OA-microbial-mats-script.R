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

#### palette
my.palette <- colorRampPalette(c("#DE4064","#FFAA5A","#FEEB94","#06D6A0","#4B7696"))

my.palette(4)

#### Upoload data

raw.data <- read.csv("/Users/cedric.hubas/Downloads/Table_physicochem_mesocosm - Feuille 1.csv",dec=",",h=T)
no.na.data <- na.omit(raw.data)
names(no.na.data)

##########################################################
#### Unsupervised analysis : multiple factor analysis (MFA)
###########################################################

# extraction des k tables
enviro <- no.na.data[,5:14]
pigments <- no.na.data[,15:51]
eps <- no.na.data[,52:53]
blocs <- c(dim(enviro)[2],
           dim(pigments)[2],
           dim(eps)[2]
)

# standardisation par unité d'inertie
mfa.data <- data.frame(
  enviro/dudi.pca(enviro,scannf=F,nf=2)$eig[1],
  pigments/dudi.pca(pigments,scannf=F,nf=2)$eig[1],
  eps/dudi.pca(eps,scannf=F,nf=2)$eig[1]
)

# MFA
res.mfa <- dudi.pca(mfa.data,scannf=F,nf=2)

var <- fviz_pca_var(res.mfa,
                    habillage=rep(c("enviro","pigments","eps"),blocs),
                    repel = TRUE     
)

ind1 <- fviz_pca_ind(res.mfa,
                     habillage = no.na.data$week,
                     geom = "point",
                     legend.title = "Time",
                     addEllipse=T,
                     repel=TRUE,
                     palette=c("#DE4064","#FFAA5A","#FEEB94","#06D6A0","#4B7696")
)


plot_grid(ind1,var,nrow=1,labels=c("a","b"))

##########################################################
#### Supervised analysis : analyse discriminate PPLS-DA
##########################################################

# transformation du facteur en variable indicatrice
var.ind <- dummy(paste(no.na.data$week))
fac <- as.factor(paste(no.na.data$week))

# transformation racine carrée
DATA.M <- as.matrix(sqrt(mfa.data))

# modèle Canonical Powered Partial least squared (PLS) regression
PLSDA<-cppls(var.ind ~ DATA.M, ncomp=20) 
summary (PLSDA)
attributes(PLSDA)

par(mfrow=c(1,2))
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

between.weeks <- bca(mfa.data,as.factor(paste(no.na.data$week)),scannf=F,nf=2)
between.weeks$ratio

between.treaments <- bca(mfa.data,as.factor(paste(no.na.data$treatments)),scannf=F,nf=2)
between.treaments$ratio

#Test pull reverse
# again

