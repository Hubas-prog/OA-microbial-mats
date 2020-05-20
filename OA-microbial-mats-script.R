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

#### Upoload data
raw.data <- read.csv("Table_physicochem_mesocosm - Feuille 1.csv",dec=",",h=T)
no.na.data <- na.omit(raw.data)
names(no.na.data)

##########################################################
#### Unsupervised analysis : multiple factor analysis (MFA)
###########################################################

# extraction des k tables
enviro <- no.na.data[,5:14]
pigments <- no.na.data[,15:51]
eps <- no.na.data[,52:55]
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
                     habillage = paste(no.na.data$week,no.na.data$treatment),
                     geom = "point",
                     legend.title = "Time",
                     addEllipse=T,
                     repel=TRUE,
                     palette=my.palette(24)
)


plot_grid(ind1,var,nrow=2,labels=c("a","b"))

##########################################################
#### Supervised analysis : analyse discriminate PPLS-DA
##########################################################

# transformation du facteur en variable indicatrice
var.ind <- dummy(paste(no.na.data$week,no.na.data$treatment,sep=""))
fac <- as.factor(paste(no.na.data$week,no.na.data$treatment,sep=""))

# transformation racine carrée
DATA.M <- as.matrix(decostand(mfa.data,"hellinger"))

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
between.interaction <- bca(res.mfa,as.factor(paste(no.na.data$treatment,no.na.data$week)),scannf=F,nf=2)
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
mfa.data.P8 <- data.frame(
  enviro.P8/dudi.pca(enviro.P8,scannf=F,nf=2)$eig[1],
  pigments.P8/dudi.pca(pigments.P8,scannf=F,nf=2)$eig[1],
  eps.P8/dudi.pca(eps.P8,scannf=F,nf=2)$eig[1]
)

# MFA
res.mfa.P8 <- dudi.pca(mfa.data.P8,scannf=F,nf=2)

var.P8 <- fviz_pca_var(res.mfa.P8,
                    habillage=rep(c("enviro","pigments","eps"),blocs.P8),
                    repel = TRUE     
)

ind1.P8 <- fviz_pca_ind(res.mfa.P8,
                     habillage = paste(no.na.data$week,no.na.data$treatment)[no.na.data$week=="P8"],
                     geom = "point",
                     legend.title = "Time",
                     addEllipse=T,
                     repel=TRUE,
                     palette=my.palette(4)
)


plot_grid(ind1.P8,var.P8,nrow=2,labels=c("a","b"))

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

