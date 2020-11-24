##########################################
# OA experiment - Mudflats 
# script by : C. Hubas 
# Data : C. Mazière, M. Bodo et al 
##########################################

##########################################
# Packages
##########################################

library(ade4)
library(ggplot2)
library(cowplot)
library(vegan)

##########################################
# Aesthetics
##########################################

# palettes
my.palette <- colorRampPalette(c("#DE4064","#FFAA5A","#FEEB94","#06D6A0","#4B7696"))

##########################################
# BC-MFA Analysis
##########################################

#### Upload data for MFA
all.data <- read.csv("Physico-chem_MESO_modif.csv",sep=";",h=T)
no.na.data <- na.omit(all.data)
names(no.na.data)

# split dataset according to sampling time
no.na.data.split<-split(no.na.data,factor(no.na.data$time))

#######################################
# Create a custom function for BC-MFA analysis
make.BC.MFA<-function(data){
  
  # separate matrix
  enviro <- data[,5:14]
  pigments <- data[,15:51]
  eps <- data[,52:55]
  blocs <- c(dim(enviro)[2],
             dim(pigments)[2],
             dim(eps)[2]
  )
  
  # make dataset
  mfa.data <- data.frame(
    enviro,
    pigments,
    eps)
  
  # perform Mltiple factor Analysis (MFA)
  res.mfa <- dudi.pca(mfa.data,
                      col.w = rep(c(1/dudi.pca(enviro,scannf=F,nf=2)$eig[1],
                                    rep(1/dudi.pca(pigments,scannf=F,nf=2)$eig[1]),
                                    rep(1/dudi.pca(eps,scannf=F,nf=2)$eig[1])),blocs),
                      scannf=F,
                      nf=2)
  
  # create control variable (must be factor)
  fac <- factor(data$treatment)
  
  # perform between class analysis on MFA results
  res.bca <- bca(res.mfa, fac, scannf = FALSE, nf = 2)
  plot(res.bca$ls,
       pch=21,
       col="white",
       bg=my.palette(length(levels(fac)))[fac],
       cex=3)
  ordicluster(res.bca$ls,hclust(vegdist(data[,5:55],"euclidean"),"average",))
  randtest(res.bca)
}

# End onf custom BC-MFA function
#######################################

# Apply custom BC-MFA function on each sampling date 
layout(matrix(c(1,3,5,7,2,4,6,8,9,9,9,9),ncol=3))
par(mar=c(2,1,1,1))
# plot individuals
lapply(no.na.data.split,make.BC.MFA)
# Add legend
legend("topleft",
       c("C","pH","TC","TC/pH"),
       text.col = my.palette(4),
       bty="n")
# plot variables for a chosen sampling time (here "P8")
X<-dudi.pca(no.na.data.split$P8[,5:55],scannf=F,nf=2)
Fa<-factor(no.na.data.split$P8$treatment)
coord<-bca(X,Fa,scannf=F,nf=2)$co
plot(coord,type="n")
arrows(x0=0,
       y0=0,
       x1=coord$Comp1,
       y1=coord$Comp2,
       col='lightgrey')
bloc.factor<-factor(rep(c("enviro","pigment","EPS"),c(10,37,4)))
text(x=coord$Comp1,
     y=coord$Comp2,
     rownames(coord),
     cex=1,
     col=rainbow(3)[bloc.factor])
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
title("variables-P8")



##########################################
# Pigment Analysis (short list)
##########################################

all.pig.data <- no.na.data[,c(1:4,15:51)]
names(all.pig.data)
apply(all.pig.data[,-c(1:4)],1,sum)
short.pig.data <- data.frame(
  Cy=all.pig.data$"aCy"+all.pig.data$"bCy",
  Al=all.pig.data$"Al",
  Caroten=all.pig.data$"BB.Car"+all.pig.data$"BE.Car",
  UCar=all.pig.data$"Car1"+all.pig.data$"Car2",
  Ca=all.pig.data$"Ca"+all.pig.data$"Ca.allo"+all.pig.data$"Ca.epi",
  Ca.iso=all.pig.data$"Ca.iso1"+all.pig.data$"Ca.iso2"+all.pig.data$"Ca.iso3"+all.pig.data$"Ca.iso4"+all.pig.data$"Ca.iso5"+all.pig.data$"Ca.iso7",
  BCa=all.pig.data$"Bca",
  Cb=all.pig.data$"Cb",
  Cc=all.pig.data$"Cc2",
  Cda=all.pig.data$"Cda",
  Ct=all.pig.data$"Ct"+all.pig.data$"Ct.iso",
  Ec=all.pig.data$"Ec",
  F=all.pig.data$F+all.pig.data$"F.iso1",
  L=all.pig.data$L+all.pig.data$"L.iso1"+all.pig.data$"L.iso2",
  My=all.pig.data$"My.iso1"+all.pig.data$"My.iso2"+all.pig.data$"My.iso3"+all.pig.data$"My.iso4",
  O=all.pig.data$"O",
  Pda=all.pig.data$"Pda.1",
  Pha=all.pig.data$"Pha",
  Pya=all.pig.data$"Pya",
  Z=all.pig.data$"Z"
)
short.pig.data$treatment<-factor(all.pig.data$treatment)
short.pig.data$time<-factor(all.pig.data$time)

# make lists
short.pig.list<-split(short.pig.data,short.pig.data$time)

#######################################
# Create a function for short.pig dataset
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

# End of custom function
#######################################

# make plots
layout(matrix(c(1,3,5,7,2,4,6,8,9,9,9,9),ncol=3))
par(mar=c(2,1,1,1))
names(short.pig.data)
lapply(short.pig.list,make.subplot.short.pig)
X<-dudi.pca(short.pig.data[short.pig.data$time=="P8",1:20],scannf=F,nf=2)
Fa<-factor(short.pig.data[short.pig.data$time=="P8",]$treatment)
coord<-bca(X,Fa,scannf=F,nf=2)$co
plot(coord,type="n")
arrows(x0=0,y0=0,x1=coord$Comp1,y1=coord$Comp2,col="lightgrey")
text(x=coord$Comp1,y=coord$Comp2,rownames(coord),cex=1.5)
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
title("variables-P8")


##########################################
# Focus on Chlorophyll a "derivatives"like" pigments
##########################################

# boxplots
box.ca.iso<-ggplot(short.pig.data,
       aes(x = time,
           y= Ca.iso,
           col=treatment))+
  geom_boxplot()+
  ylab("Chlorophyll a derivatives (%)")+
  xlab("Time")+
  theme_bw()

box.ca<-ggplot(short.pig.data,
       aes(x = time,
           y= Ca,
           col=treatment))+
  geom_boxplot()+
  ylab("Chlorophyll a (%)")+
  xlab("Time")+
  theme_bw()

box.pha<-ggplot(short.pig.data,
       aes(x = time,
           y= Pha,
           col=treatment))+
  geom_boxplot()+
  ylab("Pheophytin a (%)")+
  xlab("Time")+
  theme_bw()

plot_grid(box.ca,box.ca.iso,box.pha,labels=c("a","b","c"),nrow=3)

# scatter plots
time.j<-c(1,8,15,23,30,37,44,51)[short.pig.data$time]

smooth_iso<-ggplot(short.pig.data,
       aes(x = time.j,
           y= Ca.iso,
           col=treatment))+
  geom_point()+
  geom_smooth(aes(fill=treatment),method = "loess")+
  ylab("Chlorophyll a derivatives (%)")+
  xlab("Time")+
  theme_bw()

smooth_ca<-ggplot(short.pig.data,
       aes(x = time.j,
           y= Ca,
           col=treatment))+
  geom_point()+
  geom_smooth(aes(fill=treatment),method = "loess")+
  ylab("Chlorophyll a (%)")+
  xlab("Time")+
  theme_bw()

smooth_Pha<-ggplot(short.pig.data,
                  aes(x = time.j,
                      y= Pha,
                      col=treatment))+
  geom_point()+
  geom_smooth(aes(fill=treatment),method = "loess")+
  ylab("Pheophytin a (%)")+
  xlab("Time")+
  theme_bw()

plot_grid(smooth_ca,smooth_iso,smooth_Pha,labels=c("a","b","c"),nrow=3)

# Ratio Chlorophyll a like vs. true Chlorophyll a
ggplot(short.pig.data,
                  aes(x = time,
                      y= Ca.iso/Pha,
                      col=treatment))+
  geom_boxplot()+
  ylab("metalation ratio")+
  xlab("Time")+
  theme_bw()

