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
library(factoextra)
library(cowplot)
library(pls)
library(RVAideMemoire)
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

# Upload pigment data
all.pig.data <- read.table("Pourcentages.txt",h=T) # all pigments
all.pig.data$treatment<-factor(gsub(pattern = "PH","pH",all.pig.data$treatment))
short.pig.data <- read.table("POURCENT.txt",h=T) # pigments grouping
all.pig.data$treatment

# Data Homogeneization
all.pig.data2 <- all.pig.data[all.pig.data$treatment!="Ctr.bis",]
short.pig.data2 <- short.pig.data[all.pig.data$treatment!="Ctr.bis",]
short.pig.data2$time<-all.pig.data2$time
short.pig.data2$treatment<-all.pig.data2$treatment

# make lists
short.pig.list<-split(short.pig.data2,all.pig.data2$time)

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
lapply(all.pig.list,make.subplot.short.pig)
X<-dudi.pca(short.pig.data2[short.pig.data2$time=="P8",1:20],scannf=F,nf=2)
Fa<-factor(short.pig.data2[short.pig.data2$time=="P8",]$treatment)
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
box.ca.iso<-ggplot(short.pig.data2,
       aes(x = time,
           y= Ca.iso,
           col=treatment))+
  geom_boxplot()+
  ylab("Chlorophyll a derivatives (%)")+
  xlab("Time")+
  theme_bw()

box.ca<-ggplot(short.pig.data2,
       aes(x = time,
           y= Ca,
           col=treatment))+
  geom_boxplot()+
  ylab("Chlorophyll a (%)")+
  xlab("Time")+
  theme_bw()

box.pha<-ggplot(short.pig.data2,
       aes(x = time,
           y= Pha,
           col=treatment))+
  geom_boxplot()+
  ylab("Pheophytin a (%)")+
  xlab("Time")+
  theme_bw()

plot_grid(box.ca,box.ca.iso,box.pha,labels=c("a","b","c"),nrow=3)

# scatter plots
time.j<-c(1,8,15,23,30,37,44,51)[short.pig.data2$time]

smooth_iso<-ggplot(short.pig.data2,
       aes(x = time.j,
           y= Ca.iso,
           col=treatment))+
  geom_point()+
  geom_smooth(aes(fill=treatment),method = "loess")+
  ylab("Chlorophyll a derivatives (%)")+
  xlab("Time")+
  theme_bw()

smooth_ca<-ggplot(short.pig.data2,
       aes(x = time.j,
           y= Ca,
           col=treatment))+
  geom_point()+
  geom_smooth(aes(fill=treatment),method = "loess")+
  ylab("Chlorophyll a (%)")+
  xlab("Time")+
  theme_bw()

smooth_Pha<-ggplot(short.pig.data2,
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
ggplot(short.pig.data2,
                  aes(x = time,
                      y= Ca.iso/Pha,
                      col=treatment))+
  geom_boxplot()+
  ylab("metalation ratio")+
  xlab("Time")+
  theme_bw()
