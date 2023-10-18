# 5 November 2020—Lab 10
# David Mason
# Set-up workspace ####
library(raster)
library(vegan)

source("Lab 10/Code/biostats.r")

# Import data ####
data(varespec)
data(varechem)

str(varespec)
summary(varespec)

str(varechem)
summary(varechem)

?varespec
?varechem
# Data selection, transformation and standardization ####
# Selecting Species 
occur <- foa.plots(varespec)

rare <- which(occur[,2]<5)
common <- which(occur[,2]>95)

reduced <- varespec[,-c(rare,common)]

# Species transformations and standardizations ####
mapply(hist,as.data.frame(varespec[,1:44]),
			 main=colnames(varespec[,1:44]),xlab="abundance")

log.full <- log1p(varespec)
log.red <- log1p(reduced)

# Full data set:
rsum <- rowSums(log.full)
cv(rsum)

csum <- colSums(log.full)
cv(csum)

# Reduced data set:
rsumRed<-rowSums(log.red)
cv(rsumRed)

csumRed<-colSums(log.red)
cv(csumRed)

# If either the row or column sums have cv >50, standardize by the total:
cSpec <- sweep(log.full, 2, csum, "/")
cSpecRed <- sweep(log.red, 2, csumRed, "/")

# Determine Response Model (RDA vs. CCA) ####
decorana(cSpec)
decorana(cSpecRed)

Vars <- varechem[,c(1,2,7)]
env <- as.data.frame(scale(Vars))

sp.CCA <- cca(cSpec~.,data=env)

f2 <- function(x){
  plot(x~sp.CCA$CC$wa[,1],xlab="CCA AXIS 1", ylab= "Abundance ")
}

# Apply the function across all the species:

mapply(f2,varespec)
# Explanatory Variables ####
Vars <- varechem[,c(1,2,7)]
Vars
round(as.dist(cor(Vars)),2)

cv(colSums(Vars))

env <- as.data.frame(scale(Vars))

# Unconstrained Ordination (CA) ####
# Full Data
ca <- cca(cSpec)
summary(ca)
plot(ca)

#Reduced Data      
ca <- cca(cSpecRed)
summary(ca)
plot(ca)
# Constrained Ordination using CCA ####
sp.CCA <- cca(cSpec~.,data=env)

summary(sp.CCA)
# Monte Carlo testing of the significance of the constrained axis ####
anova(sp.CCA)
anova(sp.CCA,by='axis')
anova(sp.CCA,by='terms')
# Observed (F matrix) and Predicted (Z Matrix) Site Scores ###
summary(sp.CCA)

par(mfrow=c(1,2))

plot(sp.CCA$CC$wa[,1],sp.CCA$CC$wa[,2],xlab="CCA AXIS 1", ylab= "CCA AXIS 2")
plot(sp.CCA$CC$u[,1],sp.CCA$CC$u[,2],xlab="CCA AXIS 1", ylab= "CCA AXIS 2")

spenvcor(sp.CCA)
# Intra-set correlations and biplot scores for the constraining variables ####
sp.CCA$CCA$biplot
# The Tri-Plot (using the site scores from the F matrix): ####
plot(sp.CCA,choices=c(1,2),display=c('wa','sp','bp'),scaling=2)



###### Reduced ######
# Constrained Ordination using CCA ####
sp.CCA <- cca(cSpecRed~.,data=env)

summary(sp.CCA)
# Monte Carlo testing of the significance of the constrained axis ####
anova(sp.CCA)
anova(sp.CCA,by='axis')
anova(sp.CCA,by='terms')
# Observed (F matrix) and Predicted (Z Matrix) Site Scores ###
summary(sp.CCA)

par(mfrow=c(1,2))

plot(sp.CCA$CC$wa[,1],sp.CCA$CC$wa[,2],xlab="CCA AXIS 1", ylab= "CCA AXIS 2")
plot(sp.CCA$CC$u[,1],sp.CCA$CC$u[,2],xlab="CCA AXIS 1", ylab= "CCA AXIS 2")

spenvcor(sp.CCA)
# Intra-set correlations and biplot scores for the constraining variables ####
sp.CCA$CCA$biplot
# The Tri-Plot (using the site scores from the F matrix): ####
plot(sp.CCA,choices=c(1,2),display=c('wa','sp','bp'),scaling=2)

