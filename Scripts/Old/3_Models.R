## --------------- HEADER ------------------------------------------------------
## Script name: 3_Models.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-14
## Date Last modified: 2023-06-15
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for analyzing the clean feeder
## vegetation data

## What about just doing the effect size at each one

## What about the 

## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse)
library(DataExplorer)
library(fitdistrplus)
library(performance)

# Clear the decks
rm(list=ls())

# Bring in the data
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative')

invasives <- feeder.veg |> 
	filter(Invasive == 'Invasive')
unique(invasives$Species)

## --------------- CLOSE VEGETATION JUST INVASIVE ------------------------------

# Subset the data
close <- feeder.veg |>
	filter(Meter < 10) # 5

# Model for functional invasive vs plants
rm(list=ls()[! ls() %in% c("feeder.veg", "close")])

for(i in 1:nrow(close)){
	if(is.na(close$Functional[i])){
		close$Functional[i] <- 'Native'
	}
}

for(i in 1:nrow(close)){
	if(isTRUE(close$Functional[i] == 'Ruderal')){
		close$Functional[i] <- 'Browse resistant'
	}
	if(isTRUE(close$Functional[i] == 'Resistant perennials')){
		close$Functional[i] <- 'Browse resistant'
	}
}

unique(close$Functional)

# Tally invasive functional groups
func.invasive <- close |>
	filter(Functional != 'Native') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive |>
	mutate(Total = `Browse resistant`+
				 	`Browsed perennials`) |>
	mutate(`Browse resistant` = `Browse resistant`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total)

# Composition failed [try dirichlet with just invasives?]
func.invasive$Composition <- DR_data(func.invasive[4:5])
func.invasive$Time <- as_factor(func.invasive$Time)
summary(DirichReg(formula = Composition ~ Feeder * Time, model = 'common', 
									data = func.invasive, weights = Total))

fig <- func.invasive |>
	dplyr::select(-Composition, -Total) |>
	pivot_longer(cols = 4:5, names_to = 'Functional', values_to = 'Proportion')

close.fig <- ggplot(fig, aes(x = Time, y = Proportion, fill = Feeder))+
	geom_boxplot()+
	facet_wrap(~Functional)+
	ggtitle('0-10 m from feeder')

## --------------- FAR VEGETATION ----------------------------------------------

# Subset the data
far <- feeder.veg |>
	filter(Meter %in% (11:25)) # 5

# Model for functional invasive vs plants
rm(list=ls()[! ls() %in% c("feeder.veg", "far", 'close.fig')])

for(i in 1:nrow(far)){
	if(is.na(far$Functional[i])){
		far$Functional[i] <- 'Native'
	}
}

for(i in 1:nrow(far)){
	if(isTRUE(far$Functional[i] == 'Ruderal')){
		far$Functional[i] <- 'Browse resistant'
	}
	if(isTRUE(far$Functional[i] == 'Resistant perennials')){
		far$Functional[i] <- 'Browse resistant'
	}
}

unique(far$Functional)

# Tally invasive functional groups
func.invasive <- far |>
	filter(Functional != 'Native') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive |>
	mutate(Total = `Browse resistant`+
				 	`Browsed perennials`) |>
	mutate(`Browse resistant` = `Browse resistant`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total)

# Composition failed [try dirichlet with just invasives?]
func.invasive$Composition <- DR_data(func.invasive[4:5])
func.invasive$Time <- as_factor(func.invasive$Time)
summary(DirichReg(formula = Composition ~ Feeder * Time, model = 'common', 
									data = func.invasive, weights = Total))

fig <- func.invasive |>
	dplyr::select(-Composition, -Total) |>
	pivot_longer(cols = 4:5, names_to = 'Functional', values_to = 'Proportion')

far.fig <- ggplot(fig, aes(x = Time, y = Proportion, fill = Feeder))+
	geom_boxplot()+
	facet_wrap(~Functional)+
	ggtitle('11-25 m from feeder')

library(patchwork)
close.fig/far.fig

## --------------- CLOSE VEGETATION NATIVE VS INTRO ----------------------------

# Model for total invasive vs native
total.plants <- close |>
	filter(!is.na(Invasive)) |> # Bare ground and litter have no invasive status
	group_by(Plot, Time, Feeder, Transect) |>
	summarize(Plants = n())

native.plants <- close |>
	filter(Invasive == 'Native') |>
	group_by(Plot, Time, Feeder, Transect) |>
	summarize(Native = n())

invasive.plant.prob <- merge(total.plants, native.plants, all = TRUE)

invasive.plant.prob <- invasive.plant.prob |>
	mutate(Introduced = Plants-Native)

invasive.plant.prob.mod <- glmer(cbind(invasive.plant.prob$Introduced, 
															invasive.plant.prob$Native) ~ 
						Feeder*Time + (1|Plot/Transect),
						data = invasive.plant.prob, family = binomial)
Anova(invasive.plant.prob.mod)


## --------------- CLOSE VEGETATION MULTI WITH NATIVE --------------------------

# Tally invasive functional groups
func.invasive <- close |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

# Add zeroes
func.invasive <- func.invasive |>
pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive |>
	mutate(Total = `Browse resistant`+
				 	`Browsed perennials`+Native) |>
	mutate(`Browse resistant` = `Browse resistant`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total,
				 Native = Native/Total)

# Composition failed [try dirichlet with just invasives?]
func.invasive$Composition <- DR_data(func.invasive[5:7])
func.invasive$Time <- as_factor(func.invasive$Time)
func.invasive.baseline <- func.invasive |> filter(Time == '0')
summary(DirichReg(formula = Composition ~ Feeder, model = 'common', 
									data = func.invasive.baseline, weights = Total))

func.invasive.post <- func.invasive |> filter(Time %in% c('1', '2', '3'))
summary(DirichReg(formula = Composition ~ Feeder * Time, model = 'common', 
									data = func.invasive, weights = Total))

plot(func.invasive$Composition, cex = 0.1, reset_par = FALSE)
points(toSimplex(func.invasive$Composition), pch = 16, cex = 0.5, col = c('red', 'blue'))

fig <- func.invasive |>
	dplyr::select(-Composition, -Total) |>
	pivot_longer(cols = 5:7, names_to = 'Functional', values_to = 'Proportion')

ggplot(fig, aes(x = Functional, y = log(Proportion+1)))+
	geom_boxplot()+
	facet_wrap(~Feeder*Time)

## Just invasives

## --------------- Messing with a ruderal only model ---------------------------

# Model the likelihood of detecting functional role
# Summarize 
# Model the 'abundance'

# binary at site or not
# multipel trials for probability of detection at transect points where present

# Subset the data
close <- feeder.veg |>
	filter(Meter < 6) # 5

ruderal <- close |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

ruderal <- ruderal |>
	pivot_wider(names_from = Functional, values_from = Count)
ruderal[is.na(ruderal)] <- 0

ruderal <- ruderal |>
	pivot_longer(cols = 5:7, names_to = 'Functional', values_to = 'Count')

ruderal <- ruderal |>
	filter(Functional == 'Ruderal')

for (i in 1:nrow(ruderal)){
	if (ruderal$Count[i] > 0) {
		ruderal$Count[i] <- 1
	}
}

class(ruderal$Time)
m1 <- glm(Count ~ Feeder*Time + Plot, data = ruderal, family = 'binomial')
Anova(m1, type = 2)

emmeans(m1, pairwise ~ Feeder, type = 'response')

sites <- ruderal[ruderal$Count == 1,1:2]

# Subset the data
close <- feeder.veg |>
	filter(Meter < 6)

ruderal <- close |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

ruderal <- ruderal |>
	pivot_wider(names_from = Functional, values_from = Count)
ruderal[is.na(ruderal)] <- 0

# Model for total invasive vs native
total.invasive <- close |>
	filter(Invasive == 'Invasive') |> # Bare ground and litter have no invasive status
	group_by(Plot, Time, Feeder, Transect) |>
	summarize(Plants = n())

ruderal.plants <- close |>
	filter(Functional == 'Ruderal') |>
	group_by(Plot, Time, Feeder, Transect) |>
	summarize(Ruderal = n())


## --------------- Permanova approach ------------------------------------------

# Tally invasive functional groups
func.invasive <- close |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive[,-5]

func.invasive <- func.invasive |>
	mutate(Total = `Browse resistant`+
				 	`Browsed perennials`) |>
	mutate(`Browse resistant` = `Browse resistant`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total)
func.invasive[is.na(func.invasive)] <- 0

library(vegan)
env <- func.invasive[,1:4]
perm <- how(nperm = 199)
setBlocks(perm) <- with(env, Plot)

# Blocking by site
m1 <- adonis2(func.invasive[5:6] ~ Plot + Feeder * Time, # Include?
				data = env, permutations = perm, method = 'euclidean'
) 

## --------------- Effect size approach ----------------------------------------

# Subset the data
close <- feeder.veg |>
	filter(Meter < 10) # 5

# Model for functional invasive vs plants
rm(list=ls()[! ls() %in% c("feeder.veg", "close")])

for(i in 1:nrow(close)){
	if(is.na(close$Functional[i])){
		close$Functional[i] <- 'Native'
	}
}

for(i in 1:nrow(close)){
	if(isTRUE(close$Functional[i] == 'Ruderal')){
		close$Functional[i] <- 'Browse resistant'
	}
	if(isTRUE(close$Functional[i] == 'Resistant perennials')){
		close$Functional[i] <- 'Browse resistant'
	}
}

unique(close$Functional)

# Tally invasive functional groups
func.invasive <- close |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive[,-4]

# Calculate the effect size
feeder <- func.invasive |> filter(Feeder == 'Y')
feeder <- feeder[,-3]
colnames(feeder)[3] <- "Feeder.Browse.Resistant"
colnames(feeder)[4] <- "Feeder.Browsed.Perennials"

no.feeder <- func.invasive |> filter(Feeder == 'N')
no.feeder <- no.feeder[,-3]
colnames(no.feeder)[3] <- "Control.Browse.Resistant"
colnames(no.feeder)[4] <- "Control.Browsed.Perennials"
no.feeder <- no.feeder[, 3:4]

effect.df <- cbind(feeder,no.feeder)
effect.df <- effect.df |>
	mutate("Browse Resistant LLR" = log((Feeder.Browse.Resistant+1)
																			/(Control.Browse.Resistant+1)),
				 "Browsed Perennials LLR" = log((Feeder.Browsed.Perennials+1)
																			/(Control.Browsed.Perennials+1)))

effect.df <- effect.df |>
	select(Plot, Time, "Browse Resistant LLR", "Browsed Perennials LLR")

## --------------- 0-5 meters --------------------------------------------------

feeder.veg <- feeder.veg |> filter(Dataset == 'Manipulative')

# Subset the data
close <- feeder.veg |>
	filter(Meter < 6) # 5

# Model for functional invasive vs plants
rm(list=ls()[! ls() %in% c("feeder.veg", "close")])

for(i in 1:nrow(close)){
	if(is.na(close$Functional[i])){
		close$Functional[i] <- 'Native'
	}
}
unique(close$Functional)

# Tally invasive functional groups
func.invasive <- close |>
	filter(Functional != 'Native') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive |>
	mutate(Total = `Browsed perennials`+
				 	`Resistant perennials` + Ruderal) |>
	mutate(`Resistant perennials` = `Resistant perennials`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total,
				 Ruderal = Ruderal/Total)

# Composition failed [try dirichlet with just invasives?]
func.invasive$Composition <- DR_data(func.invasive[4:6])
func.invasive$Time <- as_factor(func.invasive$Time)
summary(DirichReg(formula = Composition ~ Feeder * Time, model = 'common', 
									data = func.invasive, weights = Total))

fig <- func.invasive |>
	dplyr::select(-Composition, -Total) |>
	pivot_longer(cols = 4:6, names_to = 'Functional', values_to = 'Proportion')

close.fig <- ggplot(fig, aes(x = Time, y = Proportion, fill = Feeder))+
	geom_boxplot()+
	facet_wrap(~Functional)+
	ggtitle('0-5 m from feeder')

## --------------- 0-5 meters --------------------------------------------------

# Subset the data
mid <- feeder.veg |>
	filter(Meter %in% 5:10) # 5

# Model for functional invasive vs plants
rm(list=ls()[! ls() %in% c("feeder.veg", "mid")])

for(i in 1:nrow(mid)){
	if(is.na(mid$Functional[i])){
		mid$Functional[i] <- 'Native'
	}
}

unique(mid$Functional)

# Tally invasive functional groups
func.invasive <- mid |>
	filter(Functional != 'Native') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive |>
	mutate(Total = `Browsed perennials`+
				 	`Resistant perennials` + Ruderal) |>
	mutate(`Resistant perennials` = `Resistant perennials`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total,
				 Ruderal = Ruderal/Total)

# Composition failed [try dirichlet with just invasives?]
func.invasive$Composition <- DR_data(func.invasive[4:6])
func.invasive$Time <- as_factor(func.invasive$Time)
summary(DirichReg(formula = Composition ~ Feeder * Time, model = 'common', 
									data = func.invasive, weights = Total))

fig <- func.invasive |>
	dplyr::select(-Composition, -Total) |>
	pivot_longer(cols = 4:6, names_to = 'Functional', values_to = 'Proportion')

mid.fig <- ggplot(fig, aes(x = Time, y = Proportion, fill = Feeder))+
	geom_boxplot()+
	facet_wrap(~Functional)+
	ggtitle('0-10 m from feeder')
