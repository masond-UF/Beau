## --------------- HEADER ------------------------------------------------------
## Script name: 10_Plant-detection-manipulative.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-10-13
## Date Last modified: 2023-10-13
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for analyzing the manipulative feeder
## vegetation data

## --------------- SET-UP WORKSPACE --------------------------------------------
library(tidyverse)

# Disable exponents
options(scipen = 999)

# Clear the decks
rm(list = ls())

# Retro community
retro.comm <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Retrospective") |>
  filter(Meter %in% 1:3) |>
	dplyr::select(Plot, Time, Meter, Feeder, Transect, Species) |>
	group_by(Plot, Time, Feeder, Transect, Species) |>
	summarize(Value = n())

# Retro gradient
retro.grad <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Retrospective") |>
  filter(Meter %in% 1:3) 

# Assign functional values to Native
for(i in 1:nrow(retro.grad)){
	if(isTRUE(retro.grad$Invasive[i] == 'Native')){
		retro.grad$Functional[i] <- 'Native'
	}
}

retro.grad <- retro.grad |>
	dplyr::select(Plot, Time, Meter, Feeder, Transect, Functional) |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Value = n())

# Manipulative community
manip.comm <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Manipulative") |>
  filter(Meter %in% 1:3) |>
	filter(Time == 0 | Time == 5) |>
	filter(Feeder == 'Y') |>
	dplyr::select(Plot, Time, Meter, Feeder, Direction, Species) |>
	group_by(Plot, Time, Feeder, Direction, Species) |>
	summarize(Value = n())

# Manipulative gradient
manip.grad <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Manipulative") |>
  filter(Meter %in% 1:3) |>
	filter(Time == 0 | Time == 5) |>
	filter(Feeder == 'Y')

# Assign functional values to Native
for(i in 1:nrow(manip.grad)){
	if(isTRUE(manip.grad$Invasive[i] == 'Native')){
		manip.grad$Functional[i] <- 'Native'
	}
}

## --------------- PREPARE DATA ------------------------------------------------

# Turn retro.comm into community matrix
retro.comm <- retro.comm |>
	pivot_wider(names_from = Species, values_from = Value)

retro.comm[is.na(retro.comm)] <- 0

# Calculate gradient values for retro.grad
unique(retro.grad$Functional)

# Create proportion
retro.grad <- retro.grad |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

retro.grad <- retro.grad |>
	pivot_wider(names_from = Functional, values_from = Count)

retro.grad[is.na(retro.grad)] <- 0

retro.grad <- retro.grad |>
	dplyr::select(-`NA`) |>
	mutate(Total = Native + Ruderal + `Browsed perennials`,
				 Prop.ruderal = Ruderal/Total,
				 Prop.brows = `Browsed perennials`/Total)
retro.grad$Prop.ruderal <- round(retro.grad$Prop.ruderal, 2)
retro.grad$Prop.brows <- round(retro.grad$Prop.brows, 2)

retro.grad <- retro.grad |>
	dplyr::select(Plot, Time, Feeder, Transect, Prop.ruderal, Prop.brows)

# Turn manip.comm into community matrix
manip.comm <- manip.comm |>
	pivot_wider(names_from = Species, values_from = Value)

manip.comm[is.na(manip.comm)] <- 0

# Calculate gradient values for manip.comm

unique(manip.grad$Functional)

# Create proportion
manip.grad <- manip.grad |>
	group_by(Plot, Time, Feeder, Direction, Functional) |>
	summarize(Count = n())

manip.grad <- manip.grad |>
	pivot_wider(names_from = Functional, values_from = Count)

manip.grad[is.na(manip.grad)] <- 0

manip.grad <- manip.grad |>
	mutate(Total = Native + Ruderal + `Browsed perennials`,
				 Prop.ruderal = Ruderal/Total,
				 Prop.brows = `Browsed perennials`/Total)
manip.grad$Prop.ruderal <- round(manip.grad$Prop.ruderal, 2)
manip.grad$Prop.brows <- round(manip.grad$Prop.brows, 2)

manip.grad <- manip.grad |>
	dplyr::select(Plot, Time, Feeder, Direction, Prop.ruderal, Prop.brows)

## --------------- COMMUNITY ANALYSES ------------------------------------------

# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(devtools)
library(pairwiseAdonis)
library(vegan)
library(permute)

# Retrospective
retro.env <- retro.comm[,1:4]
retro.comm <- retro.comm[,5:62]
retro.comm <- retro.comm[, colSums(retro.comm != 0) > 0]

perm <- how(nperm = 199)
setBlocks(perm) <- with(retro.env, Plot)

# Blocking by site
mod <- adonis2(retro.comm ~ Feeder, # Include?
				data = retro.env, permutations = perm, 
) # R2 0.03251


# NEED TO RUN OLD VERSION OF R

# # Check homogeneity of variance assumption
# site <- pre.env[,2]
# mod <- betadisper(pre.spec.hell.dist, site)
# anova(mod) # fails
# TukeyHSD(mod) 
# 
# # Pairwise
# pair.mod <- pairwise.adonis2(retro.comm ~ Feeder,
# 				data = retro.env, strata = 'Plot')

# Manipulative
manip.env <- manip.comm[,1:4]
manip.comm <- manip.comm[,5:106]
manip.comm <- manip.comm[, colSums(manip.comm != 0) > 0]

perm <- how(nperm = 199)
setBlocks(perm) <- with(manip.env, Plot)

manip.env$Time <- as_factor(manip.env$Time)
# Blocking by site
mod <- adonis2(manip.comm ~ Time, # Include?
				data = manip.env, permutations = perm
) # R2 0.0238

# # NEED TO RUN OLD VERSION OF R
# 
# # Check homogeneity of variance assumption
# site <- pre.env[,2]
# mod <- betadisper(pre.spec.hell.dist, site)
# anova(mod) # fails
# TukeyHSD(mod) 
# 
# # Pairwise
# pair.mod <- pairwise.adonis2(retro.comm ~ Feeder,
# 				data = retro.env, strata = 'Plot')

## --------------- RETROSPECTIVE ORDINATION ------------------------------------

rm(mod, perm)

# Create ordination
ord <- metaMDS(retro.comm, trymax = 100)

# Using the scores function from vegan to extract
# the site scores and convert to a data.frame
site.scores <- as.data.frame(vegan::scores(ord, "sites"))

# create a column of site names, from the rownames of data.scores
site.scores <- cbind(retro.env, site.scores)

# Using the scores function from vegan to extract
# the species scores and convert to a data.frame
species.scores <- as.data.frame(vegan::scores(ord, "species"))

# create a column of species, from the rownames of species.scores
species.scores$Species <- rownames(species.scores)

# Create ellipse
site.scores$Time <- as_factor(site.scores$Time)

site.scores$New.time <- NA

for(i in 1:nrow(site.scores)){
 	if(isTRUE(site.scores$Time[i] == '2' | site.scores$Time[i] == '3')){
 		 site.scores$New.time[i] <- '2-3'
 	}
 	if(isTRUE(site.scores$Time[i] == '4' | site.scores$Time[i] == '5'
 		 | site.scores$Time[i] == '5')){
 		 site.scores$New.time[i] <- '4-7'
 	}
} 
site.scores$New.time <- as_factor(site.scores$New.time)

site.scores <- site.scores |>
	ungroup() |>
	dplyr::select(-Time)

colnames(site.scores)[6] <- 'Time'

# Time 
retro.time <- ggplot(data = site.scores, aes(x = NMDS1, y = NMDS2, 
													fill = Feeder, linetype = Time, color = Feeder)) +
  scale_color_brewer(palette = "Dark2",
                       name = "") +
	scale_fill_brewer(palette = "Dark2",
                       name = "") +
	stat_ellipse(level = 0.7, linewidth = 1.5) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	geom_point(color = "black", alpha = 0.3, size = 3, shape = 21) +
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

# Ruderal
class(site.scores$Plot)
class(retro.grad$Plot)
class(site.scores$Feeder)
class(retro.grad$Feeder)
class(site.scores$Transect)
class(retro.grad$Transect)
class(site.scores$Time)
class(retro.grad$Time)
retro.grad$Time <- as_factor(retro.grad$Time)

site.scores <- merge(site.scores, retro.grad, keep.all = TRUE,
							by = c('Plot', 'Feeder', 'Transect'))

library(viridis)
retro.ruderal <- ggplot(data = site.scores, aes(x = NMDS1, y = NMDS2, 
													fill = Prop.ruderal, linetype = Feeder, color = Feeder)) +
  scale_color_brewer(palette = "Dark2",
                       name = "") +
	scale_fill_viridis(discrete = FALSE, option = 'H') +
	geom_point(color = "black", alpha = 0.3, size = 3, shape = 21) +
	stat_ellipse(level = 0.5, linewidth = 1.5) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

# Resistant perennials
retro.browse <- ggplot(data = site.scores, aes(x = NMDS1, y = NMDS2, 
													fill = Prop.brows, linetype = Feeder, color = Feeder)) +
  scale_color_brewer(palette = "Dark2",
                       name = "") +
	scale_fill_viridis(discrete = FALSE, option = 'H') +
	geom_point(color = "black", alpha = 0.3, size = 3, shape = 21) +
	stat_ellipse(level = 0.5, linewidth = 1.5) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

## --------------- MANIPULATIVE ORDINATION -------------------------------------

# # Bring in the code to clean up community data
# source("Scripts/biostats.r")
# 
# # Selecting Species 
# occur <- foa.plots(manip.comm)
# 
# rare <- which(occur[,2]<5)
# common <- which(occur[,2]>95)
# 
# red.manip.comm <- manip.comm[,-c(rare,common)]
# 
# # Drop zero rows
# red.comb <- cbind(manip.env, red.manip.comm)
# red.comb$Total <- rowSums(red.comb[ , c(5:16)], na.rm=TRUE)
# red.comb <- red.comb |>
# 	filter(Total > 0)
# 
# red.manip.comm <- red.comb[5:16]
# red.manip.env<- red.comb[1:4]

# Create ordination
ord <- metaMDS(manip.comm, trymax = 100)

# Using the scores function from vegan to extract
# the site scores and convert to a data.frame
site.scores <- as.data.frame(vegan::scores(ord, "sites"))

# create a column of site names, from the rownames of data.scores
site.scores <- cbind(manip.env, site.scores)

# Using the scores function from vegan to extract
# the species scores and convert to a data.frame
species.scores <- as.data.frame(vegan::scores(ord, "species"))

# create a column of species, from the rownames of species.scores
species.scores$Species <- rownames(species.scores)

# Create ellipse
site.scores$Time <- as_factor(site.scores$Time)

# Time 
manip.time <- ggplot(data = site.scores, aes(x = NMDS1, y = NMDS2, 
													fill = Time, linetype = Time, color = Time,
													group = Time)) +
  scale_color_brewer(palette = "Dark2",
                       name = "") +
	scale_fill_brewer(palette = "Dark2",
                       name = "") +
	geom_point(color = "black", alpha = 0.3, size = 3, shape = 21) +
	stat_ellipse(level = 0.7, linewidth = 1.5) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

# Ruderal
class(site.scores$Plot)
class(retro.grad$Plot)
class(site.scores$Feeder)
class(retro.grad$Feeder)
class(retro.grad$Transect)
class(site.scores$Time)
class(retro.grad$Time)
manip.grad$Time <- as_factor(manip.grad$Time)

site.scores <- merge(site.scores, manip.grad,
							by = c('Plot', 'Feeder', 'Direction', 'Time'))

library(viridis)
manip.ruderal <- ggplot(data = site.scores, aes(x = NMDS1, y = NMDS2, 
													fill = Prop.ruderal, linetype = Time, color = Time,
													group = Time)) +
  scale_color_brewer(palette = "Dark2",
                       name = "") +
	scale_fill_viridis(discrete = FALSE, option = 'H') +
	geom_point(color = "black", alpha = 0.3, size = 3, shape = 21) +
	stat_ellipse(level = 0.5, linewidth = 1.5) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

# Resistant perennials
manip.browse <- ggplot(data = site.scores, aes(x = NMDS1, y = NMDS2, 
													fill = Prop.brows, linetype = Time, color = Time,
													group = Time)) +
  scale_color_brewer(palette = "Dark2",
                       name = "") +
	scale_fill_viridis(discrete = FALSE, option = 'H') +
	geom_point(color = "black", alpha = 0.3, size = 3, shape = 21) +
	stat_ellipse(level = 0.5, linewidth = 1.5) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())
## --------------- COMBINE ORDINATIONS -----------------------------------------

library(patchwork)
(manip.time | manip.ruderal | manip.browse) / (retro.time | retro.ruderal | retro.browse)


