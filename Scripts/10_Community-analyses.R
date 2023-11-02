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

# Manipulative community
manip.comm <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Manipulative") |>
  filter(Meter %in% 1:3) |>
	filter(Feeder == 'Y') |>
	dplyr::select(Plot, Time, Meter, Feeder, Direction, Species) |>
	group_by(Plot, Time, Feeder, Direction, Species) |>
	summarize(Value = n())

# Manipulative gradient
manip.grad <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Manipulative") |>
  filter(Meter %in% 1:3) |>
	filter(!is.na(Invasive)) |>
	filter(Feeder == 'Y')

# Assign functional values to Native
for(i in 1:nrow(manip.grad)){
	if(isTRUE(manip.grad$Invasive[i] == 'Native')){
		manip.grad$Functional[i] <- 'Native'
	}
}

## --------------- CREATE GRADIENT ---------------------------------------------

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
manip.grad$Prop.brows <- round(manip.grad$Prop.browse, 2)

manip.grad <- manip.grad |>
	dplyr::select(Plot, Time, Feeder, Direction, Prop.ruderal, Prop.brows)

## --------------- COMMUNITY ANALYSES ------------------------------------------

# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(devtools)
library(pairwiseAdonis)
library(vegan)
library(permute)

source("Scripts/biostats.R")

# Manipulative
manip.env <- manip.comm[,1:4]
manip.comm <- manip.comm[,5:129]
manip.comm <- manip.comm[, colSums(manip.comm != 0) > 0]

# Drop species
occur <- foa.plots(manip.comm)
rare <- which(occur[,2]<1)
common <- which(occur[,2]>99)

manip.comm <- manip.comm[,-c(rare,common)]

perm <- how(nperm = 199)
setBlocks(perm) <- with(manip.env, Plot)

manip.env$Time <- as_factor(manip.env$Time)
# Blocking by site
mod <- adonis2(manip.comm ~ Time, # Include?
				data = manip.env, permutations = perm)

## --------------- ORDINATION --------------------------------------------------

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

site.scores <- merge(site.scores, manip.grad, keep.all = TRUE,
							by = c('Plot', 'Feeder', 'Time', 'Direction'))

site.scores <- site.scores |>
	filter(!(Plot == 'D2' & Time == 0 & Feeder == 'Y' & Direction == 'N')) |>
		filter(!(Plot == 'D2' & Time == 0 & Feeder == 'Y' & Direction == 'W')) |>
			filter(!(Plot == 'LV6' & Time == 1 & Feeder == 'Y' & Direction == 'E'))

ruderal.1 <- ggplot(data = site.scores |> filter(Time == 0 | Time == 1), aes(x = NMDS1, y = NMDS2, 
	fill = Prop.ruderal, linetype = Time, group = Time)) +
  # scale_color_brewer(palette = "Dark2",
  #                      name = "") +
	scale_fill_viridis_c(option = 'A',
											 direction = -1)+
	geom_point(color = "black", size = 4, shape = 21) +
	scale_x_continuous(limits = c(-5, 2),
										 breaks = c(-5,-4,-3,-2,-1, 0, 1, 2))+
	scale_y_continuous(limits = c(-1.5, 3))+
	stat_ellipse(level = 0.7, linewidth = 1) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

ruderal.2 <- ggplot(data = site.scores |> filter(Time == 0 | Time == 2), aes(x = NMDS1, y = NMDS2, 
	fill = Prop.ruderal, linetype = Time, group = Time)) +
  # scale_color_brewer(palette = "Dark2",
  #                      name = "") +
	scale_fill_viridis_c(option = 'A',
											 direction = -1)+
	geom_point(color = "black", size = 4, shape = 21) +
	scale_x_continuous(limits = c(-5, 2),
										 breaks = c(-5,-4,-3,-2,-1, 0, 1, 2))+
	scale_y_continuous(limits = c(-1.5, 3))+
	stat_ellipse(level = 0.7, linewidth = 1) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

ruderal.3 <- ggplot(data = site.scores |> filter(Time == 0 | Time == 5), aes(x = NMDS1, y = NMDS2, 
	fill = Prop.ruderal, linetype = Time, group = Time)) +
  # scale_color_brewer(palette = "Dark2",
  #                      name = "") +
	scale_fill_viridis_c(option = 'A',
											 direction = -1)+
	geom_point(color = "black", size = 4, shape = 21) +
	scale_x_continuous(limits = c(-5, 2),
										 breaks = c(-5,-4,-3,-2,-1, 0, 1, 2))+
	scale_y_continuous(limits = c(-1.5, 3))+
	stat_ellipse(level = 0.7, linewidth = 1) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

library(patchwork)
ruderal <- ruderal.1+ruderal.2+ruderal.3


browse.1 <- ggplot(data = site.scores |> filter(Time == 0 | Time == 1), aes(x = NMDS1, y = NMDS2, 
	fill = Prop.browse, linetype = Time, group = Time)) +
  # scale_color_brewer(palette = "Dark2",
  #                      name = "") +
	scale_fill_viridis_c(option = 'A',
											 direction = -1)+
	geom_point(color = "black", size = 4, shape = 21) +
	scale_x_continuous(limits = c(-5, 2),
										 breaks = c(-5,-4,-3,-2,-1, 0, 1, 2))+
	scale_y_continuous(limits = c(-1.5, 3))+
	stat_ellipse(level = 0.7, linewidth = 1) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

browse.2 <- ggplot(data = site.scores |> filter(Time == 0 | Time == 2), aes(x = NMDS1, y = NMDS2, 
	fill = Prop.browse, linetype = Time, group = Time)) +
  # scale_color_brewer(palette = "Dark2",
  #                      name = "") +
	scale_fill_viridis_c(option = 'A',
											 direction = -1)+
	geom_point(color = "black", size = 4, shape = 21) +
	scale_x_continuous(limits = c(-5, 2),
										 breaks = c(-5,-4,-3,-2,-1, 0, 1, 2))+
	scale_y_continuous(limits = c(-1.5, 3))+
	stat_ellipse(level = 0.7, linewidth = 1) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

browse.3 <- ggplot(data = site.scores |> filter(Time == 0 | Time == 5), aes(x = NMDS1, y = NMDS2, 
	fill = Prop.browse, linetype = Time, group = Time)) +
  # scale_color_brewer(palette = "Dark2",
  #                      name = "") +
	scale_fill_viridis_c(option = 'A',
											 direction = -1)+
	geom_point(color = "black", size = 4, shape = 21) +
	scale_x_continuous(limits = c(-5, 2),
										 breaks = c(-5,-4,-3,-2,-1, 0, 1, 2))+
	scale_y_continuous(limits = c(-1.5, 3))+
	stat_ellipse(level = 0.7, linewidth = 1) +
	scale_linetype_manual(values=c('solid','dashed')) + 
	theme_bw() +
	theme(legend.position = 'none',
				axis.title = element_blank())

browse <- browse.1+browse.2+browse.3

ruderal/browse
