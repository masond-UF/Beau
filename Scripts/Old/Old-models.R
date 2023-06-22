## --------------- HEADER ------------------------------------------------------
## Script name: 3_Models.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-12
## Date Last modified: 2023-06-12
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for assessing the likelihood of detecting
## plants and the likelihood that detected plants are invasive at the feeders
## and paired natural sites.

# probability of detecting plant far lower [fewer plants with feeder]
# plants you do detect are invasive showing shift in relative dominance
# probability of bare ground and herbivory with distance to transect

## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse)

# Clear the decks
rm(list=ls())

# Bring in the data
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv')

## --------------- PREPARE THE DATA --------------------------------------------

# Subset the data
feeder.veg.close <- feeder.veg |>
	filter(Meter < 6)

# Probability of detecting plants
plant.prob <- feeder.veg.close |>
	group_by(Plot, Time, Feeder, Transect, Veg.Type) |>
	dplyr::summarise(Bare.ground = n()) |>
	filter(Veg.Type == 'bare ground') |>
	mutate(Total.points = 5,
				 Plants = Total.points-Bare.ground) |>
	select(Plot, Time, Feeder, Transect, Total.points, Bare.ground, Plants)

## --------------- PROBABILITY OF DETECTING BARE GROUND ------------------------

plant.prob.mod <- glmer(cbind(plant.prob$Bare.ground, plant.prob$Plants) ~ 
						Feeder*Time + (1|Plot),
						data = plant.prob, family = binomial)
Anova(plant.prob.mod)
confint(emmeans(plant.prob.mod, pairwise~Feeder*Time, type = 'response'))
test(emtrends(plant.prob.mod, pairwise~Feeder, var = 'Time'))

## --------------- PROBABILITY OF DETECTING INVASIVE SPECIES -------------------

total.plants <- feeder.veg.close |>
	filter(Invasive != 'na') |>
	group_by(Plot, Time, Feeder, Transect) |>
	summarize(Plants = n())

native.plants <- feeder.veg.close |>
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
emmeans(invasive.plant.prob.mod, ~Feeder*Time, type = 'response')

