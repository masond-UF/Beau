## --------------- HEADER ------------------------------------------------------
## Script name: 5_Plant-detection-retrospective.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-22
## Date Last modified: 2023-06-22
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for analyzing the retrospective plant
## feeder data.

## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse) # data munging
library(lme4) # mixed effects glms
library(car) # Anova
library(fitdistrplus) # fitdistrplus
library(styler) # cleaning code
library(emmeans) # generating and comparing means from glms
library(broom) # creating clean dataframes from emmeans
library(broom.mixed) # creating clean dataframes from emmeans
library(performance) # check_model

# Disable exponents
options(scipen = 999)

# Clear the decks
rm(list = ls())

# Bring in the data
retro <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Retrospective") |>
	filter(Meter < 11)

manip <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Manipulative") |>
	filter(Meter < 11)

## --------------- RETROSPECTIVE -----------------------------------------------

retro <- retro |>
	dplyr::select(Plot, Transect, Meter, Feeder, Invasive) |>
	filter(!is.na(Invasive)) |> 
	group_by(Plot, Transect, Feeder, Invasive) |>
	summarize(Count = n())

retro <- retro |>
	pivot_wider(names_from = Invasive, values_from = Count)

retro[is.na(retro)] <- 0
colnames(retro)[5] <- 'Successes'
colnames(retro)[4] <- 'Failures'

retro.mod <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = retro)	
Anova(retro.mod)
summary(retro.mod)
retro.means <- emmeans(retro.mod, ~Feeder, type = 'response')
retro.means <- tidy(retro.means, conf.int = TRUE)

## --------------- MANIPULATIVE ------------------------------------------------

manip <- manip |>
	dplyr::select(Plot, Time, Transect, Meter, Feeder, Invasive) |>
	filter(!is.na(Invasive)) |> 
	filter(Time != 5) |> 
	group_by(Plot, Time, Transect, Feeder, Invasive) |>
	summarize(Count = n())

manip <- manip |>
	pivot_wider(names_from = Invasive, values_from = Count)

manip[is.na(manip)] <- 0
colnames(manip)[6] <- 'Successes'
colnames(manip)[5] <- 'Failures'

manip.mod <- glmer(cbind(Successes, Failures) ~ Feeder * Time + (1|Plot), 
											family = 'binomial', data = manip)	
Anova(manip.mod)
summary(retro.mod)
manip.means <- emmeans(manip.mod, ~Feeder*Time, type = 'response')
manip.means <- tidy(manip.means, conf.int = TRUE)
manip.trends <- emmeans(manip.mod, ~Feeder, var = 'Time')
manip.trends <- tidy(manip.trends, conf.int = TRUE)
emmip(manip.mod, Feeder ~ Time, cov.reduce = range)

## --------------- RETURN SAMPLES ----------------------------------------------

manip <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Manipulative") |>
	filter(Meter < 6)

manip <- manip |>
	dplyr::select(Plot, Time, Direction, Meter, Feeder, Invasive) |>
	filter(!is.na(Invasive)) |> 
	filter(Time == 0 | Time == 5) |> 
	group_by(Plot, Time, Direction, Feeder, Invasive) |>
	summarize(Count = n())

manip <- manip |>
	pivot_wider(names_from = Invasive, values_from = Count)

manip[is.na(manip)] <- 0
colnames(manip)[5] <- 'Successes'
colnames(manip)[6] <- 'Failures'

baseline <- manip |>
	filter(Time == 0)

baseline.mod <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot/Direction), 
											family = 'binomial', data = baseline)	
Anova(baseline.mod)

return <- manip |>
	filter(Time == 0 & Feeder == 'Y' | Time == 5)

# Which sites have two comparisons
plots <- return |>
	group_by(Plot, Time) |>
	summarize(Count = n())

plots <- plots |>
	pivot_wider(names_from = Time, values_from = Count)

plots[is.na(plots)] <- 0
colnames(plots)[2] <- 'Baseline'
colnames(plots)[3] <- 'Return'

plots <- plots |>
	filter(Baseline> 0,
				 Return > 0)

plots <- as.vector(plots$Plot)

return <- return |>
	filter(Plot %in% plots)
	
return.mod <- glmer(cbind(Successes, Failures) ~ as_factor(Time) + (1|Plot/Direction), 
											family = 'binomial', data = return)	
Anova(return.mod)

return.means <- emmeans(return.mod, ~Time, type = 'response')
return.means <- tidy(return.means, conf.int = TRUE)

