## --------------- HEADER ------------------------------------------------------
## Script name: 4_Plant-detection-manipulative.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-23
## Date Last modified: 2023-06-22
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for analyzing the manipulative feeder
## vegetation data

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

# Clear the decks
rm(list=ls())

# Bring in the data
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative')

# ditch.ls <- c('P2', 'P3', 'P5', 'P6', 'P7', 'P8', 'WS4', 'WS8')

## --------------- PREPARE ALL DATA --------------------------------------------

# Subset the data
feeder.veg <- feeder.veg |>
	filter(Meter %in% 1:25) # 5
	# filter(!Plot %in% ditch.ls)

# Assign a value to native species in the functional category
for(i in 1:nrow(feeder.veg)){
	if(is.na(feeder.veg$Functional[i]) & isTRUE(feeder.veg$Invasive[i] == 'Native')){
		feeder.veg$Functional[i] <- 'Native'
	}
}

# Assign a value to no nplants species in the functional category
for(i in 1:nrow(feeder.veg)){
	if(is.na(feeder.veg$Functional[i]) & is.na(feeder.veg$Invasive[i])){
		feeder.veg$Functional[i] <- 'Not plant'
	}
}

## --------------- PROPORTION CLOSE INVASIVE/NATIVE ----------------------------

# Subset the data
feeder.veg.sum <- feeder.veg |>
	filter(Meter < 11) 

# Create a df with all site/plot/time/feeder combinations
sites <- feeder.veg.sum |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Everything = n())

# Tally total points that are plants
total <- feeder.veg.sum |>
	filter(Meter < 11) |>
	filter(Species != 'Litter') |>
	filter(Species != 'litter') |>
	filter(Species != 'Bare Ground') |>
	filter(Species != 'Coarse Woody Debris') |>
	filter(Species != 'water') |>
	filter(!is.na(Invasive)) |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

# Introduce zeroes
total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Tally total counts for invasive/native
Invasive <- feeder.veg.sum |>
	filter(Meter < 11) |>
	filter(Species != 'Litter') |>
	filter(Species != 'litter') |>
	filter(Species != 'Bare Ground') |>
	filter(Species != 'Coarse Woody Debris') |>
	filter(Species != 'water') |>
	filter(!is.na(Invasive)) |>
	group_by(Plot, Time, Feeder, Transect, Invasive) |>
	summarize(Count = n())

Invasive <- Invasive |>
	pivot_wider(names_from = 'Invasive', values_from = 'Count')
Invasive[is.na(Invasive)] <- 0
Invasive <- Invasive |>
	pivot_longer(cols = 5:6, names_to = 'Invasive', values_to = 'Successes')

feeder.veg.sum <- merge(Invasive, total, all = TRUE)
feeder.veg.sum <- feeder.veg.sum |>
	filter(!is.na(Invasive)) # not sure why necessary, keep NA

feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Invasive, Successes,
								Failures, Total)

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg.sum", "feeder.veg")])

# Invasive
invasive <- feeder.veg.sum |>
	filter(Invasive == 'Invasive')
invasive.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = invasive, family = 'binomial')
Anova(invasive.mod)
check_model(invasive.mod)
invasive.means <- emmeans(invasive.mod, ~Feeder*Time, type = 'response')
invasive.means <- tidy(invasive.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

# Native
native <- feeder.veg.sum |>
	filter(Invasive == 'Native')
native.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = native, family = 'binomial')
Anova(native.mod)
check_model(native.mod)
native.means <- emmeans(native.mod, ~Feeder*Time, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Invasive <- 'Native'

close.means <- rbind(invasive.means, native.means)
close.means <- close.means |>
	mutate_if(is.numeric, round, digits = 3)
close.means$Transect <- "Close"
close.means <- close.means |>
	dplyr::select(Transect, Invasive, Time, Feeder, prob, everything())

# Clear the decks
rm(list = ls()[!ls() %in% c("close.means", "feeder.veg")])

pos <- position_dodge(.6)
ggplot(close.means, aes(x = as_factor(Time), y = prob, group = interaction(Feeder,Time),
												fill = Feeder))+
	geom_pointrange(aes(ymin = conf.low, ymax = conf.high), 
                    position = pos, size = 1)+
	geom_point(position = pos, size = 4, shape = 21)+
	scale_fill_manual(values = c('gray', 'skyblue'),
										guide = guide_legend(reverse = TRUE),
										labels = c('Control', 'Feeder'))+
	ylab("Probability")+
	xlab('Year')+
	facet_wrap(~Invasive)+
	ggtitle('Probability plant is invasive and native <10 m')+
	theme_bw()

## --------------- PROPORTION FAR INVASIVE/NATIVE ------------------------------

rm(list=ls()[! ls() %in% c("feeder.veg")])

# Subset the data
feeder.veg.sum <- feeder.veg |>
	filter(Meter %in% 11:25) 

# Create a df with all site/plot/time/feeder combinations
sites <- feeder.veg.sum |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Everything = n())

# Tally total points that are plants
total <- feeder.veg.sum |>
	filter(Meter %in% 11:25) |>
	filter(Species != 'Litter') |>
	filter(Species != 'litter') |>
	filter(Species != 'Bare Ground') |>
	filter(Species != 'Coarse Woody Debris') |>
	filter(Species != 'water') |>
	filter(!is.na(Invasive)) |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

# Introduce zeroes
total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Tally total counts for invasive/native
invasive <- feeder.veg.sum |>
	filter(Meter %in% 11:25) |>
	filter(Species != 'Litter') |>
	filter(Species != 'litter') |>
	filter(Species != 'Bare Ground') |>
	filter(Species != 'Coarse Woody Debris') |>
	filter(Species != 'water') |>
	filter(!is.na(Invasive)) |>
	group_by(Plot, Time, Feeder, Transect, Invasive) |>
	summarize(Count = n())

invasive <- invasive |>
	pivot_wider(names_from = 'Invasive', values_from = 'Count')
invasive[is.na(invasive)] <- 0
invasive <- invasive |>
	pivot_longer(cols = 5:6, names_to = 'Invasive', values_to = 'Successes')

feeder.veg.sum <- merge(invasive, total, all = TRUE)
feeder.veg.sum <- feeder.veg.sum |>
	filter(!is.na(Invasive)) # not sure why necessary, keep NA

feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Invasive, Successes,
								Failures, Total)

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg.sum", "feeder.veg")])

# Invasive
invasive <- feeder.veg.sum |>
	filter(Invasive == 'Invasive')
invasive.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = invasive, family = 'binomial')
Anova(invasive.mod)
check_model(invasive.mod)
invasive.means <- emmeans(invasive.mod, ~Feeder*Time, type = 'response')
invasive.means <- tidy(invasive.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

# Native
native <- feeder.veg.sum |>
	filter(Invasive == 'Native')
native.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = native, family = 'binomial')
Anova(native.mod)
check_model(native.mod)
native.means <- emmeans(native.mod, ~Feeder*Time, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Invasive <- 'Native'

far.means <- rbind(invasive.means, native.means)
far.means <- far.means |>
	mutate_if(is.numeric, round, digits = 3)
far.means$Transect <- "Close"
far.means <- far.means |>
	dplyr::select(Transect, Invasive, Time, Feeder, prob, everything())

# Clear the decks
rm(list = ls()[!ls() %in% c("far.means", "feeder.veg")])

pos <- position_dodge(.6)
ggplot(far.means, aes(x = as_factor(Time), y = prob, group = interaction(Feeder,Time),
												fill = Feeder))+
	geom_pointrange(aes(ymin = conf.low, ymax = conf.high), 
                    position = pos, size = 1)+
	geom_point(position = pos, size = 4, shape = 21)+
	scale_fill_manual(values = c('gray', 'skyblue'),
										guide = guide_legend(reverse = TRUE),
										labels = c('Control', 'Feeder'))+
	ylab("Probability")+
	xlab('Year')+
	facet_wrap(~Invasive)+
	ggtitle('Probability plant is invasive and native <10 m')+
	theme_bw()

## --------------- PROPORTION FUNCTIONAL MODEL ---------------------------------

# Create a df with all site/plot/time/feeder combinations
sites <- feeder.veg.sum |>
	group_by(Site, Plot, Time, Feeder, Functional, Transect) |>
	summarize(Everything = n())

# Calculate total invasive hits
total <- feeder.veg.sum |>
	filter(Meter < 11) |>
	filter(Invasive == 'Invasive') |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

# Merge total and sites to force zeroes where invasives missing
total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Tally by functional role
functional <- feeder.veg.sum |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

# Force zeroes
feeder.veg.sum <- merge(functional, total, all = TRUE)
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg.sum", "feeder.veg")])

# Drop natives
feeder.veg.sum <- feeder.veg.sum |> filter(Functional != 'Native')

feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Functional', values_from = 'Count')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0

feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 7:9, names_to = 'Functional', values_to = 'Successes')

# Calculate final dataframe
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Functional, Successes,
								Failures, Total)

# Browsed model
browsed <- feeder.veg.sum |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
browsed.means <- emmeans(browsed.mod, ~Feeder*Time, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

# Ruderal model
ruderals <- feeder.veg.sum |>
	filter(Functional == 'Ruderal')
ruderals.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = ruderals, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
ruderal.means <- emmeans(ruderals.mod, ~Feeder*Time, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderals'

# Resistant model
resistant <- feeder.veg.sum |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
check_model(resistant.mod)
resistant.means <- emmeans(resistant.mod, ~Feeder*Time, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant perennials'

# Create dataframe that summarizes the model results
means <- rbind(browsed.means, resistant.means, ruderal.means)
means <- means |>
	mutate_if(is.numeric, round, digits = 3)
means$Transect <- "Close"
means <- means |>
	select(Transect, Functional, Time, Feeder, prob, everything())

# Clear the decks
rm(list = ls()[!ls() %in% c("means", "feeder.veg")])

# Visualize results
pos <- position_dodge(.6)
ggplot(means, aes(x = as_factor(Time), y = prob, group = interaction(Feeder,Time),
												fill = Feeder))+
	geom_pointrange(aes(ymin = prob-std.error, ymax = prob+std.error), 
                    position = pos, size = 1)+
	geom_point(position = pos, size = 4, shape = 21)+
	scale_fill_manual(values = c('gray', 'skyblue'),
										guide = guide_legend(reverse = TRUE),
										labels = c('Control', 'Feeder'))+
	ylab("Probability")+
	xlab('Year')+
	facet_wrap(~Functional)+
	ggtitle('Probability invasive plant is functional group 1-25 m')+
	theme_bw()

## --------------- CLOSE PROPORTION FUNCTIONAL MODEL ---------------------------

# Subset the data
feeder.veg.sum <- feeder.veg |>
	filter(Meter < 11) #

# Create a df with all site/plot/time/feeder combinations
sites <- feeder.veg.sum |>
	group_by(Site, Plot, Time, Feeder, Functional, Transect) |>
	summarize(Everything = n())

# Calculate total invasive hits
total <- feeder.veg.sum |>
	filter(Meter < 11) |>
	filter(Invasive == 'Invasive') |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

# Merge total and sites to force zeroes where invasives missing
total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Tally by functional role
functional <- feeder.veg.sum |>
	filter(Meter < 11) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

# Force zeroes
feeder.veg.sum <- merge(functional, total, all = TRUE)
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg.sum", "feeder.veg")])

# Drop natives
feeder.veg.sum <- feeder.veg.sum |> filter(Functional != 'Native')

feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Functional', values_from = 'Count')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0

feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 7:9, names_to = 'Functional', values_to = 'Successes')

# Calculate final dataframe
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Functional, Successes,
								Failures, Total)

# Browsed model
browsed <- feeder.veg.sum |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
browsed.means <- emmeans(browsed.mod, ~Feeder*Time, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

# Ruderal model
ruderals <- feeder.veg.sum |>
	filter(Functional == 'Ruderal')
ruderals.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = ruderals, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
ruderal.means <- emmeans(ruderals.mod, ~Feeder*Time, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderals'

# Resistant model
resistant <- feeder.veg.sum |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
check_model(resistant.mod)
resistant.means <- emmeans(resistant.mod, ~Feeder*Time, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant perennials'

# Create dataframe that summarizes the model results
close.means <- rbind(browsed.means, resistant.means, ruderal.means)
close.means <- close.means |>
	mutate_if(is.numeric, round, digits = 3)
close.means$Transect <- "Close"
close.means <- close.means |>
	select(Transect, Functional, Time, Feeder, prob, everything())

# Clear the decks
rm(list = ls()[!ls() %in% c("close.means", "feeder.veg")])

# Visualize results
pos <- position_dodge(.6)
ggplot(close.means, aes(x = as_factor(Time), y = prob, group = interaction(Feeder,Time),
												fill = Feeder))+
	geom_pointrange(aes(ymin = prob-std.error, ymax = prob+std.error), 
                    position = pos, size = 1)+
	geom_point(position = pos, size = 4, shape = 21)+
	scale_fill_manual(values = c('gray', 'skyblue'),
										guide = guide_legend(reverse = TRUE),
										labels = c('Control', 'Feeder'))+
	ylab("Probability")+
	xlab('Year')+
	facet_wrap(~Functional)+
	ggtitle('Probability invasive plant is functional group <10 m')+
	theme_bw()

## --------------- FAR PROPORTION FUNCTIONAL MODEL ---------------------------

# Subset the data
feeder.veg.sum <- feeder.veg |>
	filter(Meter %in% 11:25) #

# Create a df with all site/plot/time/feeder combinations
sites <- feeder.veg.sum |>
	group_by(Site, Plot, Time, Feeder, Functional, Transect) |>
	summarize(Everything = n())

# Calculate total invasive hits
total <- feeder.veg.sum |>
	filter(Meter %in% 11:25) |>
	filter(Invasive == 'Invasive') |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

# Merge total and sites to force zeroes where invasives missing
total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Tally by functional role
functional <- feeder.veg.sum |>
	filter(Meter < 11) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

# Force zeroes
feeder.veg.sum <- merge(functional, total, all = TRUE)
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg.sum", "feeder.veg")])

# Drop natives
feeder.veg.sum <- feeder.veg.sum |> filter(Functional != 'Native')

feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Functional', values_from = 'Count')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0

feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 7:9, names_to = 'Functional', values_to = 'Successes')

# Calculate final dataframe
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Functional, Successes,
								Failures, Total)

# Browsed model
browsed <- feeder.veg.sum |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
browsed.means <- emmeans(browsed.mod, ~Feeder*Time, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

# Ruderal model
ruderals <- feeder.veg.sum |>
	filter(Functional == 'Ruderal')
ruderals.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = ruderals, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
ruderal.means <- emmeans(ruderals.mod, ~Feeder*Time, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderals'

# Resistant model
resistant <- feeder.veg.sum |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
check_model(resistant.mod)
resistant.means <- emmeans(resistant.mod, ~Feeder*Time, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant perennials'

# Create dataframe that summarizes the model results
far.means <- rbind(browsed.means, resistant.means, ruderal.means)
far.means <- far.means |>
	mutate_if(is.numeric, round, digits = 3)
far.means$Transect <- "Close"
far.means <- far.means |>
	select(Transect, Functional, Time, Feeder, prob, everything())

# Clear the decks
rm(list = ls()[!ls() %in% c("far.means", "feeder.veg")])

# Visualize results
pos <- position_dodge(.6)
ggplot(far.means, aes(x = as_factor(Time), y = prob, group = interaction(Feeder,Time),
												fill = Feeder))+
	geom_pointrange(aes(ymin = prob-std.error, ymax = prob+std.error), 
                    position = pos, size = 1)+
	geom_point(position = pos, size = 4, shape = 21)+
	scale_fill_manual(values = c('gray', 'skyblue'),
										guide = guide_legend(reverse = TRUE),
										labels = c('Control', 'Feeder'))+
	ylab("Probability")+
	xlab('Year')+
	facet_wrap(~Functional)+
	ggtitle('Probability invasive plant is functional group 11-25 m')+
	theme_bw()









## --------------- BINARY FULL FUNCTIONAL --------------------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	dplyr::select(Plot, Transect, Meter, Feeder, Functional) |>
	unique() |>
	group_by(Plot, Transect, Feeder, Functional) |>
	summarize(Successes = n())

# Introduce zeroes for missing values
feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Functional', values_from = 'Successes')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0
feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 4:8, names_to = 'Functional', values_to = 'Successes')

# Check values

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
feeder.veg.sum$Total <- 25 # Number of transect points
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total - Successes) |>
	dplyr::select(Plot, Transect, Feeder, Functional, Total, Successes, Failures)

# Ruderal
ruderal <- feeder.veg.sum |>
	filter(Functional == 'Ruderal')

mod.ruderal <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = ruderal)	
Anova(mod.ruderal)
summary(mod.ruderal)

ruderal.means <- emmeans(mod.ruderal, ~Feeder, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderal'

# Resistant perennials
resistant <- feeder.veg.sum |>
	filter(Functional == 'Resistant perennials')

mod.resistant <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = resistant)	
Anova(mod.resistant)
summary(mod.resistant)

resistant.means <- emmeans(mod.resistant, ~Feeder, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant'

# Browsed perennials
browsed <- feeder.veg.sum |>
	filter(Functional == 'Browsed perennials')

mod.browsed <- glm(cbind(Successes, Failures) ~ Feeder + Plot, # Not converging
											family = 'binomial', data = browsed) # with random effects
Anova(mod.browsed)
summary(mod.browsed)

browsed.means <- emmeans(mod.browsed, ~Feeder, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

## --------------- BINARY CLOSE FUNCTIONAL -------------------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	filter(Meter < 11) |>
	dplyr::select(Plot, Transect, Meter, Feeder, Functional) |>
	unique() |>
	group_by(Plot, Transect, Feeder, Functional) |>
	summarize(Successes = n())

# Introduce zeroes for missing values
feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Functional', values_from = 'Successes')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0
feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 4:8, names_to = 'Functional', values_to = 'Successes')

# Check values
v1 <- nrow(feeder.veg.sum) / 5 # functional groups
v2 <- v1 / 2 # feeders
v3 <- v2/4 # transects per plot
isTRUE(length(unique(feeder.veg$Plot)) == v3) # Correct

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
feeder.veg.sum$Total <- 11 # Number of transect points
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total - Successes) |>
	dplyr::select(Plot, Transect, Feeder, Functional, Total, Successes, Failures)

# Ruderal
ruderal <- feeder.veg.sum |>
	filter(Functional == 'Ruderal')
mod.ruderal <- glmer(cbind(Successes, Failures) ~ Feeder + Plot, # No convergence
											family = 'binomial', data = ruderal) # with random effects
Anova(mod.ruderal)
summary(mod.ruderal)
ruderal.means <- emmeans(mod.ruderal, ~Feeder, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderal'

# Resistant perennials
resistant <- feeder.veg.sum |>
	filter(Functional == 'Resistant perennials')

mod.resistant <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = resistant)	
Anova(mod.resistant)
summary(mod.resistant)

resistant.means <- emmeans(mod.resistant, ~Feeder, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant'

# Browsed perennials
browsed <- feeder.veg.sum |>
	filter(Functional == 'Browsed perennials')

mod.browsed <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot),
											family = 'binomial', data = browsed) 
Anova(mod.browsed)
summary(mod.browsed)

browsed.means <- emmeans(mod.browsed, ~Feeder, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

# No plant
non.plant <- feeder.veg.sum |>
	filter(Functional == 'Not plant')

mod.non.plant <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot),
											family = 'binomial', data = non.plant) 
Anova(mod.non.plant)
summary(mod.non.plant)

not.plant.means <- emmeans(mod.non.plant, ~Feeder, type = 'response')
not.plant.means <- tidy(not.plant.means, conf.int = TRUE)
not.plant.means$Functional <- 'Not plant'

## --------------- BINARY FAR FUNCTIONAL ---------------------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	filter(Meter %in% 11:25) |>
	dplyr::select(Plot, Transect, Meter, Feeder, Functional) |>
	unique() |>
	group_by(Plot, Transect, Feeder, Functional) |>
	summarize(Successes = n())

# Introduce zeroes for missing values
feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Functional', values_from = 'Successes')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0
feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 4:8, names_to = 'Functional', values_to = 'Successes')

# Check values
v1 <- nrow(feeder.veg.sum) / 5 # functional groups
v2 <- v1 / 2 # feeders
v3 <- v2/4 # transects per plot
isTRUE(length(unique(feeder.veg$Plot)) == v3) # Correct

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
feeder.veg.sum$Total <- 25 # Number of transect points
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total - Successes) |>
	dplyr::select(Plot, Transect, Feeder, Functional, Total, Successes, Failures)

# Ruderal
ruderal <- feeder.veg.sum |>
	filter(Functional == 'Ruderal')
mod.ruderal <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), # No convergence
											family = 'binomial', data = ruderal) # with random effects
Anova(mod.ruderal)
summary(mod.ruderal)
ruderal.means <- emmeans(mod.ruderal, ~Feeder, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderal'

# Resistant perennials
resistant <- feeder.veg.sum |>
	filter(Functional == 'Resistant perennials')

mod.resistant <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = resistant)	
Anova(mod.resistant)
summary(mod.resistant)

resistant.means <- emmeans(mod.resistant, ~Feeder, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant'

# Browsed perennials
browsed <- feeder.veg.sum |>
	filter(Functional == 'Browsed perennials')

mod.browsed <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot),
											family = 'binomial', data = browsed) 
Anova(mod.browsed)
summary(mod.browsed)

browsed.means <- emmeans(mod.browsed, ~Feeder, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

# No plant
non.plant <- feeder.veg.sum |>
	filter(Functional == 'Not plant')

mod.non.plant <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot),
											family = 'binomial', data = non.plant) 
Anova(mod.non.plant)
summary(mod.non.plant)

not.plant.means <- emmeans(mod.non.plant, ~Feeder, type = 'response')
not.plant.means <- tidy(not.plant.means, conf.int = TRUE)
not.plant.means$Functional <- 'Not plant'









## --------------- LINKING BAREGROUND TO INVASIVE ------------------------------

# find all transect points that are 0 that become a 1 for ruderals
# impact of distance on whether or not that becomes a 1

# same for bare ground?

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Create a df where each transect is a row
meters <- feeder.veg |>
	filter(Feeder == 'Y') |>
	dplyr::select(Site, Plot, Time, Transect) |>
	unique()

# Add a column for each meter on the transect
meters[c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
				'11', '12', '13', '14', '15', '16', '17', '18', '19',
				'20', '21', '22', '23', '24', '25')] <- 0

# Pivot longer
meters <- meters |>
	pivot_longer(cols = 5:29, names_to = 'Distance', values_to = 'Placeholder')

# Filter for the ruderal point-meter hits
ruderal <- feeder.veg |>
	filter(Functional == 'Ruderal') |>
	filter(Feeder == 'Y') |>
	dplyr::select(Site, Plot, Time, Transect, Meter)
colnames(ruderal)[5] <- "Distance"
ruderal$Present <- 1

# Combine with the meter df
ruderal <- merge(meters, ruderal, all = TRUE)
rm(meters)

# Force zeroes
ruderal[is.na(ruderal)] <- 0

# Combine the placeholder 0 and present/absent value into a single binary value
ruderal <- ruderal |>
	mutate(Binary = Placeholder + Present) |>
	dplyr::select(-Present, -Placeholder)

# Make columns the correct type
ruderal$Distance <- as.numeric(ruderal$Distance)
ruderal$Transect <- as_factor(ruderal$Transect)
ruderal$Time <- as_factor(ruderal$Time)

# Filter out the time series
baseline <- ruderal |> filter(Time == '0')
colnames(baseline)[6] <- 'Baseline'
baseline <- baseline |> dplyr::select(-Time)

year1 <- ruderal |> filter(Time == '1')
colnames(year1)[6] <- 'Year1'
year1 <- year1 |> dplyr::select(-Time)

year2 <- ruderal |> filter(Time == '2')
colnames(year2)[6] <- 'Year2'
year2 <- year2 |> dplyr::select(-Time)

# Combine the time seres in wide format
ruderal <- merge(baseline, year1, all = TRUE)
ruderal <- merge(ruderal, year2, all = TRUE)
ruderal <- ruderal |> na.omit()

# Clear the decks
rm(list=ls()[! ls() %in% c("close", "ruderal")])

# Isolate empty transect points that ruderals later colonized
ruderal <- ruderal |>
	filter(Baseline < 1) |>
	mutate(Total = Baseline + Year1 + Year2) |>
	filter(Total > 0) |>
	dplyr::select(-Total)

# Convert to long format
ruderal <- ruderal |>
	dplyr::select(-Baseline) |>
	pivot_longer(cols = 5:6, names_to = "Time", values_to = "Binary")

# Bin the distance values
ruderal$Bin <- NA
for(i in 1:nrow(ruderal)){
	if(ruderal$Distance[i] %in% 1:5){
		ruderal$Bin[i] <- "1-5"
	}
	if(ruderal$Distance[i] %in% 6:10){
		ruderal$Bin[i] <- "6-10"
	}
	if(ruderal$Distance[i] %in% 11:15){
		ruderal$Bin[i] <- "11-15"
	}
	if(ruderal$Distance[i] %in% 16:20){
		ruderal$Bin[i] <- "16-20"
	}
	if(ruderal$Distance[i] %in% 21:25){
		ruderal$Bin[i] <- "21-25"
	}
}

# Drop binary column
ruderal <- ruderal[,-5]

# Drop duplicate values 
# (we should not double count points colonized by and remaining ruderals)
ruderal <- unique(ruderal) 

# Tally colonization events by bin
ruderal.sum <- ruderal |>
	group_by(Plot, Bin) |>
	summarize(Count = sum(Binary))

# Check distribution
descdist(ruderal.sum$Count, discrete = TRUE) # use nonparametric

# Assign ordinality to bin value
ruderal.sum$Bin <- factor(ruderal.sum$Bin, order = TRUE,
													levels = c("1-5", "6-10", "11-15",
																		 "16-20", "21-25"))

# Visualize the data
ggplot(ruderal.sum)+
  aes(x = Bin, y = Count)+
  geom_jitter(width = 0.1, size = 4, shape = 21, fill = 'skyblue',
  						stroke = 1, alpha = 0.5)+
	ylab('Detections')+
	xlab('Distance from feeder (m)')+
  theme_bw()

# Test
kruskal.test(Count ~ Bin, data = ruderal.sum) # p = 0.029

# Repeat exercise at control sites