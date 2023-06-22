## --------------- HEADER ------------------------------------------------------
## Script name: 2_EDA.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-12
## Date Last modified: 2023-06-15
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for performing EDA on the clean feeder
## vegetation data

## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse)
library(DataExplorer)

# Clear the decks
rm(list=ls())

# Bring in the data
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv')

invasive <- feeder.veg |> filter(Invasive == 'Invasive')
unique(invasive$Species)

plot_str(feeder.veg, type = 'r') # Structure
plot_intro(feeder.veg) # Introduction
plot_missing(feeder.veg) # NAs from bare ground obs. in the invasive column
plot_bar(feeder.veg) # frequency distributions of discrete features
plot_histogram(feeder.veg)

ggplot(feeder.veg, aes(x = Meter)) + 
  geom_histogram(colour = 4, fill = "white", bins = 50)
# Transect points switched from 4 x 50 m to 8 x 25 m

# Even number of transects?
hits.per.transect <- feeder.veg |>
	group_by(Plot, Time, Transect, Feeder) |>
	dplyr::summarize(n=n())
# WS2 had extreme value (170+) because two transects had the same label
# This was rectified
# Other high values are accounted for in the data

transects.per.plot <- hits.per.transect |>
	group_by(Plot) |>
	dplyr::summarize(n=n())
# Sites have too many and not enough transects

transects.per.plot.yr <- hits.per.transect |>
	group_by(Plot, Time) |>
	dplyr::summarize(n=n())
# LV3 LV4 WS8 missing two transects
# WS6 missing one
# WS2 has an extra

transects.per.plot.yr.feeder <- feeder.veg |>
	group_by(Plot, Time, Feeder) |>
	dplyr::summarize(n=n_distinct(Transect))
# WS2 has an extra at no feeder
# WS6 missing one at no feeder
# LV3 LV4 No feeder in 2019
# WS8 missing one at each treatment in 2018

## --------------- VISUALIZE BAREGROUND ----------------------------------------
rm(list=ls()[! ls() %in% c("feeder.veg")])

# Get proportion bare ground at each meter regardless of plot and transect
meter.sum <- feeder.veg |>
									filter(Feeder == 'Y') |>
									group_by(Meter, Invasive) |>
									dplyr::summarize(Count = n())

bare.ground <- meter.sum |>
	filter(is.na(Invasive))
colnames(bare.ground)[2] <- 'Plant'
bare.ground$Plant <- 'N'

plants <- meter.sum |>
	filter(!is.na(Invasive)) |>
	group_by(Meter) |>
	dplyr::summarize(Count = sum(Count))
plants$Plant <- 'Y'
plants <- plants |>
	dplyr::select(Meter, Plant, Count)

bare.ground <- rbind(bare.ground,plants)
bare.ground <- bare.ground |>
	pivot_wider(names_from = Plant, values_from = Count) |>
	mutate(Bare.prop = N/(N+Y))

ggplot(bare.ground, aes(x = Meter, y = Bare.prop))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))

# Do it with time
meter.sum <- feeder.veg |>
									filter(Feeder == 'Y') |>
									group_by(Time, Meter, Invasive) |>
									dplyr::summarize(Count = n())

bare.ground <- meter.sum |>
	filter(is.na(Invasive))
colnames(bare.ground)[3] <- 'Plant'
bare.ground$Plant <- 'N'

plants <- meter.sum |>
	filter(!is.na(Invasive)) |>
	group_by(Time, Meter) |>
	dplyr::summarize(Count = sum(Count))
plants$Plant <- 'Y'
plants <- plants |>
	dplyr::select(Time, Meter, Plant, Count)

bare.ground <- rbind(bare.ground,plants)
bare.ground <- bare.ground |>
	pivot_wider(names_from = Plant, values_from = Count) |>
	mutate(Bare.prop = N/(N+Y))

ggplot(bare.ground, aes(x = Meter, y = Bare.prop))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Time)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))


## --------------- VISUALIZE HERBIVORY -----------------------------------------
rm(list=ls()[! ls() %in% c("feeder.veg")])

# Make empties NA
# feeder.veg$Browse[feeder.veg$Browse==""]<-NA

browse.fig <- feeder.veg |>
	filter(!is.na(Browse))
sum(browse.fig$Browse==1)/(sum(browse.fig$Browse==1)+sum(browse.fig$Browse==0))
# Less than 1 percent of plants are browsed

total.plants <- browse.fig |>
	group_by(Feeder, Meter) |>
	summarize(Total = n())

total.browsed <- feeder.veg |>
	filter(Browse == '1') |>
	group_by(Feeder, Meter) |>
	summarize(Browsed = n())

browse.fig <- merge(total.browsed, total.plants, all = TRUE)
browse.fig[is.na(browse.fig)] <- 0

browse.fig <- browse.fig |>
	mutate(Browse.proportion = Browsed/Total)

ggplot(browse.fig, aes(x = Meter, y = Browse.proportion, fill = Feeder))+
	geom_col(position = 'dodge')+
	scale_x_continuous(breaks = seq(0,50,1))+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

## --------------- VISUALIZE CLOSE PLANTS --------------------------------------
rm(list=ls()[! ls() %in% c("feeder.veg")])

# Subset the data
feeder.veg.close <- feeder.veg |>
	filter(Meter < 6)

# Total plants
tot.plants <- feeder.veg.close |>
	na.omit() |>
	group_by(Time) |>
	summarise(Total = n())

plants <- feeder.veg.close |>
	na.omit() |>
	group_by(Time, Feeder) |>
	summarise(Plants = n())

plants <- merge(plants, tot.plants, all = TRUE) |>
	mutate(Proportion = Plants)

# Invasive total
tot.invasive <- feeder.veg.close |>
	filter(!is.na(Invasive)) |>
	group_by(Feeder) |>
	summarise(Total = n())

invasive <- feeder.veg.close |>
	group_by(Feeder, Invasive) |>
	summarise(Plants = n()) |>
	filter(!is.na(Invasive))

invasive <- merge(invasive, tot.invasive, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(invasive, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Invasive)

ggplot(invasive, aes(x = Invasive, y = Proportion))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Feeder)

# Invasive x time
tot.invasive <- feeder.veg.close |>
	filter(!is.na(Invasive)) |>
	group_by(Time, Feeder) |>
	summarise(Total = n())

invasive <- feeder.veg.close |>
	group_by(Time, Feeder, Invasive) |>
	summarise(Plants = n()) |>
	filter(!is.na(Invasive))

invasive <- merge(invasive, tot.invasive, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(invasive, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Invasive)

ggplot(invasive, aes(x = Time, y = Proportion, fill = Invasive))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Feeder)

# Functional
tot.functional <- feeder.veg.close |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder) |>
	summarise(Total = n())

functional <- feeder.veg.close |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder, Functional) |>
	summarise(Plants = n()) |>
	na.omit()

functional <- merge(functional, tot.functional, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(functional, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Functional)

ggplot(functional, aes(x = Functional, y = Proportion))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Feeder)

# Functional x time
tot.functional <- feeder.veg.close |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder) |>
	summarise(Total = n())

functional <- feeder.veg.close |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder, Functional) |>
	summarise(Plants = n()) |>
	na.omit()

functional <- merge(functional, tot.functional, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(functional, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Functional)

ggplot(functional, aes(x = Time, y = Proportion, fill = Functional))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Feeder)

# Sites that have each functional group
ruderal <- feeder.veg.close |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n()) |>
	filter(Functional == 'Ruderal')


# Species
species <- feeder.veg.close |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder, Species) |>
	summarise(Plants = n())

ggplot(species, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Species)

# Species x time
species <- feeder.veg.close |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder, Species) |>
	summarise(Plants = n())

ggplot(species, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Species)

## --------------- VISUALIZE MID PLANTS ----------------------------------------

rm(list=ls()[! ls() %in% c("feeder.veg")])

# Subset the data
feeder.veg.mid <- feeder.veg |>
	filter(Meter %in% c(5:10))

# Total plants
tot.plants <- feeder.veg.mid |>
	na.omit() |>
	group_by(Time) |>
	summarise(Total = n())

plants <- feeder.veg.mid |>
	na.omit() |>
	group_by(Time, Feeder) |>
	summarise(Plants = n())

plants <- merge(plants, tot.plants, all = TRUE) |>
	mutate(Proportion = Plants)

# Invasive total
tot.invasive <- feeder.veg.mid |>
	filter(!is.na(Invasive)) |>
	group_by(Feeder) |>
	summarise(Total = n())

invasive <- feeder.veg.mid |>
	group_by(Feeder, Invasive) |>
	summarise(Plants = n()) |>
	filter(!is.na(Invasive))

invasive <- merge(invasive, tot.invasive, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(invasive, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Invasive)

ggplot(invasive, aes(x = Invasive, y = Proportion))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Feeder)

# Invasive x time
tot.invasive <- feeder.veg.mid |>
	filter(!is.na(Invasive)) |>
	group_by(Time, Feeder) |>
	summarise(Total = n())

invasive <- feeder.veg.mid |>
	group_by(Time, Feeder, Invasive) |>
	summarise(Plants = n()) |>
	filter(!is.na(Invasive))

invasive <- merge(invasive, tot.invasive, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(invasive, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Invasive)

ggplot(invasive, aes(x = Time, y = Proportion, fill = Invasive))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Feeder)

# Functional
tot.functional <- feeder.veg.mid |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder) |>
	summarise(Total = n())

functional <- feeder.veg.mid |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder, Functional) |>
	summarise(Plants = n()) |>
	na.omit()

functional <- merge(functional, tot.functional, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(functional, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Functional)

ggplot(functional, aes(x = Functional, y = Proportion))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Feeder)

# Functional x time
tot.functional <- feeder.veg.mid |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder) |>
	summarise(Total = n())

functional <- feeder.veg.mid |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder, Functional) |>
	summarise(Plants = n()) |>
	na.omit()

functional <- merge(functional, tot.functional, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(functional, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Functional)

ggplot(functional, aes(x = Time, y = Proportion, fill = Functional))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Feeder)

# Species
species <- feeder.veg.mid |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder, Species) |>
	summarise(Plants = n())

ggplot(species, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Species)

# Species x time
species <- feeder.veg.mid |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder, Species) |>
	summarise(Plants = n())

ggplot(species, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Species)

## --------------- VISUALIZE FAR PLANTS ----------------------------------------

rm(list=ls()[! ls() %in% c("feeder.veg")])

# Subset the data
feeder.veg.far <- feeder.veg |>
	filter(Meter %in% c(15:25))

# Total plants
tot.plants <- feeder.veg.far |>
	na.omit() |>
	group_by(Time) |>
	summarise(Total = n())

plants <- feeder.veg.far |>
	na.omit() |>
	group_by(Time, Feeder) |>
	summarise(Plants = n())

plants <- merge(plants, tot.plants, all = TRUE) |>
	mutate(Proportion = Plants)

# Invasive total
tot.invasive <- feeder.veg.far |>
	filter(!is.na(Invasive)) |>
	group_by(Feeder) |>
	summarise(Total = n())

invasive <- feeder.veg.far |>
	group_by(Feeder, Invasive) |>
	summarise(Plants = n()) |>
	filter(!is.na(Invasive))

invasive <- merge(invasive, tot.invasive, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(invasive, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Invasive)

ggplot(invasive, aes(x = Invasive, y = Proportion))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Feeder)

# Invasive x time
tot.invasive <- feeder.veg.far |>
	filter(!is.na(Invasive)) |>
	group_by(Time, Feeder) |>
	summarise(Total = n())

invasive <- feeder.veg.far |>
	group_by(Time, Feeder, Invasive) |>
	summarise(Plants = n()) |>
	filter(!is.na(Invasive))

invasive <- merge(invasive, tot.invasive, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(invasive, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Invasive)

ggplot(invasive, aes(x = Time, y = Proportion, fill = Invasive))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Feeder)

# Functional
tot.functional <- feeder.veg.far |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder) |>
	summarise(Total = n())

functional <- feeder.veg.far |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder, Functional) |>
	summarise(Plants = n()) |>
	na.omit()

functional <- merge(functional, tot.functional, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(functional, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Functional)

ggplot(functional, aes(x = Functional, y = Proportion))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Feeder)

# Functional x time
tot.functional <- feeder.veg.far |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder) |>
	summarise(Total = n())

functional <- feeder.veg.far |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder, Functional) |>
	summarise(Plants = n()) |>
	na.omit()

functional <- merge(functional, tot.functional, all = TRUE) |>
	mutate(Proportion = Plants/Total)

ggplot(functional, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Functional)

ggplot(functional, aes(x = Time, y = Proportion, fill = Functional))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Feeder)

# Species
species <- feeder.veg.far |>
	filter(Invasive == 'Invasive') |>
	group_by(Feeder, Species) |>
	summarise(Plants = n())

ggplot(species, aes(x = Feeder, y = Plants))+
	geom_bar(stat = 'identity')+
	facet_wrap(~Species)

# Species x time
species <- feeder.veg.far |>
	filter(Invasive == 'Invasive') |>
	group_by(Time, Feeder, Species) |>
	summarise(Plants = n())

ggplot(species, aes(x = Time, y = Plants, fill = Feeder))+
	geom_bar(stat = 'identity', position = 'dodge')+
	facet_wrap(~Species)

## --------------- MARCUS PLANTS -----------------------------------------------

plants <- feeder.veg |>
	filter(!is.na(Invasive)) |>
	filter(Dataset == 'Manipulative')
34679/54042

# plants <- plants |>
	# filter(Time == 1)

total.pair <- plants |>
	group_by(Plot, Feeder, Time) |>
	summarize(Count = n()) |>
	pivot_wider(names_from = Time, values_from = Count)

summary)to