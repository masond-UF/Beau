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

ggplot(feeder.veg, aes(x = Meter)) + 
  geom_histogram(colour = 4, fill = "white", bins = 50)
# Transect points switched from 4 x 50 m to 8 x 25 m

# Even number of transects?
feeder.veg$Direction <- toupper(feeder.veg$Direction)
hits.per.transect <- feeder.veg |>
	group_by(Plot, Time, Direction, Feeder) |>
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
	dplyr::summarize(n=n_distinct(Direction))
# WS2 has an extra at no feeder
# WS6 missing one at no feeder
# LV3 LV4 No feeder in 2019
# WS8 missing one at each treatment in 2018

## --------------- VISUALIZE OVERALL -------------------------------------------
rm(list=ls()[! ls() %in% c("feeder.veg")])
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv')

for(i in 1:nrow(feeder.veg)){
	if(isTRUE(feeder.veg$Species[i] == 'Bare Ground')){
		feeder.veg$Functional[i] <- 'Bare Ground'
	}
	if(isTRUE(feeder.veg$Invasive[i] == 'Native')){
		feeder.veg$Functional[i] <- 'Native'
	}
}

# Get proportion bare ground at each meter regardless of plot and transect
meter.sum <- feeder.veg |>
									group_by(Feeder, Meter, Functional) |>
									dplyr::summarize(Count = n())

meter.sum <- meter.sum |>
	pivot_wider(names_from = Functional, values_from = Count) |>
	select(Feeder, Meter, Native, everything())

meter.sum[is.na(meter.sum)] <- 0

colnames(meter.sum)[8] <- 'Other'

meter.sum <- meter.sum |>
	mutate(Total.pts = Native + `Bare Ground` + `Browsed perennials` +
				 	`Resistant perennials` + Ruderal + Other,
				 Plant.pts = Native + `Bare Ground` + `Browsed perennials` +
				 	`Resistant perennials` + Ruderal,
				 Per.bg = `Bare Ground`/Total.pts,
				 Per.browsed = `Browsed perennials`/Plant.pts,
				 Per.ruderal = Ruderal/Plant.pts,
				 Per.resist = `Resistant perennials`/Plant.pts)

# Percent bare ground
ggplot(meter.sum, aes(x = Meter, y = Per.bg))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Time*Feeder, nrow = 4)

# Percent browsed
ggplot(meter.sum, aes(x = Meter, y = Per.browsed))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Time*Feeder, nrow = 4)

# Percent ruderal
ggplot(meter.sum, aes(x = Meter, y = Per.ruderal))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Time*Feeder, nrow = 4)

# Percent resistance
ggplot(meter.sum, aes(x = Meter, y = Per.resist))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Time*Feeder, nrow = 4)

## --------------- VISUALIZE MANIPULATIVE --------------------------------------

rm(list=ls()[! ls() %in% c("feeder.veg")])
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative')

for(i in 1:nrow(feeder.veg)){
	if(isTRUE(feeder.veg$Species[i] == 'Bare Ground')){
		feeder.veg$Functional[i] <- 'Bare Ground'
	}
	if(isTRUE(feeder.veg$Invasive[i] == 'Native')){
		feeder.veg$Functional[i] <- 'Native'
	}
}

# Get proportion bare ground at each meter regardless of plot and transect
meter.sum <- feeder.veg |>
									group_by(Time, Feeder, Meter, Functional) |>
									dplyr::summarize(Count = n())

meter.sum <- meter.sum |>
	pivot_wider(names_from = Functional, values_from = Count) |>
	select(Feeder, Meter, Native, everything())

meter.sum[is.na(meter.sum)] <- 0

colnames(meter.sum)[8] <- 'Other'

meter.sum <- meter.sum |>
	mutate(Total.pts = Native + `Bare Ground` + `Browsed perennials` +
				 	`Resistant perennials` + Ruderal + Other,
				 Plant.pts = Native + `Bare Ground` + `Browsed perennials` +
				 	`Resistant perennials` + Ruderal,
				 Per.bg = `Bare Ground`/Total.pts,
				 Per.browsed = `Browsed perennials`/Plant.pts,
				 Per.ruderal = Ruderal/Plant.pts,
				 Per.resist = `Resistant perennials`/Plant.pts)

# Percent bare ground
ggplot(meter.sum, aes(x = Meter, y = Per.bg))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Feeder)

# Percent browsed
ggplot(meter.sum, aes(x = Meter, y = Per.browsed))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Feeder)

# Percent ruderal
ggplot(meter.sum, aes(x = Meter, y = Per.ruderal))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Feeder)

## --------------- VISUALIZE RETROSPECTIVE -------------------------------------

rm(list=ls()[! ls() %in% c("feeder.veg")])
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Retrospective')

for(i in 1:nrow(feeder.veg)){
	if(isTRUE(feeder.veg$Species[i] == 'Bare Ground')){
		feeder.veg$Functional[i] <- 'Bare Ground'
	}
	if(isTRUE(feeder.veg$Invasive[i] == 'Native')){
		feeder.veg$Functional[i] <- 'Native'
	}
}

# Get proportion bare ground at each meter regardless of plot and transect
meter.sum <- feeder.veg |>
									group_by(Feeder, Meter, Functional) |>
									dplyr::summarize(Count = n())

meter.sum <- meter.sum |>
	pivot_wider(names_from = Functional, values_from = Count) |>
	select(Feeder, Meter, Native, everything())

meter.sum[is.na(meter.sum)] <- 0

colnames(meter.sum)[7] <- 'Other'

meter.sum <- meter.sum |>
	mutate(Total.pts = Native + `Bare Ground` + `Browsed perennials` +
				 	`Resistant perennials` + Ruderal + Other,
				 Plant.pts = Native + `Bare Ground` + `Browsed perennials` +
				 	`Resistant perennials` + Ruderal,
				 Per.bg = `Bare Ground`/Total.pts,
				 Per.browsed = `Browsed perennials`/Plant.pts,
				 Per.ruderal = Ruderal/Plant.pts,
				 Per.resist = `Resistant perennials`/Plant.pts)

# Percent bare ground
ggplot(meter.sum, aes(x = Meter, y = Per.bg))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Feeder)

# Percent browsed
ggplot(meter.sum, aes(x = Meter, y = Per.browsed))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Feeder)

# Percent ruderal
ggplot(meter.sum, aes(x = Meter, y = Per.ruderal))+
	geom_col()+
	scale_x_continuous(breaks = seq(1, 50, 1))+
	facet_wrap(~Feeder)

## --------------- VISUALIZE HERBIVORY -----------------------------------------
rm(list=ls()[! ls() %in% c("feeder.veg")])
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative')

# Make empties NA
# feeder.veg$Browse[feeder.veg$Browse==""]<-NA

browse.fig <- feeder.veg |>
	filter(!is.na(Browse))
sum(browse.fig$Browse==1)/(sum(browse.fig$Browse==1)+sum(browse.fig$Browse==0))
# Less than 1 percent of plants are browsed

total.plants <- browse.fig |>
	filter(Veg.Type != 'na') |>
	group_by(Feeder, Meter) |>
	summarize(Total = n())

total.browsed <- feeder.veg |>
	filter(Browse == '1') |>
	group_by(Time,Feeder, Meter) |>
	summarize(Browsed = n())

browse.fig <- merge(total.browsed, total.plants, all = TRUE)
browse.fig[is.na(browse.fig)] <- 0

browse.fig <- browse.fig |>
	mutate(Browse.proportion = Browsed/Total)

ggplot(browse.fig, aes(x = Meter, y = Browse.proportion, fill = Feeder))+
	geom_col(position = 'dodge')+
	scale_x_continuous(breaks = seq(0,50,1))+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Overall higher with feeder 

ggplot(browse.fig, aes(x = Meter, y = Browse.proportion, fill = Feeder))+
	geom_col(position = 'dodge')+
	scale_x_continuous(breaks = seq(0,50,1))+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))+
	facet_wrap(~Time)
# Make sure to drop 0 and 5

## --------------- TAKE HOME MESSAGE -------------------------------------------

# feeders increase soil disturbance (small zone of impact)
# feeders increase browse overall (no gradient)
# as a result, feeders are having dichotomous impacts on feeders
## --------------- CHECK PROPORTION OF BAREGROUND AFTER RETURN -----------------

return.sample <- feeder.veg |>
	filter(Dataset == 'Manipulative' &
				 	Time == 5 &
				 	Meter < 4) 

sites <- return.sample |>
	dplyr::select(Plot, Direction) |>
	unique()

95*3 # transects x points (1, 2, 3 m) = 285 points

bareground <- return.sample |>
	group_by(Plot, Direction, Species) |>
	summarize(Count = n()) |>
	filter(Species %in% c('Bare Ground', 'Litter'))

sum(bareground$Count) # 117 

117/285 # 41 %

## --------------- OLD VISUALIZE CLOSE PLANTS --------------------------------------
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

## --------------- OLD VISUALIZE MID PLANTS ----------------------------------------

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

## --------------- OLD VISUALIZE FAR PLANTS ----------------------------------------

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

## --------------- OLD MARCUS PLANTS -----------------------------------------------

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