## --------------- HEADER ------------------------------------------------------
## Script name: 1_Clean-date.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-12
## Date Last modified: 2023-10-04
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for combining and cleaning the 
## manipulative and retrospective feeder vegetation data
 
## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(data.table)

rm(list=ls())

retro <- read.csv('Raw-data/retrospective-transects.csv')
# retro <- retro[,1:9]

manipulative <- fread('Raw-data/manipulative-transects.csv')
manipulative <- manipulative[,c(1:13,17)]

## --------------- COMBINE DATAFRAMES ------------------------------------------

# Treatment column
colnames(retro)[3] <- 'Feeder'
colnames(manipulative)[3] <- 'Feeder'
colnames(retro)[7] <- 'Veg\nType'
colnames(retro)[10] <- 'Browse'
colnames(retro)[11] <- 'Notes'

# Add transect value to retrospective
vals <- paste(retro$Plot, retro$Feeder, retro$Direction)
retro$Transect <- match(vals, unique(vals))
retro <- retro |> 
	mutate("Browse percentage" = 'not measured',
				 "Bite count" = "not measured") |>
	dplyr::select(Plot, Time, Feeder, Direction, Transect, everything())

# Combine data
manipulative$Dataset <- "Manipulative"
retro$Dataset <- "Retrospective"
comb.dat <- rbind(manipulative, retro, use.names=FALSE)

nrow(manipulative) + nrow(retro) # 5042 matches

# Need to add an ID column because there are ~30 instances where the same plant
# was identified at the same point. Merge is erasing this as a duplicate value
comb.dat$ID <- seq(1, 55442, 1)

# Fix direction


## --------------- FIX INVASIVE SPECIES NAMES ----------------------------------

# Make sure invasive and native are correct
for(i in 1:nrow(comb.dat)){
	if(isTRUE(comb.dat$Invasive[i] == 'native')==TRUE){
		comb.dat$Invasive[i] <- 'Native'
	}
	if(isTRUE(comb.dat$Invasive[i] == 'invasive')==TRUE){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

# Check
for(i in 1:nrow(comb.dat)){
	if(isTRUE(comb.dat$Invasive[i] == '')==TRUE){
		comb.dat$Invasive[i] <- 'na'
	}
	if(isTRUE(comb.dat$Invasive[i] == 'X')==TRUE){
		comb.dat$Invasive[i] <- 'na'
	}
	if(isTRUE(comb.dat$Invasive[i] == 'Forb')==TRUE){
		comb.dat$Invasive[i] <- 'Native'
	}
}
unique(comb.dat$Invasive)

# Fix NA values
for(i in 1:nrow(comb.dat)){
	if(isTRUE(comb.dat$Invasive[i] == 'na')==TRUE){
		comb.dat$Invasive[i] <- NA
	}
	
}

# Adjust time value
for(i in 1:nrow(comb.dat)){
	if(isTRUE(comb.dat$Time[i] == '2018')==TRUE){
		comb.dat$Time[i] <- '0' # or 1
	}
	if(isTRUE(comb.dat$Time[i] == '2019')==TRUE){
		comb.dat$Time[i] <- '1' # or 2
	}
	if(isTRUE(comb.dat$Time[i] == '2020')==TRUE){
		comb.dat$Time[i] <- '2' # or 3
	}
	if(isTRUE(comb.dat$Time[i] == '2023')==TRUE){
		comb.dat$Time[i] <- '5' # or 3
	}
}

# Check
unique(comb.dat$Time)

# Set order
comb.dat$Time <- factor(comb.dat$Time, order = TRUE, 
                                    levels = c("0", "1", "2", "3", "4",
                                    					 "5", "7"))
comb.dat |>
	filter(Dataset == 'Retrospective') |>
	group_by(Time) |>
	ggplot(aes(x = Time))+
	geom_histogram(stat = 'count')

comb.dat |>
	filter(Dataset == 'Manipulative') |>
	group_by(Time) |>
	ggplot(aes(x = Time))+
	geom_histogram(stat = 'count')	
	
# Fix feeders
for(i in 1:nrow(comb.dat)){
	if(isTRUE(comb.dat$Feeder[i] == 'F')==TRUE){
		comb.dat$Feeder[i] <- 'Y'
	}
	if(isTRUE(comb.dat$Feeder[i] == 'NF')==TRUE){
		comb.dat$Feeder[i] <- 'N'
	}
}
unique(comb.dat$Feeder)

# Fix the invasive species
for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == 'crabgrass'){
		 comb.dat$Species[i] <- 'Crabgrass'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == 'chinese privet'){
		 comb.dat$Species[i] <- 'Chinese Privet'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == 'japan grass'){
		 comb.dat$Species[i] <- 'Japan Grass'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == 'japan grass '){
		 comb.dat$Species[i] <- 'Japan Grass'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == 'japanese climbing fern'){
		 comb.dat$Species[i] <- 'Japanese Climbing Fern'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == 'japanese honeysuckle'){
		 comb.dat$Species[i] <- 'Japanese Honeysuckle'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == 'sericea lespedeza'){
		 comb.dat$Species[i] <- 'Sericea Lespedeza'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] ==  "johnson grass"){
		 comb.dat$Species[i] <- "Johnson Grass"
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] ==  "brazilian vervain"){
		 comb.dat$Species[i] <- "Brazilian Vervain"  
	}
}

comb.dat$Species <- as.character(comb.dat$Species)
comb.dat$Species[comb.dat$Species=="vasey's grass"]<-"Vasey's Grass"
comb.dat$Species[comb.dat$Species=="Vasey's grass"]<-"Vasey's Grass"
comb.dat$Species[comb.dat$Species=="Vasey'S\xcaGrass"]<-"Vasey's Grass"
comb.dat$Species[comb.dat$Species=="Vasey'S Grass"]<-"Vasey's Grass"
comb.dat$Species[comb.dat$Species=="chinaberry"]<-"Chinaberry"
comb.dat$Species[comb.dat$Species=="dallis grass"]<-"Dallis Grass"
comb.dat$Species[comb.dat$Species=="rye grass"]<-"Rye Grass"
comb.dat$Species[comb.dat$Species=="hemp sesbania"]<-"Hemp Sesbania"
comb.dat$Species[comb.dat$Species=="bahia"]<-"Bahia Grass"
comb.dat$Species[comb.dat$Species=="Bermuda"]<-"Bermudagrass"
comb.dat$Species[comb.dat$Species=="Japan Grass "]<-"Japan Grass"

# Check values
unique(comb.dat$Species[comb.dat$Invasive == 'Invasive'])

unique(comb.dat$Species[comb.dat$Invasive == 'Native'])

# Need to check
# Spanish needles / Bidens
# Dallis Grass 
# Lespedeza
# Thistle
# Ludwigia
# bullrush

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "Dallis Grass"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "White Clover"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "Millet"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "Carpet Weed"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "Barnyard Grass"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "Spider Lily"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "Ivy Leaf Morning Glory"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "False Yam"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "Foxtail"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] == "Purselane"){
		comb.dat$Invasive[i] <- 'Invasive'
	}
}


# Check values
unique(comb.dat$Species[comb.dat$Invasive == 'Invasive'])

# Some of these species are not actually invasive
not.invasive <- c("Hemp Sesbania", 'Jointvetch', 'Dodder')
for(i in 1:nrow(comb.dat)){
	if(comb.dat$Species[i] %in% not.invasive){
		comb.dat$Invasive[i] <- 'Native'
	}
}

## --------------- ADD FUNCTIONAL ROLE TO INVASIVES ----------------------------

# Climbing fern disturbed and undisturbed and resistant to deer

ruderals <- tibble(c("Japan Grass", "Crabgrass", "Japanese Climbing Fern",
							"Johnson Grass", "Rye Grass", "Sickle Pod",
							"Brazilian Vervain",  "Leaf Morning Glory", "Millet",
							"Carpet Weed", "Barnyard Grass", "Bermudagrass", 'Redtop', 'Bahia',
							"Purselane", "Foxtail", "False Yam", "Dallis Grass",
							"Vasey's Grass", "Sericea Lespedeza"))
ruderals$Functional <- 'Ruderal'
colnames(ruderals)[1] <- 'Species'

browsed.perennials <- tibble(c("Japanese Honeysuckle", "Chinese Wisteria",
												"Chinese Privet", "Callery Pear", "Multiflora Rose",
												"White Clover", "Redtop"))
browsed.perennials$Functional <- 'Browsed perennials'
colnames(browsed.perennials)[1] <- 'Species'

resistant.perennials <- tibble(c("Goutweed", "Crepe Myrtle", "Chinaberry", 
																 "Spider Lily"))
resistant.perennials$Functional <- 'Resistant perennials'
colnames(resistant.perennials)[1] <- 'Species'

index <- rbind(ruderals, browsed.perennials)
index <- rbind(index, resistant.perennials)

invasive <- comb.dat |> filter(Invasive == 'Invasive')
unique(invasive$Species)
nrow(ruderals)+nrow(browsed.perennials)+nrow(resistant.perennials)
length(unique(comb.dat$Species[comb.dat$Invasive == 'Invasive']))
# extra is NA

comb.dat <- merge(comb.dat, index, by = 'Species', all.x = TRUE)

comb.dat <- comb.dat |>
	dplyr::select(ID, Plot, Time, Feeder, Direction, Transect, Meter, 'Veg Type',
				 Invasive, Functional, Species, everything())

## --------------- FIX HERBIVORY -----------------------------------------------

rm(list=ls()[! ls() %in% c("comb.dat")])

# There are different schemes for recording browse and a lot of missing values

# Filter out the first experiment
manip.herb <- comb.dat |>
	filter(Dataset == 'Manipulative') 

# Anything with invasive = NA is a non plant, make the value for browse NA
for(i in 1:nrow(manip.herb)){
	if(is.na(manip.herb$Invasive[i])){
		manip.herb$Browse[i] <- NA
	}
}

# Anything blank remaining is likely missing a zero
manip.herb$Browse[manip.herb$Browse==""] <- 0

# Fix all the other weird values
manip.herb.weird.labels <- c('X', '1 of 1', '1o1', '1o4', '3 of 6', '6 of 9')
# Some of these may be leaves or branches eaten?

for(i in 1:nrow(manip.herb)){
	if(manip.herb$Browse[i] %in% manip.herb.weird.labels){
		manip.herb$Browse[i] <- 1
	}
}

# Fix the retrospective
retro.herb <- comb.dat |>
	filter(Dataset == 'Retrospective') 

retro.herb$Browse[retro.herb$Browse==""] <- 0
retro.herb$Browse[retro.herb$Browse=="X"] <- 1
retro.herb$Browse[retro.herb$Browse=="x"] <- 1
retro.herb$Browse[retro.herb$Browse=="y"] <- 1

# Recombine the data again
comb.dat <- rbind(manip.herb, retro.herb)

## --------------- ADD A COLUMN FOR SITE ---------------------------------------

comb.dat$Site <- NA

manip.sites <- comb.dat |> filter(Dataset == 'Manipulative')
length(unique(manip.sites$Plot))

D <- c("D1", "D2", "D3")
LV <- c("LV1", "LV2", "LV3", "LV4", "LV5", "LV6")
P <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")
WS <- c("WS1", "WS2", "WS3", "WS4", "WS5", "WS6", "WS7", "WS8")

for(i in 1:nrow(comb.dat)){
	if(comb.dat$Plot[i] %in% D){
		comb.dat$Site[i] <- 'D'
	}
	if(comb.dat$Plot[i] %in% LV){
		comb.dat$Site[i] <- 'LV'
	}
	if(comb.dat$Plot[i] %in% P){
		comb.dat$Site[i] <- 'P'
	}
	if(comb.dat$Plot[i] %in% WS){
		comb.dat$Site[i] <- 'WS'
	}
}

comb.dat <- comb.dat |> dplyr::select(ID, Site, everything())

## --------------- LAST DETAILS ------------------------------------------------

# Need to fix a lowercase w
comb.dat$Direction <- toupper(comb.dat$Direction)

## --------------- EXPORT THE DATA ---------------------------------------------

# Export
write.csv(comb.dat,
					'Clean-data/Feeder-veg.csv', row.names = FALSE)

