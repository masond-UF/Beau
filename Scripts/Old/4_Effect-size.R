## --------------- HEADER ------------------------------------------------------
## Script name: 4_Effect-size.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-14
## Date Last modified: 2023-06-15
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for generating effect sizes

## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse) # Everything
library(fitdistrplus) # Checking distribution
library(pwr) # Calculating Cohen's H
library(lme4) # Conducting mixed effects linear model

# Clear the decks
rm(list=ls())

# Bring in the manipulative experiment data
d <- read.csv('Clean-data/Feeder-veg.csv') |> 
	filter(Dataset == 'Manipulative')

## --------------- OVERALL INVASIVES -------------------------------------------

# Tally natives and invasives at each plot*time
total <- d |>
	filter(Meter < 26) |>
	filter(!is.na(Invasive)) |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

manip.invasive <- d |>
	filter(Meter < 26) |>
	filter(!is.na(Invasive)) |>
	group_by(Plot, Time, Feeder, Invasive) |>
	summarize(Count = n())

manip.invasive <- manip.invasive |>
	pivot_wider(names_from = Invasive, values_from = Count)
manip.invasive[is.na(manip.invasive)] <- 0

manip.invasive <- manip.invasive |>
	pivot_longer(cols = 4:5, names_to = 'Invasive', values_to = 'Count')

# Calculate proportion
manip <- merge(manip.invasive, total, all = TRUE)
manip <- manip |>
	mutate(Proportion = Count/Total)

# Clear the decks
rm(total, manip.invasive)

# Subset by year
baseline <- manip |>
	filter(Time == '0')
baseline <- baseline[, -c(2,5:6)]
colnames(baseline)[4] <- 'Baseline.proportion'

year1 <- manip |>
	filter(Time == '1')
year1 <- year1[, -c(2,5:6)]
colnames(year1)[4] <- 'Year1.proportion'

year2 <- manip |>
	filter(Time == '2')
year2 <- year2[, -c(2,5:6)]
colnames(year2)[4] <- 'Year2.proportion'

# Merge
manip <- merge(year2, year1, all.x = TRUE)
manip <- merge(manip, baseline, all.x = TRUE)

# Clear the decks
rm(baseline, year1, year2)

# Introduce zeroes
manip[is.na(manip)] <- 0

# Calculate the delta from baseline to subsequent years
manip <- manip |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

manip.sites <- manip |>
	dplyr::select(Plot, Feeder, Invasive, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H") |>
	filter(Invasive != 'Native')

# Visualize by sites
ggplot(manip.sites, aes(x = Year, y = cohen.H, color = Feeder))+
	geom_point()+
	facet_wrap(~Plot)
	

descdist(manip.sites$cohen.H, discrete = FALSE)

m1 <- lmer(cohen.H~Feeder*Year + (1|Plot), data = manip.sites)
Anova(m1, type = 3)

cohen.H.df <- manip.sites |>
	group_by(Year, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))

control <- manip |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- manip |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	filter(Invasive == 'Invasive')

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Delta.1, Delta.2) |>
	pivot_longer(2:3, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()

## --------------- FUNCTIONAL --------------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d")])

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter < 26) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter < 26) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, test, sites)

# Subset by year
baseline <- functional |>
	filter(Time == '0')
baseline <- baseline[, -c(2,5:6)]
colnames(baseline)[4] <- 'Baseline.proportion'

year1 <- functional |>
	filter(Time == '1')
year1 <- year1[, -c(2,5:6)]
colnames(year1)[4] <- 'Year1.proportion'

year2 <- functional |>
	filter(Time == '2')
year2 <- year2[, -c(2,5:6)]
colnames(year2)[4] <- 'Year2.proportion'

# Merge
functional.wd <- merge(year2, year1, all.x = TRUE)
functional.wd <- merge(functional.wd, baseline, all.x = TRUE)

# Clear the decks
rm(baseline, year1, year2)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Off the map. Not good.

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- FUNCTIONAL NO ZEROES BAD ------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d")])

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter < 26) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter < 26) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Drop the zeroes for each [BAD WAY OF DOING IT]
browsed <- functional |>
	filter(Functional == 'Browsed perennials')
browsed.ls <- browsed |>
	filter(Time == 0 & Proportion == 0)
browsed.ls <- unique(browsed.ls$Plot)
browsed <- browsed |>
	filter(!Plot %in% browsed.ls)

ruderal <- functional |>
	filter(Functional == 'Ruderal')
ruderal.ls <- ruderal |>
	filter(Time == 0 & Proportion == 0)
ruderal.ls <- unique(ruderal.ls$Plot)
ruderal <- ruderal |>
	filter(!Plot %in% ruderal.ls)

resistant <- functional |>
	filter(Functional == 'Resistant perennials')
resistant.ls <- resistant |>
	filter(Time == 0 & Proportion == 0)
resistant.ls <- unique(resistant.ls$Plot)
resistant <- resistant |>
	filter(!Plot %in% resistant.ls)

# Subset by year browsed
browsed.baseline <- browsed |>
	filter(Time == '0')
browsed.baseline <- browsed.baseline[, -c(2,5:6)]
colnames(browsed.baseline)[4] <- 'Baseline.proportion'

browsed.year1 <- browsed |>
	filter(Time == '1')
browsed.year1 <- browsed.year1[, -c(2,5:6)]
colnames(browsed.year1)[4] <- 'Year1.proportion'

browsed.year2 <- browsed |>
	filter(Time == '2')
browsed.year2 <- browsed.year2[, -c(2,5:6)]
colnames(browsed.year2)[4] <- 'Year2.proportion'

# Merge
browsed.wd <- merge(browsed.year2, browsed.year1, all.x = TRUE)
browsed.wd <- merge(browsed.wd, browsed.baseline, all.x = TRUE)

# Subset by year ruderal
ruderal.baseline <- ruderal |>
	filter(Time == '0')
ruderal.baseline <- ruderal.baseline[, -c(2,5:6)]
colnames(ruderal.baseline)[4] <- 'Baseline.proportion'

ruderal.year1 <- ruderal |>
	filter(Time == '1')
ruderal.year1 <- ruderal.year1[, -c(2,5:6)]
colnames(ruderal.year1)[4] <- 'Year1.proportion'

ruderal.year2 <- ruderal |>
	filter(Time == '2')
ruderal.year2 <- ruderal.year2[, -c(2,5:6)]
colnames(ruderal.year2)[4] <- 'Year2.proportion'

# Merge
ruderal.wd <- merge(ruderal.year2, ruderal.year1, all.x = TRUE)
ruderal.wd <- merge(ruderal.wd, ruderal.baseline, all.x = TRUE)

# Subset by year resistant
resistant.baseline <- resistant |>
	filter(Time == '0')
resistant.baseline <- resistant.baseline[, -c(2,5:6)]
colnames(resistant.baseline)[4] <- 'Baseline.proportion'

resistant.year1 <- resistant |>
	filter(Time == '1')
resistant.year1 <- resistant.year1[, -c(2,5:6)]
colnames(resistant.year1)[4] <- 'Year1.proportion'

resistant.year2 <- resistant |>
	filter(Time == '2')
resistant.year2 <- resistant.year2[, -c(2,5:6)]
colnames(resistant.year2)[4] <- 'Year2.proportion'

# Merge
resistant.wd <- merge(resistant.year2, resistant.year1, all.x = TRUE)
resistant.wd <- merge(resistant.wd, resistant.baseline, all.x = TRUE)

# Combine
functional.wd <- rbind(browsed.wd, ruderal.wd, resistant.wd)

# Clear the decks
rm(browsed.baseline, browsed.year1, browsed.year2, browsed.wd,
	 ruderal.baseline, ruderal.year1, ruderal.year2, ruderal.wd,
	 resistant.baseline, resistant.year1, resistant.year2, resistant.wd,
	 browsed, resistant, ruderal)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Baseline.proportion, Year1.proportion),
				 cohen.H.year2 = ES.h(Baseline.proportion, Year2.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Off the map. Not good.

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)

## --------------- OVERALL INVASIVES WHEN PRESENT ------------------------------

rm(list=ls()[! ls() %in% c("d")])

# Tally natives and invasives at each plot*time
total <- d |>
	filter(Meter < 26) |>
	filter(!is.na(Invasive)) |>
	group_by(Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

manip.invasive <- d |>
	filter(Meter < 26) |>
	filter(!is.na(Invasive)) |>
	group_by(Plot, Time, Feeder, Invasive, Transect) |>
	summarize(Count = n())

manip.invasive <- manip.invasive |>
	pivot_wider(names_from = Invasive, values_from = Count)
manip.invasive[is.na(manip.invasive)] <- 0

manip.invasive <- manip.invasive |>
	pivot_longer(cols = 5:6, names_to = 'Invasive', values_to = 'Count')

# Calculate proportion
manip <- merge(manip.invasive, total, all = TRUE)
manip <- manip |>
	mutate(Proportion = Count/Total)

# Clear the decks
rm(total, manip.invasive)

# Subset by year
baseline <- manip |>
	filter(Time == '0')
baseline <- baseline[, -c(2,6:7)]
colnames(baseline)[5] <- 'Baseline.proportion'

year1 <- manip |>
	filter(Time == '1')
year1 <- year1[, -c(2,6:7)]
colnames(year1)[5] <- 'Year1.proportion'

year2 <- manip |>
	filter(Time == '2')
year2 <- year2[, -c(2,6:7)]
colnames(year2)[5] <- 'Year2.proportion'

# Merge
manip <- merge(year2, year1, all.x = TRUE)
manip <- merge(manip, baseline, all.x = TRUE)

# Clear the decks
rm(baseline, year1, year2)

# Introduce zeroes
manip[is.na(manip)] <- 0

# Pivot wider 
manip <- manip |>
	pivot_longer(cols = 5:7, names_to = 'Time', values_to = 'Proportion')

## --------------- FUNCTIONAL <5 m ---------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d")])

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter < 6) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter < 6) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Subset by year
baseline <- functional |>
	filter(Time == '0')
baseline <- baseline[, -c(2,5:6)]
colnames(baseline)[4] <- 'Baseline.proportion'

year1 <- functional |>
	filter(Time == '1')
year1 <- year1[, -c(2,5:6)]
colnames(year1)[4] <- 'Year1.proportion'

year2 <- functional |>
	filter(Time == '2')
year2 <- year2[, -c(2,5:6)]
colnames(year2)[4] <- 'Year2.proportion'

# Merge
functional.wd <- merge(year2, year1, all.x = TRUE)
functional.wd <- merge(functional.wd, baseline, all.x = TRUE)

# Clear the decks
rm(baseline, year1, year2)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Off the map. Not good.

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)+
		ylab("Mean Cohen's H")

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- DITCH THE ZEROES	--------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d", "functional.wd")])

# Figure out which sites are losers
sites.ls <- functional.wd |>
	mutate(Total = Year2.proportion+Year1.proportion+Baseline.proportion+
				 	cohen.H.year1+cohen.H.year2) |>
	dplyr::select(Plot, Feeder, Functional, Total) |>
	filter(Total == 0)

ruderal.ls <- sites.ls |> filter(Functional == 'Ruderal')
ruderal.ls <- unique(ruderal.ls$Plot)

browsed.ls <- sites.ls |> filter(Functional == 'Browsed perennials')
browsed.ls <- unique(browsed.ls$Plot)

resistant.ls <- sites.ls |> filter(Functional == 'Resistant perennials')
resistant.ls <- unique(resistant.ls$Plot)

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter < 6) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter < 6) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Subset the data based on sites with no changes
browsed <- functional |>
	filter(!Plot %in% browsed.ls & Functional == 'Browsed perennials')
ruderal <- functional |>
	filter(!Plot %in% ruderal.ls & Functional == 'Ruderal')
resistant <- functional |>
	filter(!Plot %in% resistant.ls & Functional == 'Resistant perennials')

# Subset by year browsed
browsed.baseline <- browsed |>
	filter(Time == '0')
browsed.baseline <- browsed.baseline[, -c(2,5:6)]
colnames(browsed.baseline)[4] <- 'Baseline.proportion'

browsed.year1 <- browsed |>
	filter(Time == '1')
browsed.year1 <- browsed.year1[, -c(2,5:6)]
colnames(browsed.year1)[4] <- 'Year1.proportion'

browsed.year2 <- browsed |>
	filter(Time == '2')
browsed.year2 <- browsed.year2[, -c(2,5:6)]
colnames(browsed.year2)[4] <- 'Year2.proportion'

# Merge
browsed.wd <- merge(browsed.year2, browsed.year1, all.x = TRUE)
browsed.wd <- merge(browsed.wd, browsed.baseline, all.x = TRUE)

# Subset by year ruderal
ruderal.baseline <- ruderal |>
	filter(Time == '0')
ruderal.baseline <- ruderal.baseline[, -c(2,5:6)]
colnames(ruderal.baseline)[4] <- 'Baseline.proportion'

ruderal.year1 <- ruderal |>
	filter(Time == '1')
ruderal.year1 <- ruderal.year1[, -c(2,5:6)]
colnames(ruderal.year1)[4] <- 'Year1.proportion'

ruderal.year2 <- ruderal |>
	filter(Time == '2')
ruderal.year2 <- ruderal.year2[, -c(2,5:6)]
colnames(ruderal.year2)[4] <- 'Year2.proportion'

# Merge
ruderal.wd <- merge(ruderal.year2, ruderal.year1, all.x = TRUE)
ruderal.wd <- merge(ruderal.wd, ruderal.baseline, all.x = TRUE)

# Subset by year resistant
resistant.baseline <- resistant |>
	filter(Time == '0')
resistant.baseline <- resistant.baseline[, -c(2,5:6)]
colnames(resistant.baseline)[4] <- 'Baseline.proportion'

resistant.year1 <- resistant |>
	filter(Time == '1')
resistant.year1 <- resistant.year1[, -c(2,5:6)]
colnames(resistant.year1)[4] <- 'Year1.proportion'

resistant.year2 <- resistant |>
	filter(Time == '2')
resistant.year2 <- resistant.year2[, -c(2,5:6)]
colnames(resistant.year2)[4] <- 'Year2.proportion'

# Merge
resistant.wd <- merge(resistant.year2, resistant.year1, all.x = TRUE)
resistant.wd <- merge(resistant.wd, resistant.baseline, all.x = TRUE)

# Combine
functional.wd <- rbind(browsed.wd, resistant.wd) # No possible ruderal sites

# Clear the decks
rm(browsed.baseline, browsed.year1, browsed.year2, browsed.wd,
	 ruderal.baseline, ruderal.year1, ruderal.year2, ruderal.wd,
	 resistant.baseline, resistant.year1, resistant.year2, resistant.wd,
	 browsed, resistant, ruderal)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Off the map. Not good.

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)+
		ylab("Mean Cohen's H")

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- FUNCTIONAL 6-10 m ------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d")])

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter %in% 6:10) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter %in% 6:10) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Subset by year
baseline <- functional |>
	filter(Time == '0')
baseline <- baseline[, -c(2,5:6)]
colnames(baseline)[4] <- 'Baseline.proportion'

year1 <- functional |>
	filter(Time == '1')
year1 <- year1[, -c(2,5:6)]
colnames(year1)[4] <- 'Year1.proportion'

year2 <- functional |>
	filter(Time == '2')
year2 <- year2[, -c(2,5:6)]
colnames(year2)[4] <- 'Year2.proportion'

# Merge
functional.wd <- merge(year2, year1, all.x = TRUE)
functional.wd <- merge(functional.wd, baseline, all.x = TRUE)

# Clear the decks
rm(baseline, year1, year2)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Off the map. Not good.

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)+
		ylab("Mean Cohen's H")

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- DITCH THE ZEROES --------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d", "functional.wd")])

# Figure out which sites are losers
sites.ls <- functional.wd |>
	mutate(Total = Year2.proportion+Year1.proportion+Baseline.proportion+
				 	cohen.H.year1+cohen.H.year2) |>
	dplyr::select(Plot, Feeder, Functional, Total) |>
	filter(Total == 0)

ruderal.ls <- sites.ls |> filter(Functional == 'Ruderal')
ruderal.ls <- unique(ruderal.ls$Plot)

browsed.ls <- sites.ls |> filter(Functional == 'Browsed perennials')
browsed.ls <- unique(browsed.ls$Plot)

resistant.ls <- sites.ls |> filter(Functional == 'Resistant perennials')
resistant.ls <- unique(resistant.ls$Plot)

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter %in% 6:10) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter %in% 6:10) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Subset the data based on sites with no changes
browsed <- functional |>
	filter(!Plot %in% browsed.ls & Functional == 'Browsed perennials')
ruderal <- functional |>
	filter(!Plot %in% ruderal.ls & Functional == 'Ruderal')
resistant <- functional |>
	filter(!Plot %in% resistant.ls & Functional == 'Resistant perennials')

# Subset by year browsed
browsed.baseline <- browsed |>
	filter(Time == '0')
browsed.baseline <- browsed.baseline[, -c(2,5:6)]
colnames(browsed.baseline)[4] <- 'Baseline.proportion'

browsed.year1 <- browsed |>
	filter(Time == '1')
browsed.year1 <- browsed.year1[, -c(2,5:6)]
colnames(browsed.year1)[4] <- 'Year1.proportion'

browsed.year2 <- browsed |>
	filter(Time == '2')
browsed.year2 <- browsed.year2[, -c(2,5:6)]
colnames(browsed.year2)[4] <- 'Year2.proportion'

# Merge
browsed.wd <- merge(browsed.year2, browsed.year1, all.x = TRUE)
browsed.wd <- merge(browsed.wd, browsed.baseline, all.x = TRUE)

# Subset by year ruderal
ruderal.baseline <- ruderal |>
	filter(Time == '0')
ruderal.baseline <- ruderal.baseline[, -c(2,5:6)]
colnames(ruderal.baseline)[4] <- 'Baseline.proportion'

ruderal.year1 <- ruderal |>
	filter(Time == '1')
ruderal.year1 <- ruderal.year1[, -c(2,5:6)]
colnames(ruderal.year1)[4] <- 'Year1.proportion'

ruderal.year2 <- ruderal |>
	filter(Time == '2')
ruderal.year2 <- ruderal.year2[, -c(2,5:6)]
colnames(ruderal.year2)[4] <- 'Year2.proportion'

# Merge
ruderal.wd <- merge(ruderal.year2, ruderal.year1, all.x = TRUE)
ruderal.wd <- merge(ruderal.wd, ruderal.baseline, all.x = TRUE)

# Subset by year resistant
resistant.baseline <- resistant |>
	filter(Time == '0')
resistant.baseline <- resistant.baseline[, -c(2,5:6)]
colnames(resistant.baseline)[4] <- 'Baseline.proportion'

resistant.year1 <- resistant |>
	filter(Time == '1')
resistant.year1 <- resistant.year1[, -c(2,5:6)]
colnames(resistant.year1)[4] <- 'Year1.proportion'

resistant.year2 <- resistant |>
	filter(Time == '2')
resistant.year2 <- resistant.year2[, -c(2,5:6)]
colnames(resistant.year2)[4] <- 'Year2.proportion'

# Merge
resistant.wd <- merge(resistant.year2, resistant.year1, all.x = TRUE)
resistant.wd <- merge(resistant.wd, resistant.baseline, all.x = TRUE)

# Combine
functional.wd <- rbind(browsed.wd, ruderal.wd, resistant.wd)

# Clear the decks
rm(browsed.baseline, browsed.year1, browsed.year2, browsed.wd,
	 ruderal.baseline, ruderal.year1, ruderal.year2, ruderal.wd,
	 resistant.baseline, resistant.year1, resistant.year2, resistant.wd,
	 browsed, resistant, ruderal)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Normal!

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)+
		ylab("Mean Cohen's H")

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- FUNCTIONAL 11-15 m ------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d")])

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter %in% 11:15) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter %in% 11:15) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Subset by year
baseline <- functional |>
	filter(Time == '0')
baseline <- baseline[, -c(2,5:6)]
colnames(baseline)[4] <- 'Baseline.proportion'

year1 <- functional |>
	filter(Time == '1')
year1 <- year1[, -c(2,5:6)]
colnames(year1)[4] <- 'Year1.proportion'

year2 <- functional |>
	filter(Time == '2')
year2 <- year2[, -c(2,5:6)]
colnames(year2)[4] <- 'Year2.proportion'

# Merge
functional.wd <- merge(year2, year1, all.x = TRUE)
functional.wd <- merge(functional.wd, baseline, all.x = TRUE)

# Clear the decks
rm(baseline, year1, year2)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Off the map. Not good.

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)+
		ylab("Mean Cohen's H")

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- DITCH THE ZEROES --------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d", "functional.wd")])

# Figure out which sites are losers
sites.ls <- functional.wd |>
	mutate(Total = Year2.proportion+Year1.proportion+Baseline.proportion+
				 	cohen.H.year1+cohen.H.year2) |>
	dplyr::select(Plot, Feeder, Functional, Total) |>
	filter(Total == 0)

ruderal.ls <- sites.ls |> filter(Functional == 'Ruderal')
ruderal.ls <- unique(ruderal.ls$Plot)

browsed.ls <- sites.ls |> filter(Functional == 'Browsed perennials')
browsed.ls <- unique(browsed.ls$Plot)

resistant.ls <- sites.ls |> filter(Functional == 'Resistant perennials')
resistant.ls <- unique(resistant.ls$Plot)

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter %in% 11:15) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter %in% 11:15) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Subset the data based on sites with no changes
browsed <- functional |>
	filter(!Plot %in% browsed.ls & Functional == 'Browsed perennials')
ruderal <- functional |>
	filter(!Plot %in% ruderal.ls & Functional == 'Ruderal')
resistant <- functional |>
	filter(!Plot %in% resistant.ls & Functional == 'Resistant perennials')

# Subset by year browsed
browsed.baseline <- browsed |>
	filter(Time == '0')
browsed.baseline <- browsed.baseline[, -c(2,5:6)]
colnames(browsed.baseline)[4] <- 'Baseline.proportion'

browsed.year1 <- browsed |>
	filter(Time == '1')
browsed.year1 <- browsed.year1[, -c(2,5:6)]
colnames(browsed.year1)[4] <- 'Year1.proportion'

browsed.year2 <- browsed |>
	filter(Time == '2')
browsed.year2 <- browsed.year2[, -c(2,5:6)]
colnames(browsed.year2)[4] <- 'Year2.proportion'

# Merge
browsed.wd <- merge(browsed.year2, browsed.year1, all.x = TRUE)
browsed.wd <- merge(browsed.wd, browsed.baseline, all.x = TRUE)

# Subset by year ruderal
ruderal.baseline <- ruderal |>
	filter(Time == '0')
ruderal.baseline <- ruderal.baseline[, -c(2,5:6)]
colnames(ruderal.baseline)[4] <- 'Baseline.proportion'

ruderal.year1 <- ruderal |>
	filter(Time == '1')
ruderal.year1 <- ruderal.year1[, -c(2,5:6)]
colnames(ruderal.year1)[4] <- 'Year1.proportion'

ruderal.year2 <- ruderal |>
	filter(Time == '2')
ruderal.year2 <- ruderal.year2[, -c(2,5:6)]
colnames(ruderal.year2)[4] <- 'Year2.proportion'

# Merge
ruderal.wd <- merge(ruderal.year2, ruderal.year1, all.x = TRUE)
ruderal.wd <- merge(ruderal.wd, ruderal.baseline, all.x = TRUE)

# Subset by year resistant
resistant.baseline <- resistant |>
	filter(Time == '0')
resistant.baseline <- resistant.baseline[, -c(2,5:6)]
colnames(resistant.baseline)[4] <- 'Baseline.proportion'

resistant.year1 <- resistant |>
	filter(Time == '1')
resistant.year1 <- resistant.year1[, -c(2,5:6)]
colnames(resistant.year1)[4] <- 'Year1.proportion'

resistant.year2 <- resistant |>
	filter(Time == '2')
resistant.year2 <- resistant.year2[, -c(2,5:6)]
colnames(resistant.year2)[4] <- 'Year2.proportion'

# Merge
resistant.wd <- merge(resistant.year2, resistant.year1, all.x = TRUE)
resistant.wd <- merge(resistant.wd, resistant.baseline, all.x = TRUE)

# Combine
functional.wd <- rbind(browsed.wd, ruderal.wd, resistant.wd)

# Clear the decks
rm(browsed.baseline, browsed.year1, browsed.year2, browsed.wd,
	 ruderal.baseline, ruderal.year1, ruderal.year2, ruderal.wd,
	 resistant.baseline, resistant.year1, resistant.year2, resistant.wd,
	 browsed, resistant, ruderal)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Normal!

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)+
		ylab("Mean Cohen's H")

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- FUNCTIONAL 0-10 m ------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d")])

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter < 11) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter < 11) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Subset by year
baseline <- functional |>
	filter(Time == '0')
baseline <- baseline[, -c(2,5:6)]
colnames(baseline)[4] <- 'Baseline.proportion'

year1 <- functional |>
	filter(Time == '1')
year1 <- year1[, -c(2,5:6)]
colnames(year1)[4] <- 'Year1.proportion'

year2 <- functional |>
	filter(Time == '2')
year2 <- year2[, -c(2,5:6)]
colnames(year2)[4] <- 'Year2.proportion'

# Merge
functional.wd <- merge(year2, year1, all.x = TRUE)
functional.wd <- merge(functional.wd, baseline, all.x = TRUE)

# Clear the decks
rm(baseline, year1, year2)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Off the map. Not good.

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)+
		ylab("Mean Cohen's H")

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- DITCH THE ZEROES --------------------------------------------

# Clear the decks
rm(list=ls()[! ls() %in% c("d", "functional.wd")])

# Figure out which sites are losers
sites.ls <- functional.wd |>
	mutate(Total = Year2.proportion+Year1.proportion+Baseline.proportion+
				 	cohen.H.year1+cohen.H.year2) |>
	dplyr::select(Plot, Feeder, Functional, Total) |>
	filter(Total == 0)

ruderal.ls <- sites.ls |> filter(Functional == 'Ruderal')
ruderal.ls <- unique(ruderal.ls$Plot)

browsed.ls <- sites.ls |> filter(Functional == 'Browsed perennials')
browsed.ls <- unique(browsed.ls$Plot)

resistant.ls <- sites.ls |> filter(Functional == 'Resistant perennials')
resistant.ls <- unique(resistant.ls$Plot)

# Tally natives and invasives at each plot*time
sites <- d |>
	group_by(Plot, Time, Feeder) |>
	summarize(Everything = n())

total <- d |>
	filter(Meter < 11) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- d |>
	filter(Meter < 11) |>
#	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = Functional, values_from = Count)

functional <- functional |>
	pivot_longer(cols = 4:7, names_to = 'Functional', values_to = 'Count')

# Calculate proportion
functional <- merge(functional, total, all = TRUE)

# Drop the NAs in functional group
functional <- functional |> filter(Functional != 'NA')

# Address zero count
functional[is.na(functional)] <- 0 # Address zero counts

# Calculate proportion
functional <- functional |>
	mutate(Proportion = Count/Total)
functional[is.na(functional)] <- 0 # Address NANs created by division

# Reorganize
functional <- functional |>
	dplyr::select(Plot, Time, Feeder, Functional, Count, Total, Proportion)

# Clear the decks
rm(total, sites)

# Subset the data based on sites with no changes
browsed <- functional |>
	filter(!Plot %in% browsed.ls & Functional == 'Browsed perennials')
ruderal <- functional |>
	filter(!Plot %in% ruderal.ls & Functional == 'Ruderal')
resistant <- functional |>
	filter(!Plot %in% resistant.ls & Functional == 'Resistant perennials')

# Subset by year browsed
browsed.baseline <- browsed |>
	filter(Time == '0')
browsed.baseline <- browsed.baseline[, -c(2,5:6)]
colnames(browsed.baseline)[4] <- 'Baseline.proportion'

browsed.year1 <- browsed |>
	filter(Time == '1')
browsed.year1 <- browsed.year1[, -c(2,5:6)]
colnames(browsed.year1)[4] <- 'Year1.proportion'

browsed.year2 <- browsed |>
	filter(Time == '2')
browsed.year2 <- browsed.year2[, -c(2,5:6)]
colnames(browsed.year2)[4] <- 'Year2.proportion'

# Merge
browsed.wd <- merge(browsed.year2, browsed.year1, all.x = TRUE)
browsed.wd <- merge(browsed.wd, browsed.baseline, all.x = TRUE)

# Subset by year ruderal
ruderal.baseline <- ruderal |>
	filter(Time == '0')
ruderal.baseline <- ruderal.baseline[, -c(2,5:6)]
colnames(ruderal.baseline)[4] <- 'Baseline.proportion'

ruderal.year1 <- ruderal |>
	filter(Time == '1')
ruderal.year1 <- ruderal.year1[, -c(2,5:6)]
colnames(ruderal.year1)[4] <- 'Year1.proportion'

ruderal.year2 <- ruderal |>
	filter(Time == '2')
ruderal.year2 <- ruderal.year2[, -c(2,5:6)]
colnames(ruderal.year2)[4] <- 'Year2.proportion'

# Merge
ruderal.wd <- merge(ruderal.year2, ruderal.year1, all.x = TRUE)
ruderal.wd <- merge(ruderal.wd, ruderal.baseline, all.x = TRUE)

# Subset by year resistant
resistant.baseline <- resistant |>
	filter(Time == '0')
resistant.baseline <- resistant.baseline[, -c(2,5:6)]
colnames(resistant.baseline)[4] <- 'Baseline.proportion'

resistant.year1 <- resistant |>
	filter(Time == '1')
resistant.year1 <- resistant.year1[, -c(2,5:6)]
colnames(resistant.year1)[4] <- 'Year1.proportion'

resistant.year2 <- resistant |>
	filter(Time == '2')
resistant.year2 <- resistant.year2[, -c(2,5:6)]
colnames(resistant.year2)[4] <- 'Year2.proportion'

# Merge
resistant.wd <- merge(resistant.year2, resistant.year1, all.x = TRUE)
resistant.wd <- merge(resistant.wd, resistant.baseline, all.x = TRUE)

# Combine
functional.wd <- rbind(browsed.wd, ruderal.wd, resistant.wd)

# Clear the decks
rm(browsed.baseline, browsed.year1, browsed.year2, browsed.wd,
	 ruderal.baseline, ruderal.year1, ruderal.year2, ruderal.wd,
	 resistant.baseline, resistant.year1, resistant.year2, resistant.wd,
	 browsed, resistant, ruderal)

# Calculate the delta from baseline to subsequent years
functional.wd <- functional.wd |>
	mutate(cohen.H.year1 = ES.h(Year1.proportion, Baseline.proportion),
				 cohen.H.year2 = ES.h(Year2.proportion, Baseline.proportion))

functional <- functional.wd |>
	dplyr::select(Plot, Feeder, Functional, cohen.H.year1, cohen.H.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "cohen.H")

descdist(functional$cohen.H, discrete = FALSE) # Normal!

m1 <- lmer(cohen.H~Feeder * Year * Functional + (1|Plot), data = functional)
Anova(m1, type = 3)

cohen.H.df <- functional |>
	group_by(Year, Functional, Feeder) |>
	summarize(mean = mean(cohen.H),
						n = n(),
						sd = sd(cohen.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(cohen.H.df, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = LCL, ymax = UCL),
        width = 0.1,
        position=position_dodge(width=0.5))+
		facet_wrap(~Functional)+
		ylab("Mean Cohen's H")

control <- functional.wd |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.H.year1'
colnames(control)[5] <- 'Control.H.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- functional.wd |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.H.year1'
colnames(feeder)[5] <- 'Feeder.H.year2'
feeder <- feeder[, -2]

Delta <- cbind(feeder, control)

Delta <- Delta |>
	mutate(Delta.1 = Feeder.H.year1-Control.H.year1,
				 Delta.2 = Feeder.H.year2-Control.H.year2)

Delta <- Delta |>
	dplyr::select(Plot, Functional, Delta.1, Delta.2) |>
	pivot_longer(3:4, names_to = "Year", values_to = "Delta.H")

Delta <- Delta |>
	group_by(Year, Functional) |>
	summarize(mean = mean(Delta.H),
						n = n(),
						sd = sd(Delta.H),
						se = sd/sqrt(n),
						margin = qt(0.975,df=n-1)*sd/sqrt(n),
						LCL = mean-margin,
						UCL = mean+margin)

ggplot(Delta, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5))+
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))+
		ylab("Mean Cohen's H Delta")+
		xlab("")+
		theme_bw()+
		facet_wrap(~Functional)

## --------------- OLD DELTA ---------------------------------------------------
# Calculate the delta from baseline to subsequent years
manip <- manip |>
	mutate(Delta.year1 = Year1.proportion-Baseline.proportion,
				 Delta.year2 = Year2.proportion-Baseline.proportion)

delta.sites <- manip |>
	dplyr::select(Plot, Feeder, Invasive, Delta.year1, Delta.year2) |>
	pivot_longer(4:5, names_to = "Year", values_to = "Delta") |>
	filter(Invasive == 'Native')

descdist(delta$Delta, discrete = FALSE)

delta <- delta |>
	group_by(Feeder, Year) |>
	summarize(mean = mean(Delta),
						n = n(),
						sd = sd(Delta),
						se = sd/sqrt(n))

ggplot(delta, aes(x = Year, y = mean, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        linetype = "dotted",
        position=position_dodge(width=0.5))

ggplot(delta.sites, aes(x = Year, y = Delta, color = Feeder))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
		facet_wrap(~Plot)

control <- manip |>
	filter(Feeder == 'N')
control <- control[, c(1,2,3,7,8)]
colnames(control)[4] <- 'Control.delta.year1'
colnames(control)[5] <- 'Control.delta.year2'
control <- control[, -2]
control <- control[, 3:4]

feeder <- manip |>
	filter(Feeder == 'Y')
feeder <- feeder[, c(1,2,3,7,8)]
colnames(feeder)[4] <- 'Feeder.delta.year1'
colnames(feeder)[5] <- 'Feeder.delta.year2'
feeder <- feeder[, -2]

LRR <- cbind(feeder, control)

LRR <- LRR |>
	filter(Invasive == "Invasive")

LRR <- LRR |>
	mutate(LRR.1 = Feeder.delta.year1-Control.delta.year1,
				 LRR.2 = Feeder.delta.year2-Control.delta.year2)

LRR <- LRR |>
	dplyr::select(Plot, LRR.1, LRR.2) |>
	pivot_longer(2:3, names_to = "Year", values_to = "LRR")

LRR <- LRR |>
	group_by(Year) |>
	summarize(mean = mean(LRR),
						n = n(),
						sd = sd(LRR),
						se = sd/sqrt(n))

ggplot(LRR, aes(x = Year, y = mean))+
    geom_point(size = 4, position=position_dodge(width=0.5)) +
    geom_errorbar(
        aes(ymin = mean-se, ymax = mean+se),
        width = 0.1,
        position=position_dodge(width=0.5))