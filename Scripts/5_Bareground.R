## --------------- HEADER ------------------------------------------------------
## Script name: 5_Bareground.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-22
## Date Last modified: 2023-010-06
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: 

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
library(DHARMa) # simulate residuals

# Clear the decks
rm(list=ls())

# Bring in the data
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Meter < 4) |>
	filter(Time != 5) 

# Create bare ground variable
feeder.veg$Bare.ground <- NA

for(i in 1:nrow(feeder.veg)){
	if(feeder.veg$Species[i] == 'Bare Ground'){
		 feeder.veg$Bare.ground[i] <- 1
	} else {
		 feeder.veg$Bare.ground[i] <- 0
	}
}

# We want to do this by point, not by observation
feeder.veg <- feeder.veg |>
	dplyr::select(Plot, Direction, Feeder, Time, Meter, Bare.ground)

feeder.veg <- unique(feeder.veg) # Drop extra rows

# Check that there are no overlapping values
check <- feeder.veg |>
	group_by(Plot, Direction, Feeder, Time, Meter) |>
	summarize(Count = n())

# Some rows have bare ground and a plant. We will not consider these disturbed
# if there are still plants (it's not disturbed enough for full plant removal)
# Drop 810 and 811
feeder.veg <- feeder.veg |>
	subset(!(Plot == 'WS2' & Direction == 'E' & Feeder == 'N' & 
				 	Time == 1 & Meter == 1 & Bare.ground == 0))

feeder.veg <- feeder.veg |>
	subset(!(Plot == 'WS2' & Direction == 'E' & Feeder == 'N' & 
				 	Time == 1 & Meter == 2 & Bare.ground == 0))

# Check that there are no overlapping values
check <- feeder.veg |>
	group_by(Plot, Direction, Feeder, Time, Meter) |>
	summarize(Count = n())

feeder.veg$Direction <- toupper(feeder.veg$Direction)

class(feeder.veg$Meter)

class(feeder.veg$Feeder)
feeder.veg$Feeder <- as_factor(feeder.veg$Feeder)

class(feeder.veg$Plot)
feeder.veg$Plot <- as_factor(feeder.veg$Plot)

class(feeder.veg$Time)
feeder.veg$Time <- as_factor(feeder.veg$Time)

class(feeder.veg$Direction)
feeder.veg$Direction <- as_factor(feeder.veg$Direction)

mod <- glmer(Bare.ground ~ Feeder * as_factor(Time) + (1|Plot/Direction),
						 data = feeder.veg, family = 'binomial')
Anova(mod)
summary(mod)

library(performance)
check_model(mod)

library(DHARMa)
sim.mod <- simulateResiduals(mod)
plot(sim.mod, quantreg=T) # outlier test significant

testOutliers(sim.mod, type = 'bootstrap') # fail
plotResiduals(sim.mod, feeder.veg$Plot, quantreg = T) # fail/pass
plotResiduals(sim.mod, feeder.veg$Time, quantreg = T) # ok 
plotResiduals(sim.mod, feeder.veg$Feeder, quantreg = T) # ok
plotResiduals(sim.mod, feeder.veg$Direction, quantreg = T) # pass/fail

# Overdispersion
testDispersion(sim.mod) # ok

# Good model fit (100% within bounds)
library(arm)
binnedplot(predict(mod, type="response", re.form=NA), 
					 resid(mod, type="response"), main='Without random effects', nclass=20)

# Influential variables
library(car)
4/length(unique(feeder.veg$Plot)) # 0.16 cook's D threshold. LV2 influential
2/sqrt(length(unique(feeder.veg$Plot))) # 0.4 for params (Belsey, Kuh, Welsch) 
# P4, LV2, P2 for time*feeder and time
influential <- influence(mod, group="Plot")
car::infIndexPlot(influential) # 
rm(influential)

4/length(unique(feeder.veg$Direction)) # 0.5 cook's D threshold, East influential
2/sqrt(length(unique(feeder.veg$Direction))) # 0.7 for params (Belsey, Kuh, Welsch) 
influential <- influence(mod, group="Direction")
car::infIndexPlot(influential) # 
rm(influential)

# Adjust for variance
lme4::VarCorr(mod)
total.sd <- sqrt(0.62827^2 + 0.99244^2) # 1.17
total.sd <- sqrt(sum(as.data.frame(lme4::VarCorr(mod))$vcov)) # 1.17

# Check contrasts of significant values
confint(emmeans(mod, revpairwise~Feeder, type = 'response'), sigma = total.sd)
emmeans(mod, pairwise~Feeder, type = 'response', sigma = total.sd)
emmeans(mod, pairwise~Feeder*Time, type = 'response', sigma = total.sd)
# Got worse every year

# Create dataframes 
feeder.means <- tibble(Feeder = c('Feeder', 'No feeder'),
								Bare.ground = c(0.2736, 0.0144),
								LCL = c(0.18964, 0.00789),
								UCL = c(0.3775, 0.0263))
write.csv(feeder.means, 'Results/Bare-ground/feeder-means.csv', row.names = FALSE)

time.means <- tibble(Time = c('Year 0', 'Year 0', 'Year 1', 
															'Year 1', 'Year 2', 'Year 2'),
								Feeder = c('Feeder', 'No feeder','Feeder', 
													 'No feeder', 'Feeder', 'No feeder'),
								Bare.ground = c(0.01770, 0.00981, 0.58409, 
															0.01354, 0.67872, 0.02259),
								LCL = c(0.00793, 0.00362, 0.46389, 
												0.00563, 0.57527, 0.01236),
								UCL = c(0.0390, 0.0263, 0.6951, 
												0.0322, 0.7672, 0.0409))
write.csv(time.means, 'Results/Bare-ground/time-means.csv', row.names = FALSE)

