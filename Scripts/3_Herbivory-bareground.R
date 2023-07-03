## --------------- HEADER ------------------------------------------------------
## Script name: 3_Herbivory-bareground.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-22
## Date Last modified: 2023-06-22
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
	filter(Meter < 26)

# Bring in the data
retro.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Retrospective') |>
	filter(Meter < 26)

# ditch.ls <- c('P2', 'P3', 'P5', 'P6', 'P7', 'P8', 'WS4', 'WS8')

## --------------- PREPARE DATA ------------------------------------------------

# Bin the distance values
feeder.veg$Bin <- NA
for(i in 1:nrow(feeder.veg)){
	if(feeder.veg$Meter[i] %in% 1:5){
		feeder.veg$Bin[i] <- "1-5"
	}
	if(feeder.veg$Meter[i] %in% 6:10){
		feeder.veg$Bin[i] <- "6-10"
	}
	if(feeder.veg$Meter[i] %in% 11:15){
		feeder.veg$Bin[i] <- "11-15"
	}
	if(feeder.veg$Meter[i] %in% 16:20){
		feeder.veg$Bin[i] <- "16-20"
	}
	if(feeder.veg$Meter[i] %in% 21:25){
		feeder.veg$Bin[i] <- "21-25"
	}
}

# Assign ordinality to bin value
feeder.veg$Bin <- factor(feeder.veg$Bin, order = TRUE,
													levels = c("1-5", "6-10", "11-15",
																		 "16-20", "21-25"))

# Drop non browse and tally
herbivory <- feeder.veg |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Bin, Browse) |>
	filter(Browse > 0) |>
	group_by(Site, Plot, Time, Feeder, Transect, Bin) |>
	summarize(Browsed = n())

ggplot(herbivory, aes(x = Bin, y = Browsed))+
	geom_jitter()+
	theme_bw()+
	facet_wrap(~Feeder)

all <- feeder.veg |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Bin, Browse) |>
	group_by(Site, Plot, Time, Feeder, Transect, Bin) |>
	summarize(Everything = n())

herbivory.all <- merge(herbivory, all, all = TRUE)
rm(all)

herbivory.all <- herbivory.all |>
	dplyr::select(-Everything)

# Model with zeroes

herbivory.all[is.na(herbivory.all)] <- 0

# Check the rows
v1 <- length(which(herbivory.all$Time==0))
v2 <- v1/2 # 2 treatments
v3 <- v2/4 # 4 transects
v4 <- v3/5 # 5 bins
print(v4)
# There are missing transects, so this is likely correct

v1 <- length(which(herbivory.all$Time==1))
v2 <- v1/2 # 2 treatments
v3 <- v2/4 # 4 transects
v4 <- v3/5 # 5 bins
print(v4)
# There are missing transects, so this is likely correct

v1 <- length(which(herbivory.all$Time==2))
v2 <- v1/2 # 2 treatments
v3 <- v2/8 # 8 transects
v4 <- v3/5 # 5 bins
print(v4)
# There are missing transects, so this is likely correct

## --------------- MANIPULATIVE BROWSE NEG BINOM ALL ZEROES --------------------

descdist(herbivory$Browsed, discrete = TRUE) # negative binomial

# Baseline has zero fitted value. Let's drop it
herbivory.all <- herbivory.all |> filter(Time != 0)
herbivory.all$Time <- as_factor(herbivory.all$Time)

null <- glm.nb(Browsed ~ 1, data = herbivory.all)
m1 <- glm.nb(Browsed ~ Feeder * Bin + Time + Site/Plot,
						data = herbivory.all)
m2 <- glmer.nb(Browsed ~ Feeder * Bin + Time + (1|Plot), # Fail
						data = herbivory.all)
m3 <- glmer.nb(Browsed ~ Feeder * Bin + Time + (1|Site),
						data = herbivory.all)
m4 <- glmer.nb(Browsed ~ Feeder * Bin + Time + (1|Site/Plot), # Fail
						data = herbivory.all)
m5 <- glmer.nb(Browsed ~ Feeder * Bin * Time + (1|Site/Plot),  # Fail
						data = herbivory.all)
m6 <- glm.nb(Browsed ~ Feeder * Bin + Time + Plot + Site,
						data = herbivory.all)
m7 <- glm.nb(Browsed ~ Feeder * Bin + Time + Plot,
						data = herbivory.all)
m8 <- glm.nb(Browsed ~ Feeder * Bin + Time + Site,
						data = herbivory.all)
m9 <- glm.nb(Browsed ~ Feeder * Bin * Time + Site,
						data = herbivory.all)
AIC(null, m1, m2, m3, m4, m5, m6, m7, m8, m9) # m1 wins
anova(null, m1) # Not working

# Clear the decks
rm(list = ls()[!ls() %in% c("herbivory.all", "feeder.veg", "m1")])

# Check model
Anova(m1)
emmeans(m1, pairwise~Feeder*Bin, type = 'response')
emmeans(m1, pairwise~Time, type = 'response')
check_zeroinflation(m1) # Pass
sim.m1 <- simulateResiduals(m1)
plot(sim.m1) # Pass
testDispersion(sim.m1) # Pass
testZeroInflation(sim.m1) # Pass
testQuantiles(sim.m1) # Pass
plotResiduals(sim.m1, herbivory.all$Feeder) # Pass
plotResiduals(sim.m1, herbivory.all$Bin) # Pass
plotResiduals(sim.m1, herbivory.all$Site) # Pass 
plotResiduals(sim.m1, herbivory.all$Plot) # Pass
plotResiduals(sim.m1, herbivory.all$Time) # Pass

# Export outcome
herb.means <- tidy(emmeans(m1, ~Feeder*Bin, type = 'response'))

herb.means$Bin <- factor(herb.means$Bin, order = TRUE,
													levels = c("1-5", "6-10", "11-15",
																		 "16-20", "21-25"))

## --------------- MANIPULATIVE VISUALIZE BROWSE DATA AND MODEL ----------------

# Model means
ggplot(herb.means, aes(x = Bin, y = response, fill = Feeder))+
	geom_point(shape = 21, size = 5,
						 position = position_dodge(width = 0.5))+
	# geom_linerange(aes(ymin = response-std.error, ymax = response+std.error),
								 # position = position_dodge(width = 0.5))+
	scale_fill_manual(values = c('gray', 'skyblue'),
										labels = c('Control', 'Feeder'),
										guide = guide_legend(reverse = TRUE),
										name = element_blank())+
	theme_bw()



# All data
for(i in 1:nrow(herbivory.all)){
	if(is.na(herbivory.all$Browsed[i])){
		herbivory.all$Browsed[i] <- 0
	}
}

herbivory.all <- herbivory.all |> filter(Time != '0')
ggplot(herbivory.all, aes(x = Bin, y = Browsed, fill = Feeder))+
	geom_jitter(shape = 21, size = 5, alpha = 0.75, stroke = 1)+
	scale_fill_manual(values = c('gray', 'skyblue'),
										labels = c('Control', 'Feeder'),
										guide = guide_legend(reverse = TRUE),
										name = element_blank())+
	scale_y_continuous(limits = c(-0.43,5))+
	ylab('Browsed plants at transect')+
	xlab('Distance from Feeder')+
	theme_bw()+
	theme(panel.grid.minor.y = element_blank())+
	facet_wrap(~Time)

## --------------- MANIPULATIVE OTHER METRICS OF HERBIVORY ---------------------

# Browse percentage
rm(list = ls()[!ls() %in% c("feeder.veg", "retro.veg")])

feeder.veg$Browse.percentage <- as.numeric(feeder.veg$Browse.percentage)

browse.percent <- feeder.veg |>
	filter(Browse == '1') |>
	filter(!is.na(Browse.percentage))

ggplot(browse.percent, aes(x = Meter, y = Browse.percentage, fill = Feeder))+
	geom_point(shape = 21)+
	facet_wrap(~Feeder)

ggplot(browse.percent, aes(x = Bin, y = Browse.percentage, fill = Feeder))+
	geom_jitter(shape = 21, size = 4, stroke = 1)+
	facet_wrap(~Feeder)+
	theme_bw()

descdist(browse.percent$Browse.percentage, discrete = FALSE) # Beta

library(betareg)

n <- nrow(browse.percent)
browse.percent <- browse.percent |>
	mutate(Browse.percentage.trans = ((Browse.percentage * (n-1) + 0.5)/n)/100)

m1 <- betareg(Browse.percentage.trans ~ Bin*Feeder, data = browse.percent)
Anova(m1)

m1 <- lm(Browse.percentage ~ Bin*Feeder, data = browse.percent)
Anova(m1)

# Bite count
rm(list = ls()[!ls() %in% c("feeder.veg", "retro.veg")])

feeder.veg$Bite.count <- as.numeric(feeder.veg$Bite.count)

bite.count <- feeder.veg |>
	filter(Bite.count > 0)

ggplot(bite.count, aes(x = Meter, y = Bite.count, fill = Feeder))+
	geom_point(shape = 21)+
	facet_wrap(~Feeder)

ggplot(bite.count, aes(x = Bin, y = Bite.count, fill = Feeder))+
	geom_jitter(shape = 21, size = 4, stroke = 1)+
	facet_wrap(~Feeder)+
	theme_bw()

descdist(bite.count$Bite.count, discrete = TRUE) # negative binomial

m1 <- glm.nb(Bite.count ~ Bin*Feeder, data = bite.count)
Anova(m1, type = 3)
check_model(m1)

## --------------- RETROSPECTIVE HERBIVORY -------------------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("retro.veg")])

# Bin the distance values
retro.veg$Bin <- NA
for(i in 1:nrow(retro.veg)){
	if(retro.veg$Meter[i] %in% 1:5){
		retro.veg$Bin[i] <- "1-5"
	}
	if(retro.veg$Meter[i] %in% 6:10){
		retro.veg$Bin[i] <- "6-10"
	}
	if(retro.veg$Meter[i] %in% 11:15){
		retro.veg$Bin[i] <- "11-15"
	}
	if(retro.veg$Meter[i] %in% 16:20){
		retro.veg$Bin[i] <- "16-20"
	}
	if(retro.veg$Meter[i] %in% 21:25){
		retro.veg$Bin[i] <- "21-25"
	}
}

retro.veg$Bin <- factor(retro.veg$Bin, order = TRUE,
													levels = c("1-5", "6-10", "11-15",
																		 "16-20", "21-25"))

retro.sum <- retro.veg |>
	group_by(Plot, Feeder, Bin) |>
	summarize(Browsed = sum(Browse))

ggplot(retro.sum, aes(x = Bin, y = Browsed, fill = Feeder))+
	geom_jitter(shape = 21, size = 5, alpha = 0.75, stroke = 1)+
	scale_fill_manual(values = c('gray', 'skyblue'),
										labels = c('Control', 'Feeder'),
										guide = guide_legend(reverse = TRUE),
										name = element_blank())+
	scale_y_continuous(limits = c(-0.43,5))+
	ylab('Browsed plants at transect')+
	xlab('Distance from Feeder')+
	theme_bw()+
	theme(panel.grid.minor.y = element_blank())+
	facet_wrap(~Feeder)

descdist(retro.sum$Browsed, discrete = TRUE)

m1 <- glmer.nb(Browsed ~ Feeder * Bin + (1|Plot),
						data = retro.sum)
Anova(m1)
check_model(m1)

## --------------- MANIPULATIVE BAREGROUND -------------------------------------

rm(list = ls()[!ls() %in% c("feeder.veg", "retro.veg")])

bare.ground <- feeder.veg |>
	filter(Time > 0) |>
	filter(Species == 'Bare Ground') |>
	group_by(Site, Plot, Time, Feeder, Transect, Bin) |>
	summarize(Successes = n()) |>
	mutate(Failures = 5-Successes)

m1 <- glmer(cbind(Successes,Failures) ~ Feeder * Bin + Time + (1|Plot), 
						data = bare.ground, family = 'binomial')

# Check model
Anova(m1)
check_model(m1)
emmeans(m1, pairwise~Feeder*Bin, type = 'response')
emmeans(m1, pairwise~Time, type = 'response')
sim.m1 <- simulateResiduals(m1)
plot(sim.m1) # Fail
testDispersion(sim.m1) # Fail
testZeroInflation(sim.m1) # Fail
testQuantiles(sim.m1) # Fail
plotResiduals(sim.m1, bare.ground$Feeder) # Fail
plotResiduals(sim.m1, bare.ground$Bin) # Fail
plotResiduals(sim.m1, bare.ground$Site) # Fail 
plotResiduals(sim.m1, bare.ground$Plot) # Fail
plotResiduals(sim.m1, bare.ground$Time) # Pass

# Export outcome
bg.means <- tidy(emmeans(m1, ~Feeder*Bin, type = 'response'))

bg.means$Bin <- factor(bg.means$Bin, order = TRUE,
													levels = c("1-5", "6-10", "11-15",
																		 "16-20", "21-25"))

bg.means$Feeder <- as_factor(bg.means$Feeder)
levels(bg.means$Feeder)[levels(bg.means$Feeder)=='N'] <- 'Control'
levels(bg.means$Feeder)[levels(bg.means$Feeder)=='Y'] <- 'Feeder'

# Visualize
ggplot(bg.means, aes(x = Bin, y = prob, fill = Feeder))+
	geom_linerange(aes(ymin = prob-std.error, ymax = prob+std.error),
								 linewidth = 1)+
	geom_point(shape = 21, stroke = 1.25, size = 5)+
	scale_fill_manual(values = c('gold', 'skyblue'))+
	scale_y_continuous(limits = c(0,0.55),
										 breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
	xlab('Distance')+
	ylab('Probability of detecting bare ground')+
	facet_wrap(~Feeder)+
	theme_bw()+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  panel.grid.minor = element_blank()
)

## --------------- RETROSPECTIVE BAREGROUND ------------------------------------

rm(list = ls()[!ls() %in% c("feeder.veg", "retro.veg")])

# Bin the distance values
retro.veg$Bin <- NA
for(i in 1:nrow(retro.veg)){
	if(retro.veg$Meter[i] %in% 1:5){
		retro.veg$Bin[i] <- "1-5"
	}
	if(retro.veg$Meter[i] %in% 6:10){
		retro.veg$Bin[i] <- "6-10"
	}
	if(retro.veg$Meter[i] %in% 11:15){
		retro.veg$Bin[i] <- "11-15"
	}
	if(retro.veg$Meter[i] %in% 16:20){
		retro.veg$Bin[i] <- "16-20"
	}
	if(retro.veg$Meter[i] %in% 21:25){
		retro.veg$Bin[i] <- "21-25"
	}
}

retro.veg$Bin <- factor(retro.veg$Bin, order = TRUE,
													levels = c("1-5", "6-10", "11-15",
																		 "16-20", "21-25"))

bare.ground <- retro.veg |>
	filter(Species == 'bare ground') |>
	group_by(Plot, Feeder, Transect, Bin) |>
	summarize(Successes = n()) |>
	mutate(Failures = 5-Successes)

m1 <- glmer(cbind(Successes,Failures) ~ Feeder * Bin + (1|Plot), 
						data = bare.ground, family = 'binomial')

# Check model
Anova(m1)
check_model(m1)
emmeans(m1, pairwise~Feeder*Bin, type = 'response')
emmeans(m1, pairwise~Time, type = 'response')
sim.m1 <- simulateResiduals(m1)
plot(sim.m1) # Fail
testDispersion(sim.m1) # Fail
testZeroInflation(sim.m1) # Fail
testQuantiles(sim.m1) # Pass
plotResiduals(sim.m1, bare.ground$Feeder) # Pass
plotResiduals(sim.m1, bare.ground$Bin) # Pass
plotResiduals(sim.m1, bare.ground$Site) # Pass
plotResiduals(sim.m1, bare.ground$Plot) # Pass

# Export outcome
bg.means <- tidy(emmeans(m1, ~Feeder*Bin, type = 'response'))

bg.means$Bin <- factor(bg.means$Bin, order = TRUE,
													levels = c("1-5", "6-10", "11-15",
																		 "16-20", "21-25"))

bg.means$Feeder <- as_factor(bg.means$Feeder)
levels(bg.means$Feeder)[levels(bg.means$Feeder)=='N'] <- 'Control'
levels(bg.means$Feeder)[levels(bg.means$Feeder)=='Y'] <- 'Feeder'

# Visualize
ggplot(bg.means, aes(x = Bin, y = prob, fill = Feeder))+
	geom_linerange(aes(ymin = prob-std.error, ymax = prob+std.error),
								 linewidth = 1)+
	geom_point(shape = 21, stroke = 1.25, size = 5)+
	scale_fill_manual(values = c('gold', 'skyblue'))+
	scale_y_continuous(limits = c(0,0.75),
										 breaks = c(0, 0.25, 0.50, 0.75))+
	xlab('Distance')+
	ylab('Probability of detecting bare ground')+
	facet_wrap(~Feeder)+
	theme_bw()+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  panel.grid.minor = element_blank()
)
