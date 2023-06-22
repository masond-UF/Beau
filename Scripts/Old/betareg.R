library(tidyverse)
library(lme4)
library(car)
library(fitdistrplus)
## --------------- SET—UP WORKSPACE --------------------------------------------
# Clear the decks
rm(list=ls())

# Bring in the data
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative')

## --------------- PREPARE CLOSE DATA ------------------------------------------

ditch <- c('P2', 'P3', 'P5', 'P6', 'P7', 'P8', 'WS4', 'WS8')

# Subset the data
close <- feeder.veg |>
	filter(Meter %in% 1:25) |> # 5
	filter(!Plot %in% ditch)

for(i in 1:nrow(close)){
	if(is.na(close$Functional[i])){
		close$Functional[i] <- 'Native'
	}
}

# Combine the ruderal and resistant?

# Tally natives and invasives at each plot*time
sites <- close |>
	group_by(Site, Plot, Time, Feeder) |>
	summarize(Everything = n())

# Tally invasive functional groups
func.invasive <- close |>
	group_by(Site, Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

func.invasive <- merge(sites, func.invasive, all = TRUE)
func.invasive[is.na(func.invasive)] <- 0
func.invasive <- func.invasive |> dplyr::select(-Everything)

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

# Calculate the proportion of invasive
# Calculate the proportion of invasive
func.invasive <- func.invasive |>
	mutate(Total = `Resistant perennials`+
				 	`Browsed perennials`+Ruderal,
				 `Resistant perennials` = `Resistant perennials`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total,
				 Ruderal = Ruderal/Total,
				 Total.inv = (`Resistant perennials`+
				 	`Browsed perennials`+Ruderal)/(Total+Native)) |>
	dplyr::select(-Native, -Total)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive |>
	pivot_longer(5:8, names_to = 'Functional', values_to = 'Proportion')

func.invasive <- func.invasive |>
	mutate(Proportion.trans = (Proportion * (794-1) + 0.5)/794)

total <- func.invasive |>
	filter(Functional == 'Total.inv')
total.mod <- betareg(Proportion.trans ~ Feeder*as_factor(Time) + Site,
											 data = total, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
total.mod <- lmer(logit(Proportion) ~ Feeder*as_factor(Time) + (1|Site),
											 data = total)
Anova(total.mod)
check_model(total.mod)
emmeans(total.mod, ~Feeder*Time)

browsed <- func.invasive |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- betareg(Proportion.trans ~ Feeder*as_factor(Time) + Site,
											 data = browsed, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
browsed.mod <- lmer(logit(Proportion) ~ Feeder*as_factor(Time) + (1|Site),
											 data = browsed)
Anova(browsed.mod)
check_model(browsed.mod)
emmeans(browsed.mod, ~Feeder*Time)

ruderals <- func.invasive |>
	filter(Functional == 'Ruderal')
ruderals.mod <- betareg(Proportion.trans ~ Feeder*as_factor(Time) + Site,
											 data = ruderals, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
ruderals.mod <- lmer(logit(Proportion) ~ Feeder*as_factor(Time) + (1|Site),
											 data = ruderals)
Anova(ruderals.mod)
check_model(ruderals.mod)
emmeans(ruderals.mod, ~Feeder*Time)

resistant <- func.invasive |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- betareg(Proportion.trans ~ Feeder*as_factor(Time) + Site,
											 data = resistant, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
resistant.mod <- lmer(logit(Proportion) ~ Feeder*as_factor(Time) + (1|Site),
											 data = resistant)
Anova(resistant.mod)
check_model(resistant.mod)
emmeans(resistant.mod, ~Feeder*Time)

# Try binomial
for(i in 1:nrow(func.invasive)){
	if(func.invasive$Proportion[i] > 0){
		func.invasive$Proportion[i] <- 1
	}
}

total <- func.invasive |>
	filter(Functional == 'Total.inv')
total.mod <- glmer(Proportion ~ Feeder*as_factor(Time) + (1|Site),
									 data = total, family = 'binomial')
Anova(total.mod)
check_model(total.mod)
emmeans(total.mod, ~Feeder*Time)

browsed <- func.invasive |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- glmer(Proportion ~ Feeder*as_factor(Time) + (1|Site),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
emmeans(browsed.mod, ~Feeder*Time)

ruderals <- func.invasive |>
	filter(Functional == 'Ruderal')
ruderals.mod <- glmer(Proportion ~ Feeder*as_factor(Time) + (1|Site),
									 data = ruderals, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
emmeans(ruderals.mod, ~Feeder*Time)

resistant <- func.invasive |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- glmer(Proportion ~ Feeder*as_factor(Time) + (1|Site),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
check_model(resistant.mod)
emmeans(resistant.mod, ~Feeder*Time)

## --------------- BINOMIAL CLOSE ----------------------------------------------

# Try it successes and failures
# Clear the decks
rm(list=ls()[! ls() %in% c("close")])

# Count total number of invasive rows
# Get 1s and 0s for each transect point
# Tally natives and invasives at each plot*time

# Tally natives and invasives at each plot*time
sites <- close |>
	group_by(Site, Plot, Time, Feeder, Functional, Transect) |>
	summarize(Everything = n())

total <- close |>
	filter(Meter < 11) |>
	filter(Invasive == 'Invasive') |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- close |>
	filter(Meter < 11) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

functional <- merge(functional, total, all = TRUE)
functional[is.na(functional)] <- 0

functional <- functional |> filter(Functional != 'Native')

functional <- functional |>
	pivot_wider(names_from = 'Functional', values_from = 'Count')
functional[is.na(functional)] <- 0

functional <- functional |>
	pivot_longer(cols = 7:9, names_to = 'Functional', values_to = 'Successes')

functional <- functional |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Functional, Successes,
								Failures, Total)

library(broom)
options(scipen=999)
# binomial mods
browsed <- functional |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
browsed.means <- emmeans(browsed.mod, ~Feeder*Time, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

ruderals <- functional |>
	filter(Functional == 'Ruderal')
ruderals.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = ruderals, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
ruderal.means <- emmeans(ruderals.mod, ~Feeder*Time, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderals'

resistant <- functional |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
check_model(resistant.mod)
resistant.means <- emmeans(resistant.mod, ~Feeder*Time, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant perennials'

close.means <- rbind(browsed.means, resistant.means, ruderal.means)
close.means <- close.means |>
	mutate_if(is.numeric, round, digits = 3)
close.means$Transect <- "Close"
close.means <- close.means |>
	select(Transect, Functional, Time, Feeder, prob, everything())

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

# Ruderals are significant throughout
# Browsed hammered close to feeder
# Resistant gone flip flop from feeder to non feeder

## --------------- BINOMIAL FAR ------------------------------------------------

# Try it successes and failures
# Clear the decks
rm(list=ls()[! ls() %in% c("close")])

# Count total number of invasive rows
# Get 1s and 0s for each transect point
# Tally natives and invasives at each plot*time

# Tally natives and invasives at each plot*time
sites <- close |>
	group_by(Site, Plot, Time, Feeder, Functional, Transect) |>
	summarize(Everything = n())

total <- close |>
	filter(Meter %in% 11:25) |>
	filter(Invasive == 'Invasive') |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- close |>
	filter(Meter %in% 11:25) |>
	filter(Invasive == 'Invasive') |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

functional <- merge(functional, total, all = TRUE)
functional[is.na(functional)] <- 0

functional <- functional |> filter(Functional != 'Native')

functional <- functional |>
	pivot_wider(names_from = 'Functional', values_from = 'Count')
functional[is.na(functional)] <- 0

functional <- functional |>
	pivot_longer(cols = 7:9, names_to = 'Functional', values_to = 'Successes')

functional <- functional |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Functional, Successes,
								Failures, Total)

library(broom)
options(scipen=999)
# binomial mods
browsed <- functional |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
browsed.means <- emmeans(browsed.mod, ~Feeder*Time, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

ruderals <- functional |>
	filter(Functional == 'Ruderal')
ruderals.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Site/Plot),
									 data = ruderals, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
summary(ruderals.mod)
ruderal.means <- emmeans(ruderals.mod, ~Feeder*Time, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderals'

resistant <- functional |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
summary(resistant.mod)
check_model(resistant.mod)
resistant.means <- emmeans(resistant.mod, ~Feeder*Time, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant perennials'

close.means <- rbind(browsed.means, resistant.means, ruderal.means)
close.means <- close.means |>
	mutate_if(is.numeric, round, digits = 3)
close.means$Transect <- "Close"
close.means <- close.means |>
	select(Transect, Functional, Time, Feeder, prob, everything())

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
	ggtitle('Probability invasive plant is functional group 11-25 m')+
	theme_bw()

## --------------- NATIVE VS INVASIVE <10 --------------------------------------

rm(list=ls()[! ls() %in% c("close")])

# Count total number of invasive rows
# Get 1s and 0s for each transect point
# Tally natives and invasives at each plot*time

# Tally natives and invasives at each plot*time
sites <- close |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Everything = n())

total <- close |>
	filter(Meter < 11) |>
	filter(Species != 'Litter') |>
	filter(Species != 'litter') |>
	filter(Species != 'Bare Ground') |>
	filter(Species != 'Coarse Woody Debris') |>
	filter(Species != 'water') |>
	filter(!is.na(Invasive)) |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- close |>
	filter(Meter < 11) |>
	filter(Species != 'Litter') |>
	filter(Species != 'litter') |>
	filter(Species != 'Bare Ground') |>
	filter(Species != 'Coarse Woody Debris') |>
	filter(Species != 'water') |>
	filter(!is.na(Invasive)) |>
	group_by(Plot, Time, Feeder, Transect, Invasive) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = 'Invasive', values_from = 'Count')
functional[is.na(functional)] <- 0
functional <- functional |>
	pivot_longer(cols = 5:6, names_to = 'Invasive', values_to = 'Successes')

functional <- merge(functional, total, all = TRUE)
functional <- functional |>
	filter(!is.na(Invasive)) # not sure why necessary, keep NA

functional <- functional |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Invasive, Successes,
								Failures, Total)

library(broom)
options(scipen=999)
# binomial mods
invasive <- functional |>
	filter(Invasive == 'Invasive')
invasive.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = invasive, family = 'binomial')
Anova(invasive.mod)
check_model(invasive.mod)
invasive.means <- emmeans(invasive.mod, ~Feeder*Time, type = 'response')
invasive.means <- tidy(invasive.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

native <- functional |>
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
	select(Transect, Invasive, Time, Feeder, prob, everything())

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

# Ruderals are significant throughout
# Browsed hammered close to feeder
# Resistant gone flip flop from feeder to non feeder



## --------------- NATIVE VS INVASIVE 10-25 ------------------------------------

rm(list=ls()[! ls() %in% c("close")])

# Count total number of invasive rows
# Get 1s and 0s for each transect point
# Tally natives and invasives at each plot*time

# Tally natives and invasives at each plot*time
sites <- close |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Everything = n())

total <- close |>
	filter(Meter %in% 11:25) |>
	filter(Species != 'Litter') |>
	filter(Species != 'litter') |>
	filter(Species != 'Bare Ground') |>
	filter(Species != 'Coarse Woody Debris') |>
	filter(Species != 'water') |>
	filter(!is.na(Invasive)) |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

functional <- close |>
	filter(Meter %in% 11:25) |>
	filter(Species != 'Litter') |>
	filter(Species != 'litter') |>
	filter(Species != 'Bare Ground') |>
	filter(Species != 'Coarse Woody Debris') |>
	filter(Species != 'water') |>
	filter(!is.na(Invasive)) |>
	group_by(Plot, Time, Feeder, Transect, Invasive) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = 'Invasive', values_from = 'Count')
functional[is.na(functional)] <- 0
functional <- functional |>
	pivot_longer(cols = 5:6, names_to = 'Invasive', values_to = 'Successes')

functional <- merge(functional, total, all = TRUE)
functional <- functional |>
	filter(!is.na(Invasive)) # not sure why necessary, keep NA

functional <- functional |>
	mutate(Failures = Total-Successes) |>
	dplyr::select(Site, Plot, Time, Feeder, Transect, Invasive, Successes,
								Failures, Total)

library(broom)
options(scipen=999)
# binomial mods
invasive <- functional |>
	filter(Invasive == 'Invasive')
invasive.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = invasive, family = 'binomial')
Anova(invasive.mod)
check_model(invasive.mod)
invasive.means <- emmeans(invasive.mod, ~Feeder*Time, type = 'response')
invasive.means <- tidy(invasive.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

native <- functional |>
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
	select(Transect, Invasive, Time, Feeder, prob, everything())

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

# Ruderals are significant throughout
# Browsed hammered close to feeder
# Resistant gone flip flop from feeder to non feeder

## --------------- PROBABILITY OF EVERYTHING <10 -------------------------------

rm(list=ls()[! ls() %in% c("close")])

for(i in 1:nrow(close)){
	if(is.na(close$Functional[i]) & isTRUE(close$Invasive[i] == 'Native')){
		close$Functional[i] <- 'Native'
	}
}

for(i in 1:nrow(close)){
	if(close$Species[i] == 'Litter' | close$Species[i] == 'litter' |
		 close$Species[i] == 'Bare Ground'){
		close$Functional[i] <- 'Bare ground'
	}
}

functional <- close |>
	filter(!is.na(Functional)) |>
	filter(Meter < 11) |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = 'Functional', values_from = 'Count')
functional[is.na(functional)] <- 0

functional$Total <- 11

bareground <- functional[, c(1,2,3,4,5,10)]
colnames(bareground)[5] <- 'Successes'
bareground <- bareground |>
	mutate(Failures = Total-Successes)

native <- functional[, c(1,2,3,4,6,10)]
colnames(native)[5] <- 'Successes'
# Cap out at 11
for(i in 1:nrow(native)){
	if(native$Successes[i] > 11){
		native$Successes[i] <- 11
	}
}
native <- native |>
	mutate(Failures = Total-Successes)

browsed <- functional[, c(1,2,3,4,7,10)]
colnames(browsed)[5] <- 'Successes'
browsed <- browsed |>
	mutate(Failures = Total-Successes)

resistant <- functional[, c(1,2,3,4,8,10)]
colnames(resistant)[5] <- 'Successes'
resistant <- resistant |>
	mutate(Failures = Total-Successes)

ruderal <- functional[, c(1,2,3,4,9,10)]
colnames(ruderal)[5] <- 'Successes'
ruderal <- ruderal |>
	mutate(Failures = Total-Successes)

library(broom)
options(scipen=999)
# binomial mods

browsed.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
browsed.means <- emmeans(browsed.mod, ~Feeder*Time, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

ruderals.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = ruderal, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
ruderal.means <- emmeans(ruderals.mod, ~Feeder*Time, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderals'

resistant.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
check_model(resistant.mod)
resistant.means <- emmeans(resistant.mod, ~Feeder*Time, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant perennials'

native.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = native, family = 'binomial')
Anova(native.mod)
check_model(native.mod)
native.means <- emmeans(native.mod, ~Feeder*Time, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Functional <- 'Native'

bareground.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = bareground, family = 'binomial')
Anova(bareground.mod)
check_model(bareground.mod)
bareground.means <- emmeans(bareground.mod, ~Feeder*Time, type = 'response')
bareground.means <- tidy(bareground.means, conf.int = TRUE)
bareground.means$Functional <- 'Bareground'

close.means <- rbind(browsed.means, resistant.means, ruderal.means,
										 native.means, bareground.means)
close.means <- close.means |>
	mutate_if(is.numeric, round, digits = 3)
close.means$Transect <- "Close"
close.means <- close.means |>
	select(Transect, Functional, Time, Feeder, prob, everything())

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

close.means |>
	filter(Functional == 'Bareground' | Functional == 'Native') |>
	ggplot(aes(x = as_factor(Time), y = prob, group = interaction(Feeder,Time),
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

## --------------- PROBABILITY OF EVERYTHING 10-25 -----------------------------

rm(list=ls()[! ls() %in% c("close")])

for(i in 1:nrow(close)){
	if(is.na(close$Functional[i]) & isTRUE(close$Invasive[i] == 'Native')){
		close$Functional[i] <- 'Native'
	}
}

for(i in 1:nrow(close)){
	if(close$Species[i] == 'Litter' | close$Species[i] == 'litter' |
		 close$Species[i] == 'Bare Ground'){
		close$Functional[i] <- 'Bare ground'
	}
}

functional <- close |>
	filter(!is.na(Functional)) |>
	filter(Meter %in% 11:25) |>
	group_by(Plot, Time, Feeder, Transect, Functional) |>
	summarize(Count = n())

functional <- functional |>
	pivot_wider(names_from = 'Functional', values_from = 'Count')
functional[is.na(functional)] <- 0

functional$Total <- 15

bareground <- functional[, c(1,2,3,4,5,10)]
colnames(bareground)[5] <- 'Successes'
bareground <- bareground |>
	mutate(Failures = Total-Successes)

native <- functional[, c(1,2,3,4,6,10)]
colnames(native)[5] <- 'Successes'
# Cap out at 11
for(i in 1:nrow(native)){
	if(native$Successes[i] > 11){
		native$Successes[i] <- 11
	}
}
native <- native |>
	mutate(Failures = Total-Successes)

browsed <- functional[, c(1,2,3,4,7,10)]
colnames(browsed)[5] <- 'Successes'
browsed <- browsed |>
	mutate(Failures = Total-Successes)

resistant <- functional[, c(1,2,3,4,8,10)]
colnames(resistant)[5] <- 'Successes'
resistant <- resistant |>
	mutate(Failures = Total-Successes)

ruderal <- functional[, c(1,2,3,4,9,10)]
colnames(ruderal)[5] <- 'Successes'
ruderal <- ruderal |>
	mutate(Failures = Total-Successes)

library(broom)
options(scipen=999)
# binomial mods

browsed.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
browsed.means <- emmeans(browsed.mod, ~Feeder*Time, type = 'response')
browsed.means <- tidy(browsed.means, conf.int = TRUE)
browsed.means$Functional <- 'Browsed perennials'

ruderals.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = ruderal, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
ruderal.means <- emmeans(ruderals.mod, ~Feeder*Time, type = 'response')
ruderal.means <- tidy(ruderal.means, conf.int = TRUE)
ruderal.means$Functional <- 'Ruderals'

resistant.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
check_model(resistant.mod)
resistant.means <- emmeans(resistant.mod, ~Feeder*Time, type = 'response')
resistant.means <- tidy(resistant.means, conf.int = TRUE)
resistant.means$Functional <- 'Resistant perennials'

native.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = native, family = 'binomial')
Anova(native.mod)
check_model(native.mod)
native.means <- emmeans(native.mod, ~Feeder*Time, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Functional <- 'Native'

bareground.mod <- glmer(cbind(Successes,Failures) ~ Feeder*as_factor(Time) + (1|Plot),
									 data = bareground, family = 'binomial')
Anova(bareground.mod)
check_model(bareground.mod)
bareground.means <- emmeans(bareground.mod, ~Feeder*Time, type = 'response')
bareground.means <- tidy(bareground.means, conf.int = TRUE)
bareground.means$Functional <- 'Bareground'

close.means <- rbind(browsed.means, resistant.means, ruderal.means,
										 native.means, bareground.means)
close.means <- close.means |>
	mutate_if(is.numeric, round, digits = 3)
close.means$Transect <- "Close"
close.means <- close.means |>
	select(Transect, Functional, Time, Feeder, prob, everything())

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

close.means |>
	filter(Functional == 'Bareground' | Functional == 'Native') |>
	ggplot(aes(x = as_factor(Time), y = prob, group = interaction(Feeder,Time),
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


## --------------- LINKING BAREGROUND TO INVASIVE ------------------------------

# find all transect points that are 0 that become a 1 for ruderals
# impact of distance on whether or not that becomes a 1

# same for bare ground?

rm(list=ls()[! ls() %in% c("close")])

meters <- close |>
	filter(Feeder == 'Y') |>
	dplyr::select(Site, Plot, Time, Transect) |>
	unique()

meters[c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
				'11', '12', '13', '14', '15', '16', '17', '18', '19',
				'20', '21', '22', '23', '24', '25')] <- 0

meters <- meters |>
	pivot_longer(cols = 5:29, names_to = 'Distance', values_to = 'Placeholder')

ruderal <- close |>
	filter(Functional == 'Ruderal') |>
	filter(Feeder == 'Y') |>
	dplyr::select(Site, Plot, Time, Transect, Meter)
colnames(ruderal)[5] <- "Distance"
ruderal$Present <- 1

ruderal <- merge(meters, ruderal, all = TRUE)
rm(meters)

ruderal[is.na(ruderal)] <- 0
ruderal <- ruderal |>
	mutate(Binary = Placeholder + Present) |>
	dplyr::select(-Present, -Placeholder)

ruderal$Distance <- as.numeric(ruderal$Distance)
ruderal$Transect <- as_factor(ruderal$Transect)
ruderal$Time <- as_factor(ruderal$Time)

plot(ruderal$Distance, ruderal$Binary)

baseline <- ruderal |> filter(Time == '0')
colnames(baseline)[6] <- 'Baseline'
baseline <- baseline |> dplyr::select(-Time)

year1 <- ruderal |> filter(Time == '1')
colnames(year1)[6] <- 'Year1'
year1 <- year1 |> dplyr::select(-Time)

year2 <- ruderal |> filter(Time == '2')
colnames(year2)[6] <- 'Year2'
year2 <- year2 |> dplyr::select(-Time)

ruderal <- merge(baseline, year1, all = TRUE)
ruderal <- merge(ruderal, year2, all = TRUE)
ruderal <- ruderal |> na.omit()

rm(list=ls()[! ls() %in% c("close", "ruderal")])

ruderal <- ruderal |>
	filter(Baseline < 1) |>
	mutate(Total = Baseline + Year1 + Year2) |>
	filter(Total > 0) |>
	dplyr::select(-Total)

ruderal <- ruderal |>
	dplyr::select(-Baseline) |>
	pivot_longer(cols = 5:6, names_to = "Time", values_to = "Binary")

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

ruderal <- ruderal[,-5]
ruderal <- unique(ruderal)

ruderal.sum <- ruderal |>
	group_by(Plot, Bin) |>
	summarize(Count = sum(Binary))

descdist(ruderal.sum$Count, discrete = TRUE)
ruderal.sum$Bin <- factor(ruderal.sum$Bin, order = TRUE,
													levels = c("1-5", "6-10", "11-15",
																		 "16-20", "21-25"))

ggplot(ruderal.sum)+
  aes(x = Bin, y = Count)+
  geom_jitter(width = 0.1, size = 4, shape = 21, fill = 'skyblue',
  						stroke = 1, alpha = 0.5)+
	ylab('Detections')+
	xlab('Distance from feeder (m)')+
  theme_bw()


kruskal.test(Count ~ Bin, data = ruderal.sum)
