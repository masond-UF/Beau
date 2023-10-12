## --------------- BINARY FULL INVASIVE/NATIVE  --------------------------------

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	dplyr::select(Plot, Transect, Meter, Feeder, Invasive) |>
	unique() |>
	group_by(Plot, Transect, Feeder, Invasive) |>
	summarize(Successes = n())

# Introduce zeroes for missing values
feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Invasive', values_from = 'Successes')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0
feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 4:6, names_to = 'Invasive', values_to = 'Successes')

# Check values
v1 <- nrow(feeder.veg.sum) / 3 # functional groups
v2 <- v1 / 2 # feeders
v3 <- v2/4 # transects per plot
isTRUE(length(unique(feeder.veg$Plot)) == v3) # Correct

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
feeder.veg.sum$Total <- 25 # Number of transect points
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total - Successes) |>
	dplyr::select(Plot, Transect, Feeder, Invasive, Total, Successes, Failures)

# Invasive
invasive <- feeder.veg.sum |>
	filter(Invasive == 'Invasive')

mod.invasive <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = invasive)	
Anova(mod.invasive)
summary(mod.invasive)

invasive.means <- emmeans(mod.invasive, ~Feeder, type = 'response')
invasive.means <- tidy(invasive.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

# Native
native <- feeder.veg.sum |>
	filter(Invasive == 'Native')

mod.native <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = native)	
Anova(mod.native)
summary(mod.native)
native.means <- emmeans(mod.native, ~Feeder, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

# Not plant
not.plant <- feeder.veg.sum |>
	filter(Invasive == 'Not plant')

mod.not.plant <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = not.plant)	
Anova(mod.not.plant)
summary(mod.not.plant)
not.plant.means <- emmeans(mod.not.plant, ~Feeder, type = 'response')
not.plant.means <- tidy(not.plant.means, conf.int = TRUE)
not.plant.means$Invasive <- 'Not plant'


## --------------- BINARY CLOSE INVASIVE/NATIVE --------------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	filter(Meter < 11) |>
	dplyr::select(Plot, Transect, Meter, Feeder, Invasive) |>
	unique() |>
	group_by(Plot, Transect, Feeder, Invasive) |>
	summarize(Successes = n())

# Introduce zeroes for missing values
feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Invasive', values_from = 'Successes')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0
feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 4:6, names_to = 'Invasive', values_to = 'Successes')


# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
feeder.veg.sum$Total <- 11 # Number of transect points
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total - Successes) |>
	dplyr::select(Plot, Transect, Feeder, Invasive, Total, Successes, Failures)

# Invasive
invasive <- feeder.veg.sum |>
	filter(Invasive == 'Invasive')
mod.invasive <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = invasive)	
Anova(mod.invasive)
summary(mod.invasive)
invasive.means <- emmeans(mod.invasive, ~Feeder, type = 'response')
invasive.means <- tidy(invasive.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

# Native
native <- feeder.veg.sum |>
	filter(Invasive == 'Native')
mod.native <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = native)	
Anova(mod.native)
summary(mod.native)
native.means <- emmeans(mod.native, ~Feeder, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Invasive <- 'Native'

# Not plant
not.plant <- feeder.veg.sum |>
	filter(Invasive == 'Not plant')

mod.not.plant <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = not.plant)	
Anova(mod.not.plant)
summary(mod.not.plant)
not.plant.means <- emmeans(mod.not.plant, ~Feeder, type = 'response')
not.plant.means <- tidy(not.plant.means, conf.int = TRUE)
not.plant.means$Invasive <- 'Not plant'


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



## --------------- PROPORTION INVASIVE/NATIVE ----------------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Tally natives and invasives at each plot*time
sites <- feeder.veg |>
	group_by(Plot, Feeder, Transect) |>
	summarize(Everything = n())

total <- feeder.veg |>
	filter(!Invasive == 'Not plant') |>
	group_by(Plot, Feeder, Transect) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	dplyr::select(Plot, Transect, Meter, Feeder, Invasive) |>
	group_by(Plot, Transect, Feeder, Invasive) |>
	summarize(Successes = n())

# Introduce zeroes for missing values
feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Invasive', values_from = 'Successes')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0
feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 4:6, names_to = 'Invasive', values_to = 'Successes')

# Check values
v1 <- nrow(feeder.veg.sum) / 3 # functional groups
v2 <- v1 / 2 # feeders
v3 <- v2/4 # transects per plot
isTRUE(length(unique(feeder.veg$Plot)) == v3) # Correct

feeder.veg.sum <- merge(feeder.veg.sum, total, all = TRUE)
feeder.veg.sum <- feeder.veg.sum |>
	filter(!Invasive == 'Not plant') # not sure why necessary, keep NA

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total - Successes) |>
	dplyr::select(Plot, Transect, Feeder, Invasive, Total, Successes, Failures)

# Invasive
invasive <- feeder.veg.sum |>
	filter(Invasive == 'Invasive')

mod.invasive <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = invasive)	
Anova(mod.invasive)
summary(mod.invasive)
invasive.means <- emmeans(mod.invasive, ~Feeder, type = 'response')
invasive.means <- tidy(invasive.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

# Native
native <- feeder.veg.sum |>
	filter(Invasive == 'Native')

mod.native <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = native)	
Anova(mod.native)
summary(mod.native)
native.means <- emmeans(mod.native, ~Feeder, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Invasive <- 'Native'



## --------------- CLOSE PROPORTION INVASIVE/NATIVE ----------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Tally natives and invasives at each plot*time
sites <- feeder.veg |>
	group_by(Plot, Feeder, Transect) |>
	summarize(Everything = n())

total <- feeder.veg |>
	filter(!Invasive == 'Not plant') |>
	group_by(Plot, Feeder, Transect) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	filter(Meter < 11) |>
	dplyr::select(Plot, Transect, Meter, Feeder, Invasive) |>
	group_by(Plot, Transect, Feeder, Invasive) |>
	summarize(Successes = n())

# Introduce zeroes for missing values
feeder.veg.sum <- feeder.veg.sum |>
	pivot_wider(names_from = 'Invasive', values_from = 'Successes')
feeder.veg.sum[is.na(feeder.veg.sum)] <- 0
feeder.veg.sum <- feeder.veg.sum |>
	pivot_longer(cols = 4:6, names_to = 'Invasive', values_to = 'Successes')

# Check values
v1 <- nrow(feeder.veg.sum) / 3 # functional groups
v2 <- v1 / 2 # feeders
v3 <- v2/4 # transects per plot
isTRUE(length(unique(feeder.veg$Plot)) == v3) # Correct

feeder.veg.sum <- merge(feeder.veg.sum, total, all = TRUE)
feeder.veg.sum <- feeder.veg.sum |>
	filter(!Invasive == 'Not plant') # not sure why necessary, keep NA

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
feeder.veg.sum <- feeder.veg.sum |>
	mutate(Failures = Total - Successes) |>
	dplyr::select(Plot, Transect, Feeder, Invasive, Total, Successes, Failures)

# Invasive
invasive <- feeder.veg.sum |>
	filter(Invasive == 'Invasive')

mod.invasive <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = invasive)	
Anova(mod.invasive)
summary(mod.invasive)
invasive.means <- emmeans(mod.invasive, ~Feeder, type = 'response')
invasive.means <- tidy(invasive.means, conf.int = TRUE)
invasive.means$Invasive <- 'Invasive'

# Native
native <- feeder.veg.sum |>
	filter(Invasive == 'Native')

mod.native <- glmer(cbind(Successes, Failures) ~ Feeder + (1|Plot), 
											family = 'binomial', data = native)	
Anova(mod.native)
summary(mod.native)
native.means <- emmeans(mod.native, ~Feeder, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Invasive <- 'Native'



## --------------- PROPORTION FUNCTIONAL ---------------------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Tally natives and invasives at each plot*time
sites <- feeder.veg |>
	group_by(Plot, Feeder, Transect) |>
	summarize(Everything = n())

total <- feeder.veg |>
	filter(!Invasive == 'Not plant') |>
	group_by(Plot, Feeder, Transect) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	dplyr::select(Plot, Transect, Meter, Feeder, Functional) |>
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

feeder.veg.sum <- merge(feeder.veg.sum, total, all = TRUE)

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
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

# Natives
native <- feeder.veg.sum |>
	filter(Functional == 'Native')

mod.native <- glm(cbind(Successes, Failures) ~ Feeder + Plot, # Not converging
											family = 'binomial', data = native) # with random effects
Anova(mod.native)
summary(mod.native)

native.means <- emmeans(mod.native, ~Feeder, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Functional <- 'Native'

## --------------- CLOSE PROPORTION FUNCTIONAL  --------------------------------

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg")])

# Tally natives and invasives at each plot*time
sites <- feeder.veg |>
	group_by(Plot, Feeder, Transect) |>
	summarize(Everything = n())

total <- feeder.veg |>
	filter(!Invasive == 'Not plant') |>
	group_by(Plot, Feeder, Transect) |>
	summarize(Total = n())

total <- merge(sites, total, all = TRUE)
total[is.na(total)] <- 0
total <- total |> dplyr::select(-Everything)

# Calculate the total number of hits for functional roles
# at each time and transect
feeder.veg.sum <- feeder.veg |>
	filter(Meter < 11) |>
	dplyr::select(Plot, Transect, Meter, Feeder, Functional) |>
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

feeder.veg.sum <- merge(feeder.veg.sum, total, all = TRUE)

# Clear the decks
rm(list = ls()[!ls() %in% c("feeder.veg", "feeder.veg.sum")])

# Calculate the number of failures 
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

mod.resistant <- glm(cbind(Successes, Failures) ~ Feeder + Plot, 
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

# Natives
native <- feeder.veg.sum |>
	filter(Functional == 'Native')

mod.native <- glm(cbind(Successes, Failures) ~ Feeder + Plot, # Not converging
											family = 'binomial', data = native) # with random effects
Anova(mod.native)
summary(mod.native)

native.means <- emmeans(mod.native, ~Feeder, type = 'response')
native.means <- tidy(native.means, conf.int = TRUE)
native.means$Functional <- 'Native'

## --------------- CHECK INVASIVES ---------------------------------------------

spec.list <- feeder.veg |>
	dplyr::select(Invasive, Species) |>
	unique()

