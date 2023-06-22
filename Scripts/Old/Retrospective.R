## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse)
library(lme4)
library(car)
library(emmeans)

rm(list=ls())

d <- read.csv('Raw-data/retrospective.csv')
d <- d[1:8659,]

## --------------- Prob of detecting invasive if transect point is plant -------

# Subset the data
d <- d |>
	filter(Meter < 6)
	
# Probability of introduced when detecting a plant
Totals <- d |>
	filter(Invasive != 'na') |>
	group_by(Plot, Time, Treatment) |>
	summarize(Total = n())
# Should be done as a proportion 
# Doesn't look like much 

Native <- d |>
	filter(Invasive == 'Native') |>
	group_by(Plot, Time, Treatment) |>
	summarize(Native = n())

Prob.df <- merge(Totals, Native) |>
	mutate(Introduced = Total-Native,
				 Proportion = Introduced/Total)

m1 <- glmer(cbind(Prob.df$Introduced, Prob.df$Native) ~ Treatment * Time + (1|Plot),
						data = Prob.df, family = binomial)
Anova(m1)

emmeans(m1, pairwise~Treatment, type = 'response')
emtrends(m1, pairwise~Treatment, var = 'Time', type = ' response')

# There is a difference in the probability of detecting introduced species
# but probably less counts overall due to trampling
emmip(m1, Treatment ~ Time, cov.reduce = range)

# Opposite effect than expected 

ggplot(Prob.df, aes(x = Time, y = Proportion, color = Treatment))+
	geom_point()


## --------------- Prob of detecting fruit invasive if transect point is plant -
d <- read.csv('Raw-data/retrospective.csv')
d <- d[1:8659,]


Intro.animal <- d |>
	filter(Species %in% c('japanese honeysuckle', 'sericea lespedeza',
											'chinese privet'))

# Could also go probability that an ANIMAL dispersed species encountered is non-native

d$Status <- NA
for(i in 1:nrow(d)){
	if(d$Invasive[i] == 'Invasive' & d$Species[i] %in% c('japanese honeysuckle', 
							'sericea lespedeza', 'chinese privet')){
		d$Status[i] <- 'Yes'
	} else {
		d$Status[i] <- 'No'
	}
}

# Subset the data
d <- d |>
	filter(Meter < 6)

# Probability of introduced when detecting a plant
Totals <- d |>
	filter(Invasive != 'na') |>
	group_by(Plot, Time, Treatment) |>
	summarize(Total = n())
# Should be done as a proportion 
# Doesn't look like much 

Other <- d |>
	filter(Invasive != 'na' & Status == 'No') |>
	group_by(Plot, Time, Treatment) |>
	summarize(Other = n())

Prob.df <- merge(Totals, Other) |>
	mutate(Introduced.anml = Total-Other) |>
	mutate(Proportion = Introduced.anml/Total)

m1 <- glmer(cbind(Prob.df$Introduced, Prob.df$Total) ~ Treatment*Time + (1|Plot),
						data = Prob.df, family = binomial)
Anova(m1)

emmeans(m1, pairwise~Treatment,type = ' response')
emtrends(m1, pairwise~Treatment, var = 'Time', type = ' response')

# There is a difference in the probability of detecting introduced species
# but probably less counts overall due to trampling
emmip(m1, Treatment ~ Time, cov.reduce = range)

ggplot(Prob.df, aes(x = Time, y = Proportion, color = Treatment))+
	geom_point()

## --------------- Prob of detecting invasive ----------------------------------

rm(list=ls())

d <- read.csv('Raw-data/retrospective.csv')
d <- d[1:8659,]


# Subset the data
d <- d |>
	filter(Meter < 6)
	
# Probability of introduced when detecting a plant
Totals <- d |>
	group_by(Plot, Time, Treatment) |>
	summarize(Total = n())

# Proportions screwed up by double hits [plants increase total]

Native <- d |>
	filter(Invasive == 'Native') |>
	group_by(Plot, Time, Treatment) |>
	summarize(Native = n())

Prob.df <- merge(Totals, Native) |>
	mutate(Introduced = Total-Native) |>
	mutate(Proportion = Introduced/Total)

m1 <- glmer(cbind(Prob.df$Introduced, Prob.df$Total) ~ Treatment +
							as_factor(Time) + (1|Plot),
						data = Prob.df, family = binomial)
Anova(m1)

emmeans(m1, pairwise~Treatment*Time,type = ' response',
				adjust = 'none')
emtrends(m1, pairwise~Treatment, var = 'Time', type = ' response')

# There is a difference in the probability of detecting introduced species
# but probably less counts overall due to trampling
emmip(m1, Treatment ~ Time, cov.reduce = range)

# Opposite effect than expected 

ggplot(Prob.df, aes(x = Time, y = Proportion, color = Treatment))+
	geom_point()

## -------- Prob of detecting invasive if transect point is plant group years --

rm(list=ls())

d <- read.csv('Raw-data/retrospective.csv')
d <- d[1:8659,]

# Subset the data
d <- d |>
	filter(Meter < 6)

# Year
d$Period <- NA
First <- c(2,3,4)
for(i in 1:nrow(d)){
	if(d$Time[i] %in% First){
		d$Period[i] <- 'First'
	} else {
		d$Period[i] <- 'Second'
	}
}

# Probability of introduced when detecting a plant
Totals <- d |>
	group_by(Plot, Period, Treatment) |>
	summarize(Total = n())
# Should be done as a proportion 
# Doesn't look like much 

Plants <- d |>
	filter(Invasive != 'na') |>
	group_by(Plot, Period, Treatment) |>
	summarize(Plants = n())

Native <- d |>
	filter(Invasive == 'Native') |>
	group_by(Plot, Period, Treatment) |>
	summarize(Native = n())

Prob.df <- merge(Totals, Plants, all = TRUE)
Prob.df[is.na(Prob.df)] <- 0

Prob.df <- merge(Prob.df, Native, all = TRUE) 
Prob.df[is.na(Prob.df)] <- 0

Prob.df <- Prob.df |>
	mutate(Introduced = Total-Native) |>
	mutate(Proportion = Introduced/Total)

Prob.df$Period <- as_factor(Prob.df$Period)

Prob.df <- Prob.df |>
	filter(Plants > 0)

m1 <- glmer(cbind(Prob.df$Introduced, Prob.df$Total) ~ Treatment * Period + (1|Plot),
						data = Prob.df, family = binomial)
Anova(m1)

emmeans(m1, ~Treatment*Period, type = 'response')
emtrends(m1, pairwise~Treatment, var = 'Time', type = ' response')

# There is a difference in the probability of detecting introduced species
# but probably less counts overall due to trampling
emmip(m1, Treatment ~ Time, cov.reduce = range)

# Opposite effect than expected 
ggplot(Prob.df, aes(x = Period, y = Proportion, color = Treatment))+
	geom_point()


# We could also make it just a binary and get ride of the multiples, which is 
# confusing 


## --------------- Prob of detecting animal-dispersed if transect point is plant 

rm(list=ls())

d <- read.csv('Raw-data/retrospective.csv')
d <- d[1:8659,]

# Subset the data
d <- d |>
	filter(Meter < 6)
	
# Probability of introduced when detecting a plant

Totals <- d |>
	group_by(Plot, Time, Treatment) |>
	summarize(Total = n())
# Should be done as a proportion 

Plants <- d |>
	filter(Veg.Type != 'na') |>
	group_by(Plot, Time, Treatment) |>
	summarize(Plants = n())
# Should be done as a proportion 
# Doesn't look like much 

endo.ls <- c('muscadine', 'american beautyberry', 'virginia creeper',
					'blackberry', 'smilax', 'sparkleberry', 'peppervine', 
					'poison ivy', 'eastern red cedar', 'possumhaw', 'deerberry',
					'yaupon', 'winged sumac', 'snakeroot')

epi.ls <- c('bluntleaf bedstraw', 'desmodium ', 'spanish needles')

dys.ls <- c('Shumard Oak', 'post oak ', 'willow oak', 'turkey oak')

comb <- c(endo.ls, epi.ls, dys.ls)

Native.animal.dispersed <- d |>
	filter(Invasive == 'Native' & Species %in% comb)|>
	group_by(Plot, Time, Treatment) |>
	summarize(Native = n())

Prob.df <- merge(Totals, Plants, all = TRUE) 
Prob.df <- merge(Prob.df, Native.animal.dispersed, all = TRUE) 
Prob.df[is.na(Prob.df)] <- 0

Prob.df <- Prob.df |>
	mutate(Other = Plants-Native) |>
	mutate(Misses = Total-Native) |>
	mutate(Proportion = Native/Total)
Prob.df[is.na(Prob.df)] <- 0

m1 <- glmer(cbind(Prob.df$Native, Prob.df$Misses) ~ Treatment + (1|Time) + (1|Plot),
						data = Prob.df, family = binomial)
Anova(m1)

emmeans(m1, pairwise~Treatment, type = 'response')
emtrends(m1, pairwise~Treatment, var = 'Time', type = ' response')

# There is a difference in the probability of detecting introduced species
# but probably less counts overall due to trampling
emmip(m1, Treatment ~ Time, cov.reduce = range)

# Opposite effect than expected 

ggplot(Prob.df, aes(x = Time, y = Proportion, color = Treatment))+
	geom_point()

# Native are more common away from feeders

