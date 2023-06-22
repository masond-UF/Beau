library(tidyverse)
library(performance)
library(DHARMa)
library(fitdistrplus)

#

rm(list=ls())

# transects
d <- read.csv('Clean-data/Feeder-veg.csv') |> 
	filter(Dataset == 'Manipulative')

# plots
d <- read.csv('Raw-data/manipulative-plots.csv')

####################### Top species counts #####################################

d |>
	filter(Invasive == 'Invasive') |>
	group_by(Species) |>
	summarize(Count = n()) |>
	arrange(desc(Count)) |>
	top_n(10)

list <- c('Japanese Honeysuckle', 'Moss', 'Poison Ivy', 'Blackberry',
				 'Goldenrod', 'Virginia Creeper', 'Muscadine', 'Snakeroot',
				 'Trumpet Creeper', 'Dicahnthelium')

top.sum <- d |>
	filter(Species %in% list) |>
	group_by(Plot, Feeder, Time, Species) |>
	summarize(Counts = n())

for(i in 1:length(list)){
	spec <- top.sum |> filter(Species == list[i])
	fig <- ggplot(spec, aes(x = Time, y = Counts, color = Feeder))+
						geom_jitter() +
						ggtitle(list[i])
	print(fig)
}

top.mean <- top.sum |>
	group_by(Feeder, Time, Species) |>
	summarize(Mean = mean(Counts),
						n = n(),
						sd = sd(Counts),
						se = sd/sqrt(n))

for(i in 1:length(list)){
	spec <- top.mean |> filter(Species == list[i])
	fig <- ggplot(spec, aes(x = Time, y = Mean, color = Feeder))+
						geom_errorbar(aes(ymin = Mean-se, ymax = Mean+se),
													width=.2, position=position_dodge(.9))+
						geom_point(position = position_dodge(width =0.9)) +
						ggtitle(list[i])
	print(fig)
}

# Pivot wide to introduce 0s and pivot back
top.sum <- top.sum |>
	pivot_wider(names_from = Species, values_from = Counts)
top.sum[is.na(top.sum)] <- 0

top.sum <- top.sum |>
	pivot_longer(4:25, names_to = 'Species', values_to = 'Counts')

top.sum <- top.sum |>
	ungroup() |>
	tidyr::complete(Plot, Feeder, Time, Species)

Baseline <- top.sum |> filter(Time == 0)
Year1 <- top.sum |> filter(Time == 1)
Year2 <- top.sum |> filter(Time == 2)

colnames(Baseline)[5] <- 'Baseline'
colnames(Year1)[5] <- 'Year1'
colnames(Year2)[5] <- 'Year2'

Year1 <- Year1[,5]
Year2 <- Year2[,5]

comb <- cbind(Baseline, Year1, Year2)

comb <- comb |>
	mutate(Year1.change = Year1 - Baseline,
				 Year2.change = Year2 - Baseline)

comb <- comb |>
	dplyr::select(-Baseline, -Year1, -Year2)

comb <- comb |> 
	pivot_longer(5:6, names_to = 'Year', values_to = 'Change')

comb.sum <- comb |>
	drop_na() |>
	group_by(Species, Year, Feeder) |>
	summarize(Mean.change = mean(Change),
						sd = sd(Change),
						n = n(),
						se = sd/sqrt(n))
	
for(i in 1:length(list)){
	spec <- comb.sum |> filter(Species == list[i])
	fig <- ggplot(spec, aes(x = Year, y = Mean.change, color = Feeder))+
						geom_point() +
						ggtitle(list[i])
	print(fig)
}

# Effect size
comb.feeder <- comb |> filter(Feeder == 'Y')
comb.feeder <- comb.feeder |> dplyr::select(-Feeder)
colnames(comb.feeder)[5] <- 'Feeder.change'

comb.no <- comb |> filter(Feeder == 'N')
comb.no <- comb.no |> dplyr::select(-Feeder)
colnames(comb.no)[5] <- 'No.Feeder.change'

comb.no <- comb.no[,5]
comb.es <- cbind(comb.feeder, comb.no)

comb.es <- comb.es |>
	mutate(LR = log((Feeder.change+15)/(No.Feeder.change+15)))

comb.es.sum <- comb.es |>
	drop_na() |>
	group_by(Species, Year) |>
	summarize(Mean.LR = mean(LR),
						sd = sd(LR),
						n = n(),
						se = sd/sqrt(n))

for(i in 1:length(list)){
	spec <- comb.es.sum |> filter(Species == list[i])
	fig <- ggplot(spec, aes(x = Year, y = Mean.LR))+
						scale_y_continuous(limits = c(-0.2,0.2))+
						geom_point() +
						ggtitle(list[i])
	print(fig)
}

####################### Top species height #####################################

rm(list=ls())

d <- read.csv('Raw-data/Short-term.csv')


list <- c("Vine", "Forb", "Tree", "Grass","Shrub", "Moss")

d$Height <- as.factor(d$Height)
d$Height <- as.numeric(d$Height)

height.mean <- d |>
	# filter(Species %in% list) |>
	group_by(Plot, Feeder, Time, Veg.Type) |>
	summarize(Mean.ht = mean(Height))

height.mean <- height.mean |>
	ungroup() |>
	tidyr::complete(Plot, Feeder, Time, Veg.Type)

Baseline <- height.mean |> filter(Time == '2018')
Year1 <- height.mean |> filter(Time == '2019')
Year2 <- height.mean |> filter(Time == '2020')

colnames(Baseline)[5] <- 'Baseline'
colnames(Year1)[5] <- 'Year1'
colnames(Year2)[5] <- 'Year2'

Year1 <- Year1[,5]
Year2 <- Year2[,5]

comb <- cbind(Baseline, Year1, Year2)

comb <- comb |>
	mutate(Year1.change = Year1 - Baseline,
				 Year2.change = Year2 - Baseline)

comb <- comb |>
	dplyr::select(-Baseline, -Year1, -Year2)

comb <- comb |> 
	pivot_longer(5:6, names_to = 'Year', values_to = 'Change')

# Effect size
comb.feeder <- comb |> filter(Feeder == 'Y')
comb.feeder <- comb.feeder |> select(-Feeder)
colnames(comb.feeder)[5] <- 'Feeder.change'

comb.no <- comb |> filter(Feeder == 'N')
comb.no <- comb.no |> select(-Feeder)
colnames(comb.no)[5] <- 'No.Feeder.change'

comb.no <- comb.no[,5]
comb.es <- cbind(comb.feeder, comb.no)

comb.es <- comb.es |>
	mutate(LR = log((Feeder.change+100)/(No.Feeder.change+100)))

comb.es.sum <- comb.es |>
	drop_na() |>
	group_by(Veg.Type, Year) |>
	summarize(Mean.LR = mean(LR),
						sd = sd(LR),
						n = n(),
						se = sd/sqrt(n))

for(i in 1:length(list)){
	spec <- comb.es.sum |> filter(Veg.Type == list[i])
	fig <- ggplot(spec, aes(x = Year, y = Mean.LR))+
						scale_y_continuous(limits = c(-0.2,0.2))+
						geom_point() +
						ggtitle(list[i])
	print(fig)
}


list <- c('Japanese Honeysuckle', 'Moss', 'Poison Ivy', 'Blackberry',
				 'Goldenrod', 'Virginia Creeper', 'Muscadine', 'Snakeroot',
				 'Trumpet Creeper', 'Dicahnthelium')

d$Height <- as.factor(d$Height)
d$Height <- as.numeric(d$Height)

top.ht <- d |>
	filter(Species %in% list) |>
	group_by(Plot, Feeder, Time, Species) |>
	summarize(mean.ht = mean(Height),
						sd = sd(Height),
						n = n(),
						se = sd/sqrt(n))

mean.top.ht <- top.ht |>
		group_by(Feeder, Time, Species) |>
		summarize(mean.top.ht = mean(mean.ht),
						sd = sd(mean.ht),
						n = n(),
						se = sd/sqrt(n))

for(i in 1:length(list)){
	spec <- mean.top.ht |> filter(Species == list[i])
	fig <- ggplot(spec, aes(x = Time, y = mean.top.ht,
													color = Feeder))+
						geom_errorbar(aes(ymin = mean.top.ht-se, ymax = mean.top.ht+se),
													width=.2, position=position_dodge(.9))+
						geom_point(position = position_dodge(width =0.9)) +
						ggtitle(list[i])
	print(fig)
} # Feeders impact height. Could be interaction with feeder + rain

top.ht <- d |>
	filter(Species %in% list) 

descdist(top.ht$Height)

m1 <- lmer(Height ~ Feeder * as_factor(Time) * Species + (1|Plot),
					 top.ht)
Anova(m1)
check_model(m1)

emmeans(m1, pairwise~Feeder+Species, type = 'response')

####################### Total count ############################################


total.count <- d |>
	group_by(Plot, Feeder, Time, Invasive) |>
	summarize(Counts = n())

class(total.count$Time)

m1 <- glmer(Counts ~ Feeder * Invasive * as_factor(Time) + (1|Plot),
						data = total.count, family = 'poisson')
Anova(m1)

check_model(m1)
sim <- simulateResiduals(m1)
testZeroInflation(sim)

emmeans(m1, pairwise~Feeder, type = 'response') 
# Less stems at feeders, mostly due to loss of natives...

emmeans(m1, pairwise~Feeder*Invasive, type = 'response') 
# Invasive lose 1.5  with feeder, is different
# Natives gain 1, not different

emmeans(m1, pairwise~Feeder*Time, type = 'response') # Less stems at feeders

####################### Total height ###########################################

descdist(d$Height)

m1 <- lmer(Height ~ Feeder * Invasive * as_factor(Time) + (1|Plot),
						data = d) 
Anova(m1)
# Could include something to reflect density (i.e., incorporate height)
check_model(m1)

emmeans(m1, pairwise~Feeder*Invasive, type = 'response') # Less stems at feeders
emmeans(m1, pairwise~Feeder*Time, type = 'response') # Less stems at feeders

# Feeders favor invasives more than natives


# ------------------ Plots -----------------------------------------------------

plots <- read.csv('Clean-data/manipulative-plots.csv')