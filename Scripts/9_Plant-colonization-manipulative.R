## --------------- HEADER ------------------------------------------------------
## Script name: 7_Plant-colonization-manipulative.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-22
## Date Last modified: 2023-010-06
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: 

## --------------- CREATE BARE GROUND TO RUDERAL DF ----------------------------

# Clear the decks
rm(list=ls())

bare.ground <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 2) |>
	filter(Species == 'Bare Ground' | 
				 	Species == 'Litter' |
				 	Functional == 'Ruderal') |>
	dplyr::select(Plot, Direction, Transect, Meter, Species)
colnames(bare.ground)[5] <- 'Before'

bare.ground.colonized <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 5) |>
	dplyr::select(Plot, Direction, Meter, Functional)

second.sample <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Time == 2) |>
		dplyr::select(Plot, Direction, Meter, Transect)

bare.ground.colonized <- left_join(bare.ground.colonized, unique(second.sample))

bare.ground.colonized <- bare.ground.colonized |>
	dplyr::select(Plot, Transect, Meter, Functional)
colnames(bare.ground.colonized)[4] <- 'After'

for(i in 1:nrow(bare.ground.colonized)){
 if(isTRUE(bare.ground.colonized$After[i] == 'Ruderal')){
			bare.ground.colonized$After[i] <- 1
 } else {
 		  bare.ground.colonized$After[i] <- 0
 }
}

bare.ground <- bare.ground |>
		dplyr::select(Plot, Transect, Meter, Before)

bare.ground.final <- merge(bare.ground, bare.ground.colonized,
													 by = c('Plot', 'Transect', 'Meter'))
rm(list=ls()[! ls() %in% c("bare.ground.final")])

bare.ground.final$After <- as.numeric(bare.ground.final$After)
mod <- glmer(After ~ Meter + (1|Plot),
						 family = 'binomial', data = bare.ground.final)
Anova(mod)

# Less invasives at BG than with Native plants
# No linear pattern with distance
# 

## --------------- CREATE NATIVE TO RUDERAL DF ---------------------------------

not.bg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 2) |>
	filter(Invasive == 'Native') |>
	dplyr::select(Plot, Meter, Transect, Invasive)
colnames(not.bg)[4] <- 'Before'

not.bg.colonized <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 5) |>
	drop_na(Invasive) |>
	dplyr::select(Plot, Direction, Meter, Invasive)

second.sample <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Time == 2) |>
		dplyr::select(Plot, Direction, Meter, Transect)

not.bg.colonized <- left_join(not.bg.colonized, unique(second.sample))

not.bg.colonized <- not.bg.colonized |>
	dplyr::select(Plot, Transect, Meter, Invasive)
colnames(not.bg.colonized)[4] <- 'After'

for(i in 1:nrow(not.bg.colonized)){
	if(isTRUE(not.bg.colonized$After[i] == 'Invasive')){
		not.bg.colonized$After[i] <- 1
	}
		if(isTRUE(not.bg.colonized$After[i] == 'Native')){
			not.bg.colonized$After[i] <- 0
		}
			if(is.na(not.bg.colonized$After[i])){
				not.bg.colonized$After[i] <- 0
	}
}

not.bg.final <- merge(not.bg, not.bg.colonized,
													 by = c('Plot', 'Transect', 'Meter'))

## --------------- COMBINE DF --------------------------------------------------

d <- rbind(bare.ground.final, not.bg.final)

rm(list=ls()[! ls() %in% c("d")])

colnames(d)[4] <- 'Start'
colnames(d)[5] <- 'Invasion'

d$Transect.ID <- paste(d$Plot, d$Direction)

d <- d |> dplyr::select(Plot, Transect, Transect.ID, Meter, Start, Invasion)

d$Invasion <- as.numeric(d$Invasion)

for(i in 1:nrow(d)){
	if(d$Start[i] != 'Native') {
		 d$Start[i] <- 'Ruderal or Bare Ground'
	}
}

## --------------- MODEL INVASION PROBABILITY ----------------------------------

mod <- glmer(Invasion ~ Start * Meter + (1|Plot/Transect),
						 family = 'binomial', data = d)
Anova(mod)
emmeans(mod, ~Start, type = 'response')
emmip(mod, Start ~ Meter, cov.reduce = range,
			nuisance = c('Plot', 'Transect'))

d$Probability <- predict(mod, d[,1:5], type="response")

logit <- function(x) log(x)/log(1-x)

ruderal.fig <- ggplot(d, aes(x = Meter, y = Probability, fill = Start, color = Start)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ylab("Probability of Ruderal Colonization") +
	coord_trans(y="logit") +
  xlab("Distance from feeder (m)") +
	scale_x_continuous(breaks = seq(0,10,1)) +
	scale_color_brewer(palette = "Dark2", name = "") +
	theme_bw() +
	scale_fill_brewer(palette = "Dark2", name = "") +
		theme(panel.grid.minor = element_blank(),
				panel.grid.major = element_blank(),
				legend.position = 'none',
				strip.background = element_blank(),
  			strip.text.x = element_blank(),
				aspect.ratio = 2,
				text=element_text(size=20),
				axis.title = element_text(face="bold"))

## --------------- CREATE BARE GROUND TO BARE GROUND DF ------------------------

rm(list=ls()[! ls() %in% c("ruderal.fig")])

bare.ground <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 2) |>
	filter(Species == 'Bare Ground' | 
				 	Species == 'Litter') |>
	dplyr::select(Plot, Direction, Transect, Meter, Species)
colnames(bare.ground)[5] <- 'Before'

bare.ground.colonized <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 5) |>
	dplyr::select(Plot, Direction, Meter, Species)

second.sample <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Time == 2) |>
		dplyr::select(Plot, Direction, Meter, Transect)

bare.ground.colonized <- left_join(bare.ground.colonized, unique(second.sample))

bare.ground.colonized <- bare.ground.colonized |>
	dplyr::select(Plot, Transect, Meter, Species)
colnames(bare.ground.colonized)[4] <- 'After'

for(i in 1:nrow(bare.ground.colonized)){
 if(isTRUE(bare.ground.colonized$After[i] == 'Bare Ground' |
 					bare.ground.colonized$After[i] == 'Litter')){
			bare.ground.colonized$After[i] <- 1
 } else {
 		  bare.ground.colonized$After[i] <- 0
 }
}

bare.ground <- bare.ground |>
		dplyr::select(Plot, Transect, Meter, Before)

bare.ground.final <- merge(bare.ground, bare.ground.colonized,
													 by = c('Plot', 'Transect', 'Meter'))

rm(list=ls()[! ls() %in% c("bare.ground.final", "ruderal.fig")])

bare.ground.final$After <- as.numeric(bare.ground.final$After)
mod <- glmer(After ~ Meter + (1|Plot),
						 family = 'binomial', data = bare.ground.final)
Anova(mod)

bare.ground.final$Before <- 'Bare Ground or Litter'

## --------------- CREATE NATIVE DF --------------------------------------------

not.bg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 2) |>
	filter(Invasive == 'Native') |>
	dplyr::select(Plot, Meter, Transect, Invasive)
colnames(not.bg)[4] <- 'Before'

not.bg.colonized <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 5) |>
	drop_na(Invasive) |>
	dplyr::select(Plot, Direction, Meter, Species)

second.sample <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Time == 2) |>
		dplyr::select(Plot, Direction, Meter, Transect)

not.bg.colonized <- left_join(not.bg.colonized, unique(second.sample))

not.bg.colonized <- not.bg.colonized |>
	dplyr::select(Plot, Transect, Meter, Species)
colnames(not.bg.colonized)[4] <- 'After'

for(i in 1:nrow(not.bg.colonized)){
	if(isTRUE(not.bg.colonized$After[i] == 'Bare Ground' |
 					not.bg.colonized$After[i] == 'Litter')){
		not.bg.colonized$After[i] <- 1
	} else {
		not.bg.colonized$After[i] <- 0
	}
}

not.bg.final <- merge(not.bg, not.bg.colonized,
													 by = c('Plot', 'Transect', 'Meter'))

## --------------- COMBINE DF --------------------------------------------------

d <- rbind(bare.ground.final, not.bg.final)

rm(list=ls()[! ls() %in% c("d", "ruderal.fig")])

colnames(d)[4] <- 'Start'
colnames(d)[5] <- 'Status'

d$Transect.ID <- paste(d$Plot, d$Direction)

d <- d |> dplyr::select(Plot, Transect, Transect.ID, Meter, Start, Status)

d$Status <- as.numeric(d$Status)

unique(d$Start)
class(d$Meter)

## --------------- MODEL BARE GROUND PROBABILITY -------------------------------

mod <- glm(Status ~ Start * Meter + Plot,
						 family = 'binomial', data = d)
Anova(mod)
emmeans(mod, ~Start, type = 'response')
emmip(mod, Start ~ Meter, cov.reduce = range,
			nuisance = c('Plot', 'Transect'))

d$Probability <- predict(mod, d[,1:5], type="response")

logit <- function(x) log(x)/log(1-x)

d <- d |> filter(Start != 'Native')

bare.ground.fig <- ggplot(d, aes(x = Meter, y = Status, fill = Start, color = Start)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ylab("Probability of Bare Ground") +
  xlab("Distance from feeder (m)") +
	scale_x_continuous(breaks = seq(0,10,1)) +
	scale_y_continuous(breaks = seq(0,0.5,0.1)) +
	scale_color_manual(values = c('#D95F02'))+
	scale_fill_manual(values = c('#D95F02'))+
	theme_bw() +
	geom_hline(yintercept = 0, color = '#1B9E77', size = 1) +
		theme(panel.grid.minor = element_blank(),
				panel.grid.major = element_blank(),
				legend.position = 'none',
				strip.background = element_blank(),
  			strip.text.x = element_blank(),
				aspect.ratio = 2,
				text=element_text(size=20),
				axis.title = element_text(face="bold"))

library(RColorBrewer)
brewer.pal(n=5,"Dark2")
"#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E"


## --------------- CREATE COMBINED FIGURE --------------------------------------

library(patchwork)

bare.ground.fig + ruderal.fig

ggsave(filename = 'Figures/Manipulative-colonization.png')
