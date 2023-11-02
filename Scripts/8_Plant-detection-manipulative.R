## --------------- HEADER ------------------------------------------------------
## Script name: 4_Plant-detection-manipulative.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-23
## Date Last modified: 2023-06-22
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for analyzing the manipulative feeder
## vegetation data

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

## --------------- TIME 0-2 ----------------------------------------------------


# Add no plants as a point (0%)

# Clear the decks
rm(list=ls())

# Bring in the data
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Time != 5) |>
	filter(!Species %in% c('Goutweed', 'Crepe Myrtle', 'Spider Lily'))

## --------------- PREPARE ALL DATA --------------------------------------------

# We are only interested in the first half of the transect data
feeder.veg <- feeder.veg |>
  filter(Meter %in% 1:10)

# Make Native
for(i in 1:nrow(feeder.veg)){
	if(is.na(feeder.veg$Functional[i])){
		feeder.veg$Functional[i] <- 'Native'
	}
}

# Make a predictor for individual transects
feeder.veg$Transect.ID <- paste(feeder.veg$Plot, feeder.veg$Feeder, 
 																feeder.veg$Transect)

feeder.veg$Time <- as_factor(feeder.veg$Time)

# Simple model
library(nnet) 
model <- multinom(Functional ~ Feeder
                  * Meter
									* Time
									+ Plot,
									  data = feeder.veg) 
summary(model)
tidy(model, conf.int = TRUE)

library(gtsummary)
tbl_regression(model, exp = TRUE)

library(ggeffects)
ggeffects::ggeffect(model, terms = c("Meter[0:10, by = 1]", "Time")) %>%
    plot()

Anova(model)

## --------------- VISUALIZE ANALYSES ------------------------------------------

# Use model to predict probabilities
pred.df <- feeder.veg |>
	dplyr::select(Plot, Transect.ID, Feeder, Meter, Time) |>
	unique()

resp <- predict(model, type="probs", newdata=pred.df)

# predictions <- ggeffects::ggemmeans(model, terms = c('Feeder', 'Meter [all]'),
# 																		nuisance = c('Plot'))

pred.df <- cbind(pred.df, resp)
pred.df$`Browsed perennials` <- round(pred.df$`Browsed perennials`, digits = 2)
pred.df$Native <- round(pred.df$Native, digits = 2)
pred.df$Ruderal <- round(pred.df$Ruderal, digits = 2)

pred.df <- pred.df |>
	pivot_longer(cols = 6:8, names_to = 'Functional', values_to = 'Probability')

pred.df <- pred.df |>
    mutate(Functional = factor(Functional)) |>
    mutate(Functional = fct_relevel(Functional, c("Native", "Browsed perennials", 
    																							"Ruderal"))) |>
	  mutate(Feeder = factor(Feeder)) |>
    mutate(Feeder = fct_relevel(Feeder, c("Y", "N"))) |>
    mutate(Time = factor(Time)) |>
    mutate(Time = fct_relevel(Time, c("0", "1", "2")))

class(pred.df$Time)

# Visualize probability facet wrapped
ggplot(pred.df, aes(x = Meter, y = Probability,
										color = Time))+
	scale_color_brewer(palette = "Dark2",
                       name = "")+
	scale_x_continuous(breaks = seq(0, 10, 1))+
	geom_jitter(shape = 21, alpha = 0.2)+
	geom_smooth(method = 'lm', formula = y~x, se = TRUE)+
	facet_wrap(~Feeder*Functional, ncol = 3, scales = "free_y")+
	theme_bw()+
	xlab('Distance from feeder (m)')+
	ylab('Detection probability')+
	theme(panel.grid.minor = element_blank(),
				panel.grid.major = element_blank(),
				legend.position = 'right',
				#strip.background = element_blank(),
  			#strip.text.x = element_blank(),
				aspect.ratio = 2,
				text=element_text(size=20),
				axis.title = element_text(face="bold"))
ggsave(filename = 'Figures/Manipulative-veg.png')


# Visualize probability together
# Save the predicted probabilities
pprob_metet_feeder <- ggeffect(model, terms = c("Meter[0:10, by = 1]", "Feeder",
																								"Plot", "Transect.ID"))
# Build the plot with ggplot manually
# The aesthetics refer to variable names in the ggeffects saved results
# We are grouping the plot by the "response.level" of our outcome.
  ggplot(data = pprob_metet_feeder, aes(x = x, y = predicted,
               color = response.level, group = response.level)) +
    # Add line layer
    geom_line() +
    # Add point layer
    geom_point() +
    # Add confidence intervals as error bars
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high,
                    color = response.level,
                    group = response.level),
                  width = 1) +
    # Change the color of the lines, remove legend title, change the
    # labels of the lines in the legend.
    scale_color_brewer(palette = "Dark2",
                       name = "") +
    # Add titles to the x and y axis
    labs(
      x = "Group as Percentage of Country's Population",
      y = "Probability"
    ) +
    # Set the theme
    theme_minimal() +
    theme(
      legend.position = "bottom", # move legend to the bottom
      axis.title = element_text(size = 14) # increase axis title size
      )+
  	facet_wrap(~group)

## --------------- CALCULATE EFFECT SIZE ---------------------------------------
  
pred.df <- pred.df %>% separate(Transect.ID, c('Site', 'Treatment', 'Transect'))|>
  	dplyr::select(-Site, -Treatment)
    
feeder.0 <- pred.df |> filter(Feeder == 'Y' & Time == 0)
colnames(feeder.0)[7] <- 'Feeder.prob.0'
feeder.0 <- feeder.0 |> dplyr::select(-Feeder, -Time)

feeder.1 <- pred.df |> filter(Feeder == 'Y' & Time == 1)
colnames(feeder.1)[7] <- 'Feeder.prob.1'
feeder.1 <- feeder.1 |> dplyr::select(-Feeder, -Time)

feeder.2 <- pred.df |> filter(Feeder == 'Y' & Time == 2)
colnames(feeder.2)[7] <- 'Feeder.prob.2'
feeder.2 <- feeder.2 |> dplyr::select(-Feeder, -Time)

control.0 <- pred.df |> filter(Feeder == 'N' & Time == 0)
colnames(control.0)[7] <- 'Control.prob.0'
control.0 <- control.0 |> dplyr::select(-Feeder, -Time)

control.1 <- pred.df |> filter(Feeder == 'N' & Time == 1)
colnames(control.1)[7] <- 'Control.prob.1'
control.1 <- control.1 |> dplyr::select(-Feeder, -Time)

control.2 <- pred.df |> filter(Feeder == 'N' & Time == 2)
colnames(control.2)[7] <- 'Control.prob.2'
control.2 <- control.2 |> dplyr::select(-Feeder, -Time)

es.df <- merge(feeder.0, feeder.1, all = FALSE)
es.df <- merge(es.df, feeder.2, all = FALSE)
es.df <- merge(es.df, control.0, all = FALSE)
es.df <- merge(es.df, control.1, all = FALSE)
es.df <- merge(es.df, control.2, all = FALSE)

# Calculate year to year change
es.df <- es.df |>
	mutate(Feeder.delta.1 = Feeder.prob.1-Feeder.prob.0,
				 Feeder.delta.2 = Feeder.prob.2-Feeder.prob.0,
				 Control.delta.1 = Control.prob.1-Control.prob.0,
				 Control.delta.2 = Control.prob.2-Control.prob.0) |>
	dplyr::select(Plot, Transect, Meter, Functional, Feeder.delta.1, Feeder.delta.2,
								Control.delta.1, Control.delta.2)

# Pivot longer
es.df <- es.df |>
	pivot_longer(cols = 5:8, names_to = 'Comparison', values_to = 'Delta')

# Create long dataframe
es.df$Feeder <- NA
for(i in 1:nrow(es.df)){
	if(isTRUE(es.df$Comparison[i] == 'Feeder.delta.1' | 
						es.df$Comparison[i] =='Feeder.delta.2')){
	  es.df$Feeder[i] <- 'Y'
	} else {
		es.df$Feeder[i] <- 'N'
	}
}

es.df$Time <- NA
for(i in 1:nrow(es.df)){
	if(isTRUE(es.df$Comparison[i] == 'Feeder.delta.1' | 
						es.df$Comparison[i] =='Control.delta.1')){
	  es.df$Time[i] <- '1'
	} else {
	  es.df$Time[i] <- '2'
	}
}

# Select variables
es.df <- es.df |>
	dplyr::select(Plot, Transect, Feeder, Time, Meter, Functional, Delta)

# Drop zeros
es.df.no.zero <- es.df |>
	filter(Delta != 0)

## --------------- MODEL AND VISUALIZE EFFECT SIZE -----------------------------

ggplot(es.df.no.zero, aes(x = Meter, y = Delta, linetype = Feeder,
													color = Functional))+
	geom_hline(yintercept = 0, color = 'black')+
	scale_color_brewer(palette = "Dark2",
                       name = "")+
	scale_x_continuous(breaks = seq(0, 10, 1))+
	scale_y_continuous(breaks = seq(-0.07, 0.07, 0.02))+
	geom_smooth()+
	facet_wrap(~Functional*Time, ncol = 2)+
		theme_bw()+
	xlab('Distance from feeder (m)')+
	ylab('Change in detection probability')+
	theme(panel.grid.minor = element_blank(),
				panel.grid.major = element_blank(),
				legend.position = 'right',
				strip.background = element_blank(),
  			strip.text.x = element_blank(),
				aspect.ratio = 0.9,
				text=element_text(size=20),
				axis.title = element_text(face="bold"))

ggsave('Figures/Manipulative-veg-delta.png', height = 10, width = 8,
			 unit = 'in')


## --------------- TIME 0 & 5 --------------------------------------------------

# Clear the decks
rm(list=ls())

# Bring in the data
baseline <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 0) |>
	filter(Veg.Type != 'na') |>
	dplyr::select(Plot, Time, Direction, Meter, Functional) # Transect

return <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 5) |>
	filter(Veg.Type != 'na') |>
	dplyr::select(Plot, Time, Direction, Meter, Functional)

# return <- left_join(return, unique(baseline))
# 
# return <- return |>
# 	dplyr::select(Plot, Direction, Transect, Meter, Functional)

feeder.veg <- rbind(baseline, return)

# We are only interested in the first half of the transect data
feeder.veg <- feeder.veg |>
  filter(Meter %in% 1:10)

# Make Native
for(i in 1:nrow(feeder.veg)){
	if(is.na(feeder.veg$Functional[i])){
		feeder.veg$Functional[i] <- 'Native'
	}
}

rm(baseline, return) # Clean up

unique(feeder.veg$Functional)

feeder.veg$Transect.ID <- paste(feeder.veg$Plot, feeder.veg$Feeder, 
 																feeder.veg$Direction)

feeder.veg$Time <- as_factor(feeder.veg$Time)
class(feeder.veg$Time)

## --------------- CONDUCT ANALYSES --------------------------------------------

# Simple model
library(nnet) 
model <- multinom(Functional ~ Time
                  * Meter, data = feeder.veg) 
summary(model)

library(broom)
tidy(model, conf.int = TRUE)

library(gtsummary)
tbl_regression(model, exp = TRUE)

library(ggeffects)
plot <- ggeffects::ggeffect(model, terms = c("Meter[0:10, by = 1]", "Time")) %>%
    plot()
plot

# Full model
model <- multinom(Functional ~ Time
                  * Meter
									+ Plot
									+ Transect.ID,
									data = feeder.veg) 

## --------------- VISUALIZE ANALYSES ------------------------------------------

# Use model to predict probabilities
pred.df <- feeder.veg |>
	dplyr::select(Plot, Transect.ID, Meter, Time)

resp <- predict(model, type="probs", newdata=pred.df)

# predictions <- ggeffects::ggemmeans(model, terms = c('Feeder', 'Meter [all]'),
# 																		nuisance = c('Plot'))

pred.df <- cbind(pred.df, resp)
pred.df$`Browsed perennials` <- round(pred.df$`Browsed perennials`, digits = 2)
pred.df$Native <- round(pred.df$Native, digits = 2)
pred.df$Ruderal <- round(pred.df$Ruderal, digits = 2)

pred.df <- pred.df |>
	pivot_longer(cols = 5:7, names_to = 'Functional', values_to = 'Probability')

pred.df <- pred.df |>
    mutate(Functional = factor(Functional)) |>
    mutate(Functional = fct_relevel(Functional, c("Native", "Browsed perennials", 
    																							"Ruderal"))) |>
	  mutate(Feeder = factor(Time)) |>
    mutate(Time = fct_relevel(Time, c("0", "5")))

class(pred.df$Time)

# Visualize probability facet wrapped
ggplot(pred.df, aes(x = Meter, y = Probability, fill = Functional,
										color = Functional))+
	geom_jitter(shape = 21, alpha = 0.2)+
	scale_color_brewer(palette = "Dark2",
                       name = "")+
	scale_fill_brewer(palette = "Dark2",
                       name = "")+
	scale_x_continuous(breaks = seq(0, 10, 1))+
	geom_smooth(method = 'lm', formula = y~x, size = 2)+
	facet_wrap(~Feeder*Functional, ncol = 3)+
	theme_bw()+
	xlab('Distance from feeder (m)')+
	ylab('Detection probability')+
	theme(panel.grid.minor = element_blank(),
				panel.grid.major = element_blank(),
				legend.position = 'none',
				strip.background = element_blank(),
  			strip.text.x = element_blank(),
				aspect.ratio = 2,
				text=element_text(size=20),
				axis.title = element_text(face="bold"))
ggsave(filename = 'Figures/Manipulative-before-after-veg.png')
