## --------------- HEADER ------------------------------------------------------
## Script name: 7_Plant-detection-retrospective.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-22
## Date Last modified: 2023-06-22
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This is a script for analyzing the retrospective plant
## feeder data.

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

# Disable exponents
options(scipen = 999)

# Clear the decks
rm(list = ls())

# Bring in the data
feeder.veg <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Retrospective")

## --------------- PREPARE DATA ------------------------------------------------

# We are only interested in the first half of the transect data
feeder.veg <- feeder.veg |>
  filter(Meter %in% 1:10)

# Drop non plants
feeder.veg <- feeder.veg |>
  filter(Veg.Type != 'na')

# Make Native
for(i in 1:nrow(feeder.veg)){
	if(is.na(feeder.veg$Functional[i])){
		feeder.veg$Functional[i] <- 'Native'
	}
}

# Make a predictor for individual transects
feeder.veg$Transect.ID <- paste(feeder.veg$Plot, feeder.veg$Feeder, 
 																feeder.veg$Transect)
class(feeder.veg$Meter)

## --------------- CONDUCT ANALYSES --------------------------------------------

# Simple model
library(nnet) 
model <- multinom(Functional ~ Feeder
                  * Meter, data = feeder.veg) 
summary(model)
tidy(model, conf.int = TRUE)

library(gtsummary)
tbl_regression(model, exp = TRUE)

library(ggeffects)
ggeffects::ggeffect(model, terms = c("Meter[0:10, by = 1]", "Feeder")) %>%
    plot()

# Full model
model <- multinom(Functional ~ Feeder
                  * Meter
									+ Plot
									+ Transect.ID,
									data = feeder.veg) 
# ggeffects::ggeffect(model, terms = c("Meter[0:10, by = 1]", "Feeder",
#																		 "Plot", "Transect.ID")) %>%
#    plot()

# Check assumptions (carries over from logistic regression)

# Statistical tests of 
Anova(model)
# Meter and Feeder*Meter significant



linearHypothesis(model,
                 c("Native:Meter = 0",
                   "Browsed perennials:Meter = 0",
                 	 "Ruderal:Meter = 0"))

## --------------- VISUALIZE ANALYSES ------------------------------------------

# Use model to predict probabilities
pred.df <- feeder.veg |>
	dplyr::select(Plot, Transect.ID, Feeder, Meter)

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
	  mutate(Feeder = factor(Feeder)) |>
    mutate(Feeder = fct_relevel(Feeder, c("Y", "N")))


# Visualize probability facet wrapped
ggplot(pred.df, aes(x = Meter, y = Probability, fill = Functional,
										color = Functional))+
	geom_jitter(shape = 21)+
	scale_color_brewer(palette = "Dark2",
                       name = "")+
	scale_fill_brewer(palette = "Dark2",
                       name = "")+
	scale_x_continuous(breaks = seq(0, 10, 1))+
	geom_smooth(method = 'lm', formula = y~x)+
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
ggsave(filename = 'Figures/Retrospective-veg.png')


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
