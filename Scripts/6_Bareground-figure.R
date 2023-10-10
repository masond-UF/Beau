## --------------- HEADER ------------------------------------------------------
## Script name: 6_Bareground-figure.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-22
## Date Last modified: 2023-010-06
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: Creates a figure for the bareground data for the 
## feeder project.

rm(list=ls())
means <- read.csv('Results/Bare-ground/feeder-means.csv')

ggplot(means, aes(x = Feeder, y = Bare.ground))+
	geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0,
								size = 1)+
	geom_point(size = 7, shape = 21, aes(fill = Feeder),
						 color = 'black', stroke = 1.5)+
	scale_fill_manual(values = c('black', 'white'))+
	scale_y_continuous(limits = c(0, 0.4))+
	ylab('Mean bare ground probability')+
	xlab('')+
	theme_bw()+
	theme(text = element_text(size = 30),
				axis.text.x = element_text(face="bold"),
				legend.position = 'none',
				panel.grid.minor.y = element_blank(),
				panel.grid.major.x = element_blank(),
				panel.grid.major.y = element_line(color = 'black',
																					size = 0.1))+
	theme(aspect.ratio = 1.5)+
	annotate('text', x = 0.75, y = 0.4, label = 'p <0.001', size = 8)+
	theme(plot.margin = margin(0, 0, 0, 0, "pt"),
				axis.title = element_text(face = 'bold'))+
	theme(axis.title.y = element_text(margin = margin(t = 0, r = 1, b = 0, l = 3)))

ggsave(filename = "Figures/Bare-ground.png",
			 height = 10.5, width = 15, dpi = 400, units = 'in')

# Figure with time by year
rm(list=ls())
time <- read.csv('Results/Bare-ground/time-means.csv')

ggplot(time, aes(x = Feeder, y = Bare.ground))+
	geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0,
								size = 1)+
	geom_point(size = 7, shape = 21, aes(fill = Feeder),
						 color = 'black', stroke = 1.5)+
	scale_fill_manual(values = c('black', 'white'))+
	ylab('Mean bare ground probability')+
	xlab('')+
	theme_bw()+
	theme(text = element_text(size = 30),
				axis.text.x = element_text(face="bold"),
				legend.position = 'none',
				panel.grid.minor.y = element_blank(),
				panel.grid.major.x = element_blank(),
				panel.grid.major.y = element_line(color = 'black',
																					size = 0.1))+
	theme(aspect.ratio = 1.5)+
	theme(plot.margin = margin(0, 0, 0, 0, "pt"),
				axis.title = element_text(face = 'bold'))+
	theme(axis.title.y = element_text(margin = margin(t = 0, r = 1, b = 0, l = 3)))+
	facet_wrap(~Time)

ggsave(filename = "Figures/Bare-ground-time.png",
			 height = 10.5, width = 15, dpi = 400, units = 'in')
