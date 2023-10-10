## --------------- HEADER ------------------------------------------------------
## Script name: 4_Herbivory-figure.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-06-22
## Date Last modified: 2023-10-06
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: Creates a figure for the herbivory data for the 
## feeder project.

rm(list=ls())

means <- read.csv('Results/Herbivory/feeder-means.csv')
  
# Contrast  Y / N
# Odds ratio 2.38 
# SE 0.262 
# df Inf    
# null 1   
# z.ratio 7.859
# p <.0001

ggplot(means, aes(x = Feeder, y = Herbivory))+
	geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0,
								size = 1)+
	geom_point(size = 7, shape = 21, aes(fill = Feeder),
						 color = 'black', stroke = 1.5)+
	scale_fill_manual(values = c('black', 'white'))+
	scale_y_continuous(limits = c(0, 0.02))+
	ylab('Mean herbivory probability')+
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
	annotate('text', x = 0.75, y = 0.02, label = 'p <0.001', size = 8)+
	theme(plot.margin = margin(0, 0, 0, 0, "pt"),
				axis.title = element_text(face = 'bold'))+
	theme(axis.title.y = element_text(margin = margin(t = 0, r = 1, b = 0, l = 3)))

ggsave(filename = "Figures/Herbivory.png",
			 height = 10.5, width = 15, dpi = 400, units = 'in')
