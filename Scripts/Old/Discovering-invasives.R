# Clear the decks
rm(list=ls())

# Bring in the manipulative experiment data
d <- read.csv('Clean-data/Feeder-veg.csv') |> 
	filter(Dataset == 'Manipulative')

species.ls <- tibble(unique(d$Species))

ruderals <- c('false yam', 'foxtail', 'Barnyard Grass', 
				'Ivy Leaf Morning Glory', 'Millet', 
					'Bugleweed', 'purselane')
