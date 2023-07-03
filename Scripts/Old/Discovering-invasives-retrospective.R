
# Clear the decks
rm(list = ls())

# Bring in the data
feeder.veg <- read.csv("Clean-data/Feeder-veg.csv") |>
  filter(Dataset == "Retrospective")

spec.ls <- feeder.veg |>
	dplyr::select(Invasive, Species) |>
	unique() |>
	tibble()

check.ls <- c('brome', 'carpet weed')