# Create a df with all site/plot/time/feeder combinations
sites <- close |>
	group_by(Site, Plot, Time, Feeder) |>
	summarize(Everything = n())

# Tally invasive functional groups
func.invasive <- close |>
	group_by(Site, Plot, Time, Feeder, Functional) |>
	summarize(Count = n())

func.invasive <- merge(sites, func.invasive, all = TRUE)
func.invasive[is.na(func.invasive)] <- 0
func.invasive <- func.invasive |> dplyr::select(-Everything)

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

# Calculate the proportion of invasive species that each functional role
# accounts for
func.invasive <- func.invasive |>
	mutate(Total = `Resistant perennials`+
				 	`Browsed perennials`+Ruderal,
				 `Resistant perennials` = `Resistant perennials`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total,
				 Ruderal = Ruderal/Total,
				 Total.inv = (`Resistant perennials`+
				 	`Browsed perennials`+Ruderal)/(Total+Native)) |>
	dplyr::select(-Native, -Total)
func.invasive[is.na(func.invasive)] <- 0 

# Convert to longer
func.invasive <- func.invasive |>
	pivot_longer(5:8, names_to = 'Functional', values_to = 'Proportion')

# Convert to binary
for(i in 1:nrow(func.invasive)){
	if(func.invasive$Proportion[i] > 0){
		func.invasive$Proportion[i] <- 1
	}
}

total <- func.invasive |>
	filter(Functional == 'Total.inv')
total.mod <- glmer(Proportion ~ Feeder*as_factor(Time) + (1|Site),
									 data = total, family = 'binomial')
Anova(total.mod)
check_model(total.mod)
emmeans(total.mod, ~Feeder*Time)

browsed <- func.invasive |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- glmer(Proportion ~ Feeder*as_factor(Time) + (1|Site),
									 data = browsed, family = 'binomial')
Anova(browsed.mod)
check_model(browsed.mod)
emmeans(browsed.mod, ~Feeder*Time)

ruderals <- func.invasive |>
	filter(Functional == 'Ruderal')
ruderals.mod <- glmer(Proportion ~ Feeder*as_factor(Time) + (1|Site),
									 data = ruderals, family = 'binomial')
Anova(ruderals.mod)
check_model(ruderals.mod)
emmeans(ruderals.mod, ~Feeder*Time)

resistant <- func.invasive |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- glmer(Proportion ~ Feeder*as_factor(Time) + (1|Site),
									 data = resistant, family = 'binomial')
Anova(resistant.mod)
check_model(resistant.mod)
emmeans(resistant.mod, ~Feeder*Time)