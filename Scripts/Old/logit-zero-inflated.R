
## --------------- SET—UP WORKSPACE --------------------------------------------

# Clear the decks
rm(list=ls())

# Bring in the data
feeder.veg <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative')

## --------------- PREPARE CLOSE DATA ------------------------------------------

# Subset the data
close <- feeder.veg |>
	filter(Meter < 6) # 5

for(i in 1:nrow(close)){
	if(is.na(close$Functional[i])){
		close$Functional[i] <- 'Native'
	}
}

# Tally natives and invasives at each plot*time
sites <- close |>
	group_by(Site, Plot, Time, Feeder, Transect) |>
	summarize(Everything = n())

# Tally invasive functional groups
func.invasive <- close |>
	group_by(Site, Plot, Time, Feeder, Functional, Transect) |>
	summarize(Count = n())

func.invasive <- merge(sites, func.invasive, all = TRUE)
func.invasive[is.na(func.invasive)] <- 0
func.invasive <- func.invasive |> dplyr::select(-Everything)

# Add zeroes
func.invasive <- func.invasive |>
	pivot_wider(names_from = Functional, values_from = Count)
func.invasive[is.na(func.invasive)] <- 0

# Calculate the proportion of invasive
# Calculate the proportion of invasive
func.invasive <- func.invasive |>
	mutate(Total.invasive = `Resistant perennials`+
				 	`Browsed perennials`+Ruderal,
				 Total.all = `Resistant perennials`+
				 	`Browsed perennials`+Ruderal + Native,
				 `Resistant perennials.tot` = `Resistant perennials`/Total.all,
				 `Browsed perennials.tot` = `Browsed perennials`/Total.all,
				 Ruderal.tot = Ruderal/Total.all,
				 Native.tot = Native/Total.all,
				 `Resistant perennials.inv` = `Resistant perennials`/Total.invasive,
				 `Browsed perennials.inv` = `Browsed perennials`/Total.invasive,
				 Ruderal.inv = Ruderal/Total.invasive) |>
	dplyr::select(-Total.all, -Total.invasive)
func.invasive[is.na(func.invasive)] <- 0

func.invasive <- func.invasive |>
	pivot_longer(6:16, names_to = 'Functional', values_to = 'Proportion')

## --------------- CLOSE INVASIVE MODEL ----------------------------------------

library(betareg)
library(zoib)
library(rjags)
library(lme4)
library(performance)
library(DHARMa)
library(car) # logit transformation
library(emmeans)

browsed.inv <- func.invasive |> 
	filter(Functional == "Browsed perennials.inv")

# Basic linear models won't cut it
browsed.inv.lm <- lmer(Proportion ~ Feeder * Time + (1|Site/Plot),
											 data = browsed.inv)
check_model(browsed.inv.lm)
sim.browsed.inv.lm <- simulateResiduals(browsed.inv.lm)
plot(sim.browsed.inv.lm)

# Logit transformation won't cut it
browsed.inv.logit <- lmer(logit(Proportion) ~ Feeder * Time + (1|Site/Plot),
											 data = browsed.inv)
check_model(browsed.inv.logit)
sim.browsed.inv.lm <- simulateResiduals(browsed.inv.lm)
plot(sim.browsed.inv.lm)

# try zero-inflated
browsed.inv <- browsed.inv |>
	mutate(Proportion.logit = logit(Proportion))
browsed.inv.zoib <- zoib(Proportion ~ data = browsed.inv)
zoib(Proportion ~ Feeder*Time|1|Feeder*Time|Feeder*Time|1, random = 1,    
       EUID = browsed.inv$Plot,  zero.inflation = TRUE,  
       one.inflation = TRUE, data = browsed.inv,   n.iter=50, 
       n.thin=20, n.burn=50)
# Arguments must be mcmc objects 

# There is no need to explore proportions because there are so few values
binary.invasive <- func.invasive |> 
	filter(Functional %in% c("Browsed perennials",
													 "Resistant perennials",
													 "Ruderal",
													 "Native"))

# for(i in 1:nrow(binary.invasive)){
	# if(isTRUE(binary.invasive$Binary[i] > 0)){
		# binary.invasive$Binary[i] <- 1
	# }
# }

rm(list=ls()[! ls() %in% c("binary.invasive", "close", "func.invasive", "feeder.veg")])


# Binomial
browsed <- binary.invasive |> 
	filter(Functional == "Browsed perennials")
browsed.mod <- glmer(Binary ~ Feeder * Time + (1|Site/Plot),
												 data = browsed, family = 'binomial')
Anova(browsed.inv.mod) # type 2 or 3?
emmeans(browsed.mod, pairwise~Feeder, type = 'response')


resistant <- binary.invasive |> 
	filter(Functional == "Resistant perennials")
resistant.mod <- glmer(Binary ~ Feeder * Time + (1|Site/Plot),
												 data = resistant, family = 'binomial')
Anova(resistant.inv.mod)
emmeans(resistant.mod, pairwise~Feeder, type = 'response')

ruderal <- binary.invasive |> 
	filter(Functional == "Ruderal")
ruderal.mod <- glmer(Binary ~ Feeder * Time + (1|Site/Plot),
												 data = ruderal, family = 'binomial')
Anova(ruderal.mod)
# Need to combine ruderal and resistant

# beta reg

browsed <- binary.invasive |> 
	filter(Functional == "Browsed perennials")

Anova(browsed.inv.mod) # type 2 or 3?
emmeans(browsed.mod, pairwise~Feeder, type = 'response')


resistant <- binary.invasive |> 
	filter(Functional == "Resistant perennials")
resistant.mod <- glmer(Binary ~ Feeder * Time + (1|Site/Plot),
												 data = resistant, family = 'binomial')
Anova(resistant.inv.mod)
emmeans(resistant.mod, pairwise~Feeder, type = 'response')

ruderal <- binary.invasive |> 
	filter(Functional == "Ruderal")
ruderal.mod <- glmer(Binary ~ Feeder * Time + (1|Site/Plot),
												 data = ruderal, family = 'binomial')
Anova(ruderal.mod)



## betareg
rm(list=ls()[! ls() %in% c("func.invasive")])

func.invasive <- func.invasive |>
	mutate(Proportion.trans = (Proportion * (794-1) + 0.5)/794)

browsed <- func.invasive |>
	filter(Functional == 'Browsed perennials.tot')
browsed.mod <- betareg(Proportion.trans ~ Feeder*Time + Plot,
											 data = browsed, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
Anova(browsed.mod)

resisant <- func.invasive |>
	filter(Functional == 'Resistant perennials.tot')
ruderal <- func.invasive |>
	filter(Functional == 'Ruderal.tot')





## Calculate all proportions
# Calculate the proportion of invasive
func.invasive <- func.invasive |>
	mutate(Total.invasive = `Resistant perennials`+
				 	`Browsed perennials`+Ruderal,
				 Total.all = `Resistant perennials`+
				 	`Browsed perennials`+Ruderal + Native,
				 `Resistant perennials.tot` = `Resistant perennials`/Total.all,
				 `Browsed perennials.tot` = `Browsed perennials`/Total.all,
				 Ruderal.tot = Ruderal/Total.all,
				 Native.tot = Native/Total.all,
				 `Resistant perennials.inv` = `Resistant perennials`/Total.invasive,
				 `Browsed perennials.inv` = `Browsed perennials`/Total.invasive,
				 Ruderal.inv = Ruderal/Total.invasive) |>
	dplyr::select(-Total.all, Total.invasive)
func.invasive[is.na(func.invasive)] <- 0


func.invasive <- func.invasive |>
	pivot_longer(6:17, names_to = 'Functional', values_to = 'Proportion')

# Not everything
# Calculate the proportion of invasive
func.invasive <- func.invasive |>
	mutate(Total = `Resistant perennials`+
				 	`Browsed perennials`+ Ruderal + Native,
				 `Resistant perennials` = `Resistant perennials`/Total,
				 `Browsed perennials` = `Browsed perennials`/Total,
				 Ruderal = Ruderal/Total,
				 Native = Native/Total) |>
	dplyr::select(-Total)

